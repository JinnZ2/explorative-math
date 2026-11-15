"""Durable append-only archives and cooperating-process dataset locks.

A lock is deliberately never reclaimed automatically: a leftover lock after a
hard process kill requires human inspection. JSONL appends are individually
flushed and fsynced. A run's complete outputs occupy one JSONL record.
"""

from __future__ import annotations

from contextlib import contextmanager
import errno
import json
import os
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping

from .core import (EnumerationConfig, GENERATOR_VERSION, OBSERVER_SLOTS,
                   SCHEMA_VERSION, STATUSES, ValidationError, canonical_json,
                   iter_records, require_int, stable_id)

ARCHIVE_FILES = ("records.jsonl", "observers.jsonl", "divergence.jsonl", "notes.jsonl")
REQUIRED_RUN_FIELDS = frozenset({"id", "params", "outputs", "status", *OBSERVER_SLOTS})
OPTIONAL_RUN_FIELDS = frozenset({"family", "run_id", "kind"})


class DatasetError(ValidationError):
    """The archive is unavailable, malformed, or locked by another writer."""


def _reject_constant(value: str) -> None:
    raise DatasetError(f"nonfinite JSON number is forbidden: {value}")


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise DatasetError(f"duplicate JSON object key: {key}")
        result[key] = value
    return result


def decode_json(text: str) -> Any:
    """Reject nonstandard constants and duplicate keys, not just syntax errors."""
    try:
        value = json.loads(text, parse_constant=_reject_constant,
                           object_pairs_hook=_unique_object)
        # JSON's 1e999 is legal syntax but becomes a nonfinite float in Python.
        canonical_json(value)
        return value
    except (ValueError, TypeError, RecursionError) as exc:
        raise DatasetError(f"invalid JSON: {exc}") from exc


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = decode_json(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError) as exc:
        raise DatasetError(f"cannot read {path.name}: {exc}") from exc
    if not isinstance(value, dict):
        raise DatasetError(f"{path.name} must contain an object")
    return value


def read_jsonl(path: Path) -> Iterator[dict[str, Any]]:
    """Stream complete JSONL objects; blank or unterminated lines are errors."""
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            for line_number, line in enumerate(stream, 1):
                if not line.endswith("\n"):
                    raise DatasetError(f"{path.name}:{line_number}: incomplete line")
                try:
                    value = decode_json(line)
                except DatasetError as exc:
                    raise DatasetError(f"{path.name}:{line_number}: {exc}") from exc
                if not isinstance(value, dict):
                    raise DatasetError(f"{path.name}:{line_number}: expected an object")
                yield value
    except (OSError, UnicodeError) as exc:
        raise DatasetError(f"cannot read {path.name}: {exc}") from exc


def fsync_directory(path: Path) -> None:
    """Persist directory entries where directory fsync is supported."""
    if os.name == "nt":
        return
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        try:
            os.fsync(descriptor)
        except OSError as exc:
            if exc.errno not in (errno.EINVAL, errno.ENOTSUP):
                raise
    finally:
        os.close(descriptor)


@contextmanager
def dataset_lock(dataset: str | os.PathLike[str]) -> Iterator[Path]:
    """Serialize operations using ``.write.lock`` created with O_EXCL.

    Cooperating readers also use this lock for consistent snapshots. There is
    no waiting, lock stealing, or stale-lock deletion; contention fails clearly.
    """
    path = Path(dataset)
    if not path.is_dir():
        raise DatasetError(f"dataset directory does not exist: {path}")
    lock = path / ".write.lock"
    try:
        descriptor = os.open(lock, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    except FileExistsError as exc:
        raise DatasetError(f"dataset is locked: {lock}; inspect stale locks before removal") from exc
    inode = os.fstat(descriptor).st_ino
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            stream.write(canonical_json({"pid": os.getpid()}) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        fsync_directory(path)
        yield path
    finally:
        try:
            if lock.stat().st_ino == inode:
                lock.unlink()
                fsync_directory(path)
        except FileNotFoundError:
            pass


def write_manifest(dataset: Path, manifest: Mapping[str, Any]) -> None:
    """Atomically replace and fsync the manifest while holding the lock."""
    temporary = dataset / ".manifest.tmp"
    try:
        with temporary.open("w", encoding="utf-8", newline="\n") as stream:
            stream.write(canonical_json(manifest) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, dataset / "manifest.json")
        fsync_directory(dataset)
    finally:
        if temporary.exists():
            temporary.unlink()


def append_jsonl(path: Path, record: Mapping[str, Any]) -> None:
    """Append one pre-encoded canonical line, flush, and fsync.

    Caller must hold the dataset lock. Serialization occurs before opening the
    file, so invalid JSON values cannot leave a partially serialized record.
    """
    encoded = (canonical_json(dict(record)) + "\n").encode("utf-8")
    with path.open("ab") as stream:
        stream.write(encoded)
        stream.flush()
        os.fsync(stream.fileno())


def count_outputs(record: Mapping[str, Any]) -> int:
    """Count list outputs individually, or one compound dictionary output."""
    outputs = record["outputs"]
    if isinstance(outputs, list):
        return len(outputs)
    if isinstance(outputs, dict):
        return 1
    raise DatasetError("run outputs must be a list or an object")


def _check_nested_statuses(value: Any) -> None:
    if isinstance(value, dict):
        if "status" in value and (not isinstance(value["status"], str) or value["status"] not in STATUSES):
            raise DatasetError("nested output status must be OBSERVED, DERIVED, or PROPOSED")
        for item in value.values():
            _check_nested_statuses(item)
    elif isinstance(value, list):
        for item in value:
            _check_nested_statuses(item)


def check_run_record(record: Any) -> dict[str, Any]:
    """Validate the generic per-run schema and null-only base observer slots."""
    if not isinstance(record, dict):
        raise DatasetError("each run must be a dictionary")
    if not REQUIRED_RUN_FIELDS <= set(record):
        raise DatasetError("run is missing required top-level fields")
    if set(record) - REQUIRED_RUN_FIELDS - OPTIONAL_RUN_FIELDS:
        raise DatasetError("run has undeclared top-level fields")
    if not isinstance(record["id"], str) or not record["id"]:
        raise DatasetError("run id must be a nonempty string")
    if not isinstance(record["status"], str) or record["status"] not in STATUSES:
        raise DatasetError("run status must be OBSERVED, DERIVED, or PROPOSED")
    if not isinstance(record["params"], dict):
        raise DatasetError("run params must be an object")
    if any(record[slot] is not None for slot in OBSERVER_SLOTS):
        raise DatasetError("base run observer slots must all be null")
    for name in OPTIONAL_RUN_FIELDS:
        if name in record and (not isinstance(record[name], str) or not record[name]):
            raise DatasetError(f"{name} must be a nonempty string")
    count_outputs(record)
    _check_nested_statuses(record["outputs"])
    canonical_json(record)
    return record


def load_manifest(dataset: Path) -> dict[str, Any]:
    manifest = read_json(dataset / "manifest.json")
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise DatasetError("unsupported or missing schema_version")
    if manifest.get("generator_version") != GENERATOR_VERSION:
        raise DatasetError("unsupported or missing generator_version")
    return manifest


def create_dataset(out: str | os.PathLike[str], config: Mapping[str, Any] | EnumerationConfig,
                   runs: Iterable[dict[str, Any]], expected_runs: int,
                   expected_outputs: int, *, source: str = "custom") -> dict[str, Any]:
    """Create a new archive from a streaming iterable of complete run records.

    Existing paths are refused, including empty directories. Memory use is one
    run plus a set of IDs for duplicate detection, not the whole output space.
    Every append is durable. Failed generation writes an actual OBSERVED error
    record when possible and a failed manifest with completed=false.
    """
    require_int("expected_runs", expected_runs, 0)
    require_int("expected_outputs", expected_outputs, 0)
    if not isinstance(source, str) or not source:
        raise ValidationError("source must be a nonempty string")
    if isinstance(config, EnumerationConfig):
        config_dict = config.to_dict()
    elif isinstance(config, Mapping):
        config_dict = dict(config)
    else:
        raise ValidationError("config must be an object or EnumerationConfig")
    # Round-trip to detach mutable caller containers and reject non-JSON values.
    config_dict = decode_json(canonical_json(config_dict))
    dataset_id = (EnumerationConfig.from_dict(config_dict).dataset_id if source == "finite"
                  else stable_id("dataset", {"schema_version": SCHEMA_VERSION,
                                             "generator_version": GENERATOR_VERSION,
                                             "source": source, "config": config_dict}))
    path = Path(out)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        path.mkdir()
    except FileExistsError as exc:
        raise DatasetError(f"refusing existing output dataset: {path}") from exc
    fsync_directory(path.parent)
    manifest: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION, "config_version": SCHEMA_VERSION,
        "generator_version": GENERATOR_VERSION, "dataset_id": dataset_id,
        "source": source, "config": config_dict, "state": "started", "completed": False,
        "expected_runs": expected_runs, "expected_outputs": expected_outputs,
        "total_runs": 0, "total_outputs": 0, "error_records": 0,
        "files": list(ARCHIVE_FILES), "failure": None,
    }
    with dataset_lock(path):
        write_manifest(path, manifest)
        ids: set[str] = set()
        try:
            for filename in ARCHIVE_FILES:
                with (path / filename).open("xb") as stream:
                    stream.flush()
                    os.fsync(stream.fileno())
            fsync_directory(path)
            for record in runs:
                check_run_record(record)
                if record["id"] in ids:
                    raise DatasetError(f"duplicate run id: {record['id']}")
                if record.get("kind") == "error":
                    raise DatasetError("generator may not yield archive error records")
                output_count = count_outputs(record)
                if manifest["total_runs"] >= expected_runs:
                    raise DatasetError("generator yielded more runs than declared")
                if manifest["total_outputs"] + output_count > expected_outputs:
                    raise DatasetError("generator yielded more outputs than declared")
                append_jsonl(path / "records.jsonl", record)
                ids.add(record["id"])
                manifest["total_runs"] += 1
                manifest["total_outputs"] += output_count
            if (manifest["total_runs"] != expected_runs or
                    manifest["total_outputs"] != expected_outputs):
                raise DatasetError("generator coverage differs from declared run/output counts")
            manifest.update(state="complete", completed=True)
            write_manifest(path, manifest)
        except BaseException as exc:
            manifest.update(state="failed", completed=False,
                            failure={"type": type(exc).__name__, "message": str(exc)})
            error_id = stable_id("error", {"dataset_id": dataset_id,
                                           "after_runs": manifest["total_runs"]})
            error_record = {"id": error_id, "run_id": error_id, "kind": "error",
                            "params": dict(manifest["failure"]), "outputs": [],
                            "status": "OBSERVED", **dict.fromkeys(OBSERVER_SLOTS)}
            try:
                append_jsonl(path / "records.jsonl", error_record)
                manifest["error_records"] = 1
            except (OSError, ValueError, TypeError):
                pass
            try:
                write_manifest(path, manifest)
            except OSError:
                pass
            raise
    return manifest


def enumerate_dataset(out: str | os.PathLike[str],
                      config: EnumerationConfig | None = None) -> dict[str, Any]:
    """Enumerate the finite CA, binary-substitution, and integer-map spaces."""
    config = EnumerationConfig() if config is None else config
    if not isinstance(config, EnumerationConfig):
        raise ValidationError("config must be an EnumerationConfig")
    counts = config.expected_counts()
    return create_dataset(out, config, iter_records(config), counts["runs"], counts["outputs"],
                          source="finite")


def find_record(dataset: Path, entry: str) -> dict[str, Any]:
    """Find one stable base ID without loading the whole archive."""
    if not isinstance(entry, str) or not entry:
        raise ValidationError("entry must be a nonempty ID")
    found: dict[str, Any] | None = None
    for record in read_jsonl(dataset / "records.jsonl"):
        if record.get("id") == entry:
            if found is not None:
                raise DatasetError(f"duplicate base id: {entry}")
            found = check_run_record(record)
    if found is None:
        raise DatasetError(f"entry not found: {entry}")
    return found
