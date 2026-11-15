"""Archive integrity validation, including deterministic replay when known.

Finite-family and geometry datasets are regenerated one run at a time. Custom
sources receive generic schema/count/reference validation rather than an
unsupported claim that their scientific outputs have been independently proved.
"""

from __future__ import annotations

from itertools import zip_longest
import os
from pathlib import Path
from typing import Any, Iterator

from .annotations import check_annotation, check_note, iter_expected_divergence
from .core import (EnumerationConfig, GENERATOR_VERSION, SCHEMA_VERSION,
                   canonical_json, iter_records, require_int, stable_id)
from .storage import (ARCHIVE_FILES, DatasetError, check_run_record, count_outputs,
                      dataset_lock, load_manifest, read_jsonl)

MANIFEST_FIELDS = frozenset({
    "schema_version", "config_version", "generator_version", "dataset_id", "source", "config",
    "state", "completed", "expected_runs", "expected_outputs", "total_runs", "total_outputs",
    "error_records", "files", "failure",
})


def _canonical_rows(path: Path) -> Iterator[dict[str, Any]]:
    # Parsing is strict; also verify byte-level canonical encoding and newlines.
    with path.open("r", encoding="utf-8", newline="") as raw:
        for line_number, row in enumerate(read_jsonl(path), 1):
            line = raw.readline()
            if line != canonical_json(row) + "\n":
                raise DatasetError(f"{path.name}:{line_number}: JSONL is not canonical")
            yield row


def _validate_unlocked(path: Path, report: dict[str, Any], allow_incomplete: bool) -> None:
    errors = report["errors"]
    manifest = load_manifest(path)
    report["completed"] = manifest.get("completed") is True
    if set(manifest) != MANIFEST_FIELDS:
        errors.append("manifest has missing or undeclared fields")
    if manifest.get("config_version") != SCHEMA_VERSION:
        errors.append("unsupported or missing config_version")
    if manifest.get("files") != list(ARCHIVE_FILES):
        errors.append("manifest archive-file list differs from the schema")
    state = manifest.get("state")
    if state not in ("started", "complete", "failed"):
        errors.append("manifest state must be started, complete, or failed")
    if type(manifest.get("completed")) is not bool:
        errors.append("manifest completed must be a boolean")
    if manifest.get("completed") is not (state == "complete"):
        errors.append("manifest state/completed mismatch")
    if state != "complete" and not allow_incomplete:
        errors.append(f"dataset is not complete (manifest state: {state})")
    for name in ("expected_runs", "expected_outputs", "total_runs", "total_outputs", "error_records"):
        require_int(f"manifest {name}", manifest.get(name), 0)
    if state == "failed":
        failure = manifest.get("failure")
        if (not isinstance(failure, dict) or set(failure) != {"type", "message"}
                or not isinstance(failure["type"], str) or not isinstance(failure["message"], str)):
            errors.append("failed manifest must describe the actual failure")
    elif manifest.get("failure") is not None:
        errors.append("nonfailed manifest failure must be null")
    if state == "complete" and manifest["error_records"] != 0:
        errors.append("complete dataset may not contain error records")
    if not isinstance(manifest.get("config"), dict):
        raise DatasetError("manifest config must be an object")
    source, config = manifest.get("source"), manifest["config"]
    if not isinstance(source, str) or not source:
        raise DatasetError("manifest source must be nonempty text")
    expected_iterator: Iterator[dict[str, Any]] | None = None
    declared_counts = {"runs": manifest["expected_runs"], "outputs": manifest["expected_outputs"]}
    if source == "finite":
        finite_config = EnumerationConfig.from_dict(config)
        expected_iterator = iter_records(finite_config)
        expected_counts = finite_config.expected_counts()
        expected_dataset_id = finite_config.dataset_id
    elif source == "geometry":
        from .geometry import count_geometry, iter_geometry_runs, validate_geometry_config
        normalized = validate_geometry_config(config)
        if canonical_json(config) != canonical_json(normalized):
            errors.append("geometry config is not normalized")
        expected_iterator = iter_geometry_runs(normalized)
        expected_counts = count_geometry(normalized)
        expected_dataset_id = stable_id("dataset", {
            "schema_version": SCHEMA_VERSION, "generator_version": GENERATOR_VERSION,
            "source": source, "config": config})
    elif source == "projection":
        from .projections import count_projection, iter_projection_runs, validate_projection_config
        normalized = validate_projection_config(config)
        if canonical_json(config) != canonical_json(normalized):
            errors.append("projection config is not normalized")
        expected_iterator = iter_projection_runs(normalized)
        expected_counts = count_projection(normalized)
        expected_dataset_id = stable_id("dataset", {
            "schema_version": SCHEMA_VERSION, "generator_version": GENERATOR_VERSION,
            "source": source, "config": config})
    else:
        expected_counts = declared_counts
        expected_dataset_id = stable_id("dataset", {
            "schema_version": SCHEMA_VERSION, "generator_version": GENERATOR_VERSION,
            "source": source, "config": config})
    if manifest.get("dataset_id") != expected_dataset_id:
        errors.append("manifest dataset_id does not match its exact config and versions")
    if declared_counts != expected_counts:
        errors.append("manifest expected counts differ from configured finite bounds")
    report["coverage_check"] = "deterministic_replay" if expected_iterator is not None else "counts_only"
    known: set[str] = set()
    actual_runs = actual_outputs = error_count = 0
    saw_error = False
    for row in _canonical_rows(path / "records.jsonl"):
        check_run_record(row)
        if row["id"] in known:
            errors.append(f"duplicate archive ID: {row['id']}")
        known.add(row["id"])
        if row.get("kind") == "error":
            saw_error = True
            error_count += 1
            expected_error_id = stable_id("error", {"dataset_id": manifest["dataset_id"],
                                                    "after_runs": actual_runs})
            if (row["id"] != expected_error_id or row.get("run_id") != expected_error_id
                    or row["status"] != "OBSERVED" or row["outputs"] != []
                    or row["params"] != manifest.get("failure")):
                errors.append("error record does not match the failed manifest")
            if state != "failed":
                errors.append("actual error records require a failed manifest")
            continue
        if saw_error:
            errors.append("normal run follows an archive error record")
        actual_runs += 1
        actual_outputs += count_outputs(row)
        if expected_iterator is not None:
            expected = next(expected_iterator, None)
            if expected is None or canonical_json(row) != canonical_json(expected):
                if "run/state coverage differs from deterministic configuration replay" not in errors:
                    errors.append("run/state coverage differs from deterministic configuration replay")
    report["counts"].update(runs=actual_runs, outputs=actual_outputs, error_records=error_count)
    if actual_runs != manifest["total_runs"] or actual_outputs != manifest["total_outputs"]:
        errors.append("manifest actual run/output totals do not match the archive")
    if error_count != manifest["error_records"] or error_count > 1:
        errors.append("manifest error-record count does not match the archive")
    if actual_runs > manifest["expected_runs"] or actual_outputs > manifest["expected_outputs"]:
        errors.append("archive exceeds declared finite coverage")
    if state == "complete":
        if actual_runs != manifest["expected_runs"] or actual_outputs != manifest["expected_outputs"]:
            errors.append("complete archive is missing declared runs or outputs")
        if expected_iterator is not None and next(expected_iterator, None) is not None:
            errors.append("complete archive omits deterministically expected runs")
    annotations = []
    base_ids = set(known)
    for sequence, row in enumerate(_canonical_rows(path / "observers.jsonl"), 1):
        check_annotation(row, sequence)
        if row["entry_id"] not in base_ids:
            errors.append(f"annotation references unknown base entry: {row['entry_id']}")
        if row["id"] in known:
            errors.append(f"duplicate archive ID: {row['id']}")
        known.add(row["id"])
        annotations.append(row)
    report["counts"]["annotations"] = len(annotations)
    divergence_count = 0
    missing = object()
    for actual, expected in zip_longest(_canonical_rows(path / "divergence.jsonl"),
                                       iter_expected_divergence(annotations), fillvalue=missing):
        if actual is not missing:
            divergence_count += 1
            if not isinstance(actual.get("id"), str) or actual.get("id") in known:
                errors.append("divergence entry has a missing/duplicate ID")
            else:
                known.add(actual["id"])
        if (actual is missing or expected is missing
                or canonical_json(actual) != canonical_json(expected)):
            if "divergence log differs from replay of explicit annotations" not in errors:
                errors.append("divergence log differs from replay of explicit annotations")
    report["counts"]["divergence_events"] = divergence_count
    notes = 0
    for notes, row in enumerate(_canonical_rows(path / "notes.jsonl"), 1):
        check_note(row, notes)
        if row["id"] in known:
            errors.append(f"duplicate archive ID: {row['id']}")
        if any(reference not in known for reference in row["references"]):
            errors.append(f"note has unknown or forward references: {row['id']}")
        known.add(row["id"])
    report["counts"]["notes"] = notes


def validate_dataset(dataset: str | os.PathLike[str], *, allow_incomplete: bool = False) -> dict[str, Any]:
    """Return ``{valid, completed, errors, counts, coverage_check}``.

    By default, incomplete datasets are invalid even when their recorded prefix
    is sound. ``allow_incomplete=True`` checks prefix integrity without treating
    incompletion alone as corruption. Validation never repairs or edits logs.
    """
    if type(allow_incomplete) is not bool:
        raise ValueError("allow_incomplete must be boolean")
    report: dict[str, Any] = {"valid": False, "completed": False, "errors": [], "counts": {},
                              "coverage_check": "unavailable"}
    try:
        with dataset_lock(dataset) as path:
            _validate_unlocked(path, report, allow_incomplete)
    except (OSError, ValueError, TypeError, KeyError, ImportError, RecursionError, OverflowError) as exc:
        report["errors"].append(str(exc))
    report["valid"] = not report["errors"]
    return report
