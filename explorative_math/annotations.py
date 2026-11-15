"""Explicit annotations without automated observers or interest filtering.

Base records are immutable and remain null in all observer slots. The observed
fact in an annotation is that a named observer reported a value, not that the
claim is objectively true. Only the disagreement log is mechanically DERIVED.
"""

from __future__ import annotations

from datetime import datetime, timezone
import math
import os
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping

from .core import OBSERVER_SLOTS, STATUSES, ValidationError
from .storage import (DatasetError, append_jsonl, dataset_lock, decode_json,
                      find_record, load_manifest, read_jsonl)

INTEREST_SLOTS = ("human_interest", "ai_interest")
DIVERGENCE_RULE = "explicit_nonnull_boolean_disagreement:v1"
ANNOTATION_FIELDS = frozenset({"id", "entry_id", "slot", "value", "observer", "reason",
                              "method", "status", "timestamp"})
NOTE_FIELDS = frozenset({"id", "status", "text", "references", "timestamp"})


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds").replace("+00:00", "Z")


def check_timestamp(value: Any) -> None:
    if not isinstance(value, str) or not value.endswith("Z"):
        raise DatasetError("timestamp must be an ISO 8601 UTC string ending in Z")
    try:
        parsed = datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as exc:
        raise DatasetError("invalid UTC timestamp") from exc
    if parsed.utcoffset() != timezone.utc.utcoffset(parsed):
        raise DatasetError("timestamp must be UTC")


def nonblank(name: str, value: Any) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValidationError(f"{name} must be nonblank text")
    return value


def check_observer_value(slot: str, value: Any, method: str | None = None) -> None:
    """Strict JSON-compatible values; null always means unknown, never false."""
    if slot not in OBSERVER_SLOTS:
        raise ValidationError(f"slot must be one of {', '.join(OBSERVER_SLOTS)}")
    if method is not None:
        nonblank("method", method)
    if slot in INTEREST_SLOTS:
        if value is not None and type(value) is not bool:
            raise ValidationError("interest value must be true, false, or null")
    elif slot == "compressibility":
        if value is not None:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValidationError("compressibility must be null or a nonnegative finite number")
            if value < 0 or (isinstance(value, float) and not math.isfinite(value)):
                raise ValidationError("compressibility must be null or a nonnegative finite number")
            nonblank("method (required for numeric compressibility)", method)
    elif value is not None and not isinstance(value, str):
        raise ValidationError("divergence_note must be text or null")


def parse_observer_value(slot: str, text: str) -> Any:
    """Parse CLI values; divergence_note permits plain text or a JSON string."""
    if slot == "divergence_note":
        if text.strip() == "null":
            return None
        if text.lstrip().startswith('"'):
            value = decode_json(text)
            if not isinstance(value, str):
                raise ValidationError("divergence_note must be text or null")
            return value
        return text
    return decode_json(text)


def check_annotation(row: Mapping[str, Any], sequence: int | None = None) -> None:
    if set(row) != ANNOTATION_FIELDS:
        raise DatasetError("annotation has missing or undeclared fields")
    nonblank("annotation id", row["id"])
    nonblank("annotation entry_id", row["entry_id"])
    if sequence is not None and row["id"] != f"annotation-{sequence:012d}":
        raise DatasetError("annotation ID/order is not canonical")
    if row["status"] != "OBSERVED":
        raise DatasetError("observer reports must have status OBSERVED")
    nonblank("observer", row["observer"])
    nonblank("reason", row["reason"])
    check_timestamp(row["timestamp"])
    check_observer_value(row["slot"], row["value"], row["method"])


def is_divergent(values: Mapping[str, Any]) -> bool:
    human, ai = values["human_interest"], values["ai_interest"]
    return type(human) is bool and type(ai) is bool and human != ai


def derive_divergence(annotation: Mapping[str, Any], previous: Mapping[str, Any],
                      current: Mapping[str, Any], annotation_ids: Mapping[str, Any],
                      sequence: int) -> dict[str, Any] | None:
    """Derive an event only when effective disagreement or its values change."""
    before = {slot: previous[slot] for slot in INTEREST_SLOTS}
    after = {slot: current[slot] for slot in INTEREST_SLOTS}
    if before == after:
        return None
    was_open, now_open = is_divergent(before), is_divergent(after)
    if not was_open and not now_open:
        return None
    event = "opened" if not was_open else ("changed" if now_open else "resolved")
    return {"id": f"divergence-{sequence:012d}", "entry_id": annotation["entry_id"],
            "status": "DERIVED", "event": event, "values": after,
            "previous_values": before,
            "annotation_ids": {slot: annotation_ids[slot] for slot in INTEREST_SLOTS},
            "cause_event": "explicit_annotation", "cause_annotation_id": annotation["id"],
            "derivation_rule": DIVERGENCE_RULE, "timestamp": annotation["timestamp"]}


def iter_expected_divergence(annotations: Iterable[Mapping[str, Any]]) -> Iterator[dict[str, Any]]:
    """Replay explicit reports; no base-record values or automatic observers."""
    effective: dict[str, dict[str, Any]] = {}
    provenance: dict[str, dict[str, Any]] = {}
    sequence = 1
    for annotation in annotations:
        entry = annotation["entry_id"]
        current = effective.setdefault(entry, dict.fromkeys(OBSERVER_SLOTS))
        ids = provenance.setdefault(entry, dict.fromkeys(OBSERVER_SLOTS))
        previous = dict(current)
        current[annotation["slot"]] = annotation["value"]
        ids[annotation["slot"]] = annotation["id"]
        event = derive_divergence(annotation, previous, current, ids, sequence)
        if event is not None:
            yield event
            sequence += 1


def _annotation_snapshot(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    annotations = list(read_jsonl(path / "observers.jsonl"))
    for sequence, row in enumerate(annotations, 1):
        check_annotation(row, sequence)
    expected = list(iter_expected_divergence(annotations))
    actual = list(read_jsonl(path / "divergence.jsonl"))
    if actual != expected:
        raise DatasetError("divergence log does not match explicit annotations; run validate")
    return annotations, actual


def _require_modifiable(path: Path) -> None:
    manifest = load_manifest(path)
    if manifest.get("state") not in ("complete", "failed"):
        raise DatasetError("cannot annotate a dataset with an unfinished started manifest")


def observe(dataset: str | os.PathLike[str], entry: str, slot: str, value: Any,
            observer: str, reason: str, method: str | None = None) -> dict[str, Any]:
    """Durably append an explicit OBSERVED report and any DERIVED disagreement.

    Numeric compressibility requires a named method; no metric is computed.
    Duplicate reports remain in history but unchanged effective values never
    create redundant divergence events. A crash between the two fsynced log
    appends can leave missing divergence; validation detects this, and further
    annotation is refused rather than silently rewriting history.
    """
    check_observer_value(slot, value, method)
    nonblank("observer", observer)
    nonblank("reason", reason)
    with dataset_lock(dataset) as path:
        _require_modifiable(path)
        find_record(path, entry)
        annotations, divergences = _annotation_snapshot(path)
        current = dict.fromkeys(OBSERVER_SLOTS)
        ids = dict.fromkeys(OBSERVER_SLOTS)
        for row in annotations:
            if row["entry_id"] == entry:
                current[row["slot"]] = row["value"]
                ids[row["slot"]] = row["id"]
        previous = dict(current)
        report = {"id": f"annotation-{len(annotations) + 1:012d}", "entry_id": entry,
                  "slot": slot, "value": value, "observer": observer, "reason": reason,
                  "method": method, "status": "OBSERVED", "timestamp": utc_timestamp()}
        current[slot] = value
        ids[slot] = report["id"]
        event = derive_divergence(report, previous, current, ids, len(divergences) + 1)
        append_jsonl(path / "observers.jsonl", report)
        if event is not None:
            append_jsonl(path / "divergence.jsonl", event)
    return report


def _history_unlocked(path: Path, entry: str, slot: str | None = None) -> list[dict[str, Any]]:
    result = []
    for sequence, row in enumerate(read_jsonl(path / "observers.jsonl"), 1):
        check_annotation(row, sequence)
        if row["entry_id"] == entry and (slot is None or row["slot"] == slot):
            result.append(row)
    return result


def get_annotation_history(dataset: str | os.PathLike[str], entry: str,
                           slot: str | None = None) -> list[dict[str, Any]]:
    """Retrieve all previous reports in append order, optionally by slot."""
    if slot is not None and slot not in OBSERVER_SLOTS:
        raise ValidationError("unknown observer slot")
    with dataset_lock(dataset) as path:
        load_manifest(path)
        find_record(path, entry)
        return _history_unlocked(path, entry, slot)


def show_record(dataset: str | os.PathLike[str], entry: str) -> dict[str, Any]:
    """Return a base run with the latest explicit values materialized in memory."""
    with dataset_lock(dataset) as path:
        load_manifest(path)
        record = dict(find_record(path, entry))
        for row in _history_unlocked(path, entry):
            record[row["slot"]] = row["value"]
        return record


def check_note(row: Mapping[str, Any], sequence: int | None = None) -> None:
    if set(row) != NOTE_FIELDS:
        raise DatasetError("note has missing or undeclared fields")
    nonblank("note id", row["id"])
    if sequence is not None and row["id"] != f"note-{sequence:012d}":
        raise DatasetError("note ID/order is not canonical")
    if not isinstance(row["status"], str) or row["status"] not in STATUSES:
        raise DatasetError("note status must be OBSERVED, DERIVED, or PROPOSED")
    nonblank("note text", row["text"])
    check_timestamp(row["timestamp"])
    if not isinstance(row["references"], list):
        raise DatasetError("note references must be a list")
    for reference in row["references"]:
        nonblank("reference", reference)


def add_note(dataset: str | os.PathLike[str], status: str, text: str,
             references: Iterable[str] = ()) -> dict[str, Any]:
    """Append an explicitly categorized user claim; never promote its status.

    References may target base runs, observer reports, divergence events, or
    earlier notes. Forward/unknown references are refused; repeats are retained.
    """
    if not isinstance(status, str) or status not in STATUSES:
        raise ValidationError("status must be OBSERVED, DERIVED, or PROPOSED")
    nonblank("text", text)
    if isinstance(references, (str, bytes)):
        raise ValidationError("references must be an iterable of IDs, not one string")
    refs = list(references)
    for reference in refs:
        nonblank("reference", reference)
    with dataset_lock(dataset) as path:
        _require_modifiable(path)
        known: set[str] = set()
        for filename in ("records.jsonl", "observers.jsonl", "divergence.jsonl"):
            for row in read_jsonl(path / filename):
                nonblank("log entry id", row.get("id"))
                if row["id"] in known:
                    raise DatasetError(f"duplicate archive ID: {row['id']}")
                known.add(row["id"])
        count = 0
        for count, row in enumerate(read_jsonl(path / "notes.jsonl"), 1):
            check_note(row, count)
            if row["id"] in known or any(ref not in known for ref in row["references"]):
                raise DatasetError("existing notes contain duplicate IDs or unknown references")
            known.add(row["id"])
        unknown = [ref for ref in refs if ref not in known]
        if unknown:
            raise ValidationError(f"unknown references: {', '.join(unknown)}")
        note = {"id": f"note-{count + 1:012d}", "status": status, "text": text,
                "references": refs, "timestamp": utc_timestamp()}
        append_jsonl(path / "notes.jsonl", note)
    return note
