"""Finite explorative mathematics with explicit, append-only observer reports.

Public archive APIs:
    enumerate_dataset(out, config=None) -> manifest
    create_dataset(out, config, runs, expected_runs, expected_outputs,
                   *, source='custom') -> manifest
    observe(dataset, entry, slot, value, observer, reason, method=None) -> report
    show_record(dataset, entry) -> materialized run
    get_annotation_history(dataset, entry, slot=None) -> list of reports
    add_note(dataset, status, text, references=()) -> note
    validate_dataset(dataset, *, allow_incomplete=False) -> validation report

The CLI is available with ``python -m explorative_math``. No observers or
compressibility metrics run automatically; null always means unknown.
"""

from .core import (
    EnumerationConfig, GENERATOR_VERSION, OBSERVER_SLOTS, SCHEMA_VERSION,
    STATUSES, ValidationError, canonical_json, evolve_ca, evolve_integer,
    iter_binary_words, iter_ca_rules, iter_ca_seeds, iter_ca_trajectory,
    iter_integer_rules, iter_integer_seeds, iter_integer_trajectory,
    iter_records, iter_runs, iter_substitution_rules,
    iter_substitution_trajectory, stable_id, substitute,
)
from .storage import DatasetError, create_dataset, enumerate_dataset
from .annotations import add_note, get_annotation_history, observe, show_record
from .validation import validate_dataset
from .geometry import count_geometry, iter_geometry_runs, validate_geometry_config

__version__ = GENERATOR_VERSION

__all__ = [
    "EnumerationConfig", "ValidationError", "DatasetError", "SCHEMA_VERSION",
    "GENERATOR_VERSION", "OBSERVER_SLOTS", "STATUSES", "canonical_json", "stable_id",
    "iter_ca_rules", "iter_ca_seeds", "evolve_ca", "iter_ca_trajectory",
    "iter_binary_words", "iter_substitution_rules", "substitute",
    "iter_substitution_trajectory", "iter_integer_rules", "iter_integer_seeds",
    "evolve_integer", "iter_integer_trajectory", "iter_runs", "iter_records",
    "create_dataset", "enumerate_dataset", "observe", "show_record",
    "get_annotation_history", "add_note", "validate_dataset",
    "validate_geometry_config", "count_geometry", "iter_geometry_runs",
]
