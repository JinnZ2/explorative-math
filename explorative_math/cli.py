"""Stdlib-only command-line interface for finite archives and explicit reports."""

from __future__ import annotations

import argparse
import sys
from typing import Any, Sequence

from .annotations import (add_note, get_annotation_history, observe,
                          parse_observer_value, show_record)
from .core import EnumerationConfig, OBSERVER_SLOTS, STATUSES, canonical_json
from .storage import create_dataset, enumerate_dataset
from .validation import validate_dataset


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="explorative-math", description="Keep every output within declared finite bounds; annotate explicitly.")
    parser.add_argument("--version", action="version", version="%(prog)s 1.0.0")
    commands = parser.add_subparsers(dest="command", required=True)
    enumerate_command = commands.add_parser("enumerate", help="enumerate finite CA, substitution, or integer rules")
    enumerate_command.add_argument("--family", choices=("ca", "substitution", "integer", "all"), default="all")
    enumerate_command.add_argument("--out", required=True)
    enumerate_command.add_argument("--steps", type=int, default=4)
    enumerate_command.add_argument("--width", type=int, default=5)
    enumerate_command.add_argument("--ca-seeds", choices=("all", "single"), default="all")
    enumerate_command.add_argument("--max-word-length", type=int, default=2)
    enumerate_command.add_argument("--substitution-seeds", default="0,1", help="comma-separated nonempty binary strings; repeats retained")
    enumerate_command.add_argument("--a-min", type=int, default=-1)
    enumerate_command.add_argument("--a-max", type=int, default=1)
    enumerate_command.add_argument("--b-min", type=int, default=-1)
    enumerate_command.add_argument("--b-max", type=int, default=1)
    enumerate_command.add_argument("--modulus", type=int, default=8)

    geometry = commands.add_parser("geometry", help="record seeded finite geometry samples and complete paths")
    geometry.add_argument("--out", required=True)
    geometry.add_argument("--experiment", choices=("all", "orthogonality", "shell", "superposition", "trajectories"), default="all")
    geometry.add_argument("--dimensions", default="2,8,32,128")
    geometry.add_argument("--samples", type=int, default=32)
    geometry.add_argument("--vectors", type=int, default=0, help="0 means 2*dimension; explicit k must exceed dimension for superposition")
    geometry.add_argument("--steps", type=int, default=32)
    geometry.add_argument("--seed", type=int, default=0)
    geometry.add_argument("--shell-thickness", type=float, default=0.1)
    geometry.add_argument("--map-scale", type=float, default=0.8)

    project = commands.add_parser("project", help="record finite rotation/projection/composition frames, matches and misses")
    project.add_argument("--out", required=True)
    project.add_argument("--dimension", type=int, default=5)
    project.add_argument("--shape", choices=("crosspolytope", "hypercube", "simplex"), default="crosspolytope")
    project.add_argument("--frames", type=int, default=12)
    project.add_argument("--seed", type=int, default=0)
    project.add_argument("--composition", choices=("none", "phi", "second"), default="phi")
    project.add_argument("--tolerance-degrees", type=float, default=1.0)
    project.add_argument("--max-vertices", type=int, default=4096)

    annotation = commands.add_parser("observe", help="append a named explicit report; never compute an observer")
    annotation.add_argument("--dataset", required=True)
    annotation.add_argument("--entry", required=True)
    annotation.add_argument("--slot", choices=OBSERVER_SLOTS, required=True)
    annotation.add_argument("--value", required=True, help="true/false/null, numeric/null, or divergence_note text/null")
    annotation.add_argument("--observer", required=True)
    annotation.add_argument("--reason", required=True)
    annotation.add_argument("--method", help="required for numeric compressibility")

    show = commands.add_parser("show", help="materialize latest explicit observer values without editing base records")
    show.add_argument("--dataset", required=True)
    show.add_argument("--entry", required=True)
    history = commands.add_parser("history", help="retrieve every previous explicit report")
    history.add_argument("--dataset", required=True)
    history.add_argument("--entry", required=True)
    history.add_argument("--slot", choices=OBSERVER_SLOTS)

    note = commands.add_parser("note", help="append a user claim with its explicit status unchanged")
    note.add_argument("--dataset", required=True)
    note.add_argument("--status", choices=sorted(STATUSES), required=True)
    note.add_argument("--text", required=True)
    note.add_argument("--references", help="comma-separated existing IDs; no forward references")

    validate = commands.add_parser("validate", help="validate completion, schema, coverage, IDs, and annotation derivation")
    validate.add_argument("--dataset", required=True)
    validate.add_argument("--allow-incomplete", action="store_true", help="check a sound prefix without requiring completion")
    return parser


def _config_args(args: argparse.Namespace) -> dict[str, Any]:
    return {key: value for key, value in vars(args).items() if key not in ("command", "out")}


def _run(args: argparse.Namespace) -> tuple[Any, int]:
    command = args.command
    if command == "enumerate":
        values = _config_args(args)
        values["substitution_seeds"] = tuple(seed.strip() for seed in args.substitution_seeds.split(","))
        return enumerate_dataset(args.out, EnumerationConfig(**values)), 0
    if command == "geometry":
        from .geometry import count_geometry, iter_geometry_runs, validate_geometry_config
        config = validate_geometry_config(_config_args(args))
        counts = count_geometry(config)
        return create_dataset(args.out, config, iter_geometry_runs(config), counts["runs"],
                              counts["outputs"], source="geometry"), 0
    if command == "project":
        from .projections import count_projection, iter_projection_runs, validate_projection_config
        config = validate_projection_config(_config_args(args))
        counts = count_projection(config)
        return create_dataset(args.out, config, iter_projection_runs(config), counts["runs"],
                              counts["outputs"], source="projection"), 0
    if command == "observe":
        value = parse_observer_value(args.slot, args.value)
        return observe(args.dataset, args.entry, args.slot, value,
                       args.observer, args.reason, args.method), 0
    if command == "show":
        return show_record(args.dataset, args.entry), 0
    if command == "history":
        return get_annotation_history(args.dataset, args.entry, args.slot), 0
    if command == "note":
        references = [] if args.references is None else [ref.strip() for ref in args.references.split(",")]
        return add_note(args.dataset, args.status, args.text, references), 0
    if command == "validate":
        report = validate_dataset(args.dataset, allow_incomplete=args.allow_incomplete)
        return report, 0 if report["valid"] else 1
    raise ValueError(f"unknown command: {command}")


def main(argv: Sequence[str] | None = None) -> int:
    """Run a command, emit canonical JSON, and return a concise exit status."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        result, code = _run(args)
        print(canonical_json(result))
        return code
    except BrokenPipeError:
        return 1
    except KeyboardInterrupt:
        print("error: interrupted", file=sys.stderr)
        return 130
    except (OSError, ValueError, TypeError, KeyError, ImportError, RecursionError, OverflowError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
