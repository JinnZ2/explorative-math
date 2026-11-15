"""h. three_gap_check -- the three_gap_sum probe over every logged cut run.

CC0. Co-authored. Stdlib only.

Reads every runs/*.json (except its own output), finds each
`distinct_gap_lengths` entry wherever it sits, and applies
probes.three_gap_sum: with exactly three lengths, is the longest the sum of
the other two? Each entry gets one row: pass, fail or not_applicable (any
other count of lengths). No row is dropped, and no run file is changed.

  python3 experiments/three_gap_check.py [seed]
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "three_gap_check"
KEY = "distinct_gap_lengths"
TOL = 1e-6


def find(node, path=""):
    """Yield (json path, value) for every KEY in a nested structure."""
    if isinstance(node, dict):
        for k, v in node.items():
            here = f"{path}.{k}" if path else k
            if k == KEY:
                yield here, v
            else:
                yield from find(v, here)
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from find(v, f"{path}[{i}]")


def as_lengths(value):
    """A logged entry is a list of numbers, a list of {"length": ...}, or null."""
    if value is None:
        return None
    return [v["length"] if isinstance(v, dict) else v for v in value]


def run(seed=0, runs_dir=None):
    d = Path(runs_dir) if runs_dir else emlib.RUNS
    rows = []
    for f in emlib.run_files(d):
        if f.name.startswith(NAME + "_"):
            continue
        rec = emlib.load_run(f)
        for where, value in find(rec.get("outputs")):
            ls = as_lengths(value)
            row = {"run": f.name, "where": where, "lengths": ls}
            if ls is None:
                row.update({"n_lengths": None, "verdict": "not_applicable",
                            "residual": None, "tol": TOL})
            else:
                row.update(Q.three_gap_sum(ls, TOL))
            rows.append(row)
    counts = {}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    by_n = {}
    for r in rows:
        k = str(r["n_lengths"])
        by_n[k] = by_n.get(k, 0) + 1
    return {"rows": rows, "verdict_counts": counts, "rows_by_n_lengths": by_n,
            "runs_with_an_entry": sorted({r["run"] for r in rows})}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    symbols = [{"pass": "p", "fail": "f", "not_applicable": "n"}[r["verdict"]]
               for r in out["rows"]]
    comp = emlib.compressibility(symbols, "one letter per row: p pass, f fail, n not_applicable")
    fails = [f"{r['run']} {r['where']}" for r in out["rows"] if r["verdict"] == "fail"]
    ai = (f"OBSERVED: {len(out['rows'])} gap-length entries in "
          f"{len(out['runs_with_an_entry'])} run files; verdicts {out['verdict_counts']}; "
          f"entries by number of lengths {out['rows_by_n_lengths']}; fails: {fails or 'none'}. "
          "AI reading: not interpretable atm beyond the counts.")
    p = emlib.write_run(NAME, seed, {"key": KEY, "tol": TOL, "runs_dir": "runs"},
                        out, None, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
