"""i. min_stretch_flag -- tag logged gap readings taken on too short a stretch.

CC0. Co-authored. Stdlib only.

For every `distinct_gap_lengths` entry in runs/*.json, find how many points
the reading saw and the frequency f of the third (rarest) gap length, and
apply probes.min_stretch_tag: fewer than FACTOR / f points -> "under-resolved".

Where f comes from, in this order (the row says which):

  reference   a longer logged reading of the SAME window, named in
              REFERENCES below (run file + path). Declared by hand.
  group       the longest reading among entries of one run that share a
              window (the nested stretches of window_extent).
  self        the entry's own counts, when it shows three lengths. A
              self-estimate is bounded by the same stretch it is judging.
  void        two lengths (or another count) and no reference: f is not
              known. The row is kept, not assessed, and carries no tag.

No run file is changed; the rows go to this experiment's own run file.

  python3 experiments/min_stretch_flag.py [seed]
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "min_stretch_flag"
KEY = "distinct_gap_lengths"
FACTOR = 3.0
SKIP = (NAME + "_", "three_gap_check_")      # run files that hold rows, not readings
# (run, where) of a reading -> (run, where) of a longer reading of the same window
REFERENCES = {
    ("generated_cut_weighted_seed1.json", "levels[0]"):
        ("window_extent_seed0.json", "stretches[5]"),
}
# run -> path prefix whose entries share one window; f is taken from the longest
GROUPS = {"window_extent_seed0.json": "stretches"}


def find(node, path=""):
    """Yield (path, dict) for every dict that holds KEY."""
    if isinstance(node, dict):
        if KEY in node:
            yield path, node
        for k, v in node.items():
            yield from find(v, f"{path}.{k}" if path else k)
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from find(v, f"{path}[{i}]")


def reading(entry):
    """(n_points, counts per length shortest first) from one logged entry.
    Either may be None when the entry does not hold it."""
    value = entry[KEY]
    word = entry.get("word")
    counts = None
    if value and isinstance(value[0], dict):
        counts = [v["count"] for v in value]
    elif value and word:
        letters = sorted(set(word))
        counts = [word.count(c) for c in letters]
        if len(counts) != len(value):
            counts = None
    if word:
        n = len(word) + 1
    elif counts:
        n = sum(counts) + 1
    else:
        n = entry.get("n_points")
        if n is None and entry.get("n_gaps") is not None:
            n = entry["n_gaps"] + 1
    return n, counts


def freq(counts):
    if not counts or len(counts) != 3:
        return None
    return min(counts) / sum(counts)


def run(seed=0, runs_dir=None, references=None, groups=None):
    d = Path(runs_dir) if runs_dir else emlib.RUNS
    references = REFERENCES if references is None else references
    groups = GROUPS if groups is None else groups
    entries = {}
    for f in emlib.run_files(d):
        if f.name.startswith(SKIP):
            continue
        for where, entry in find(emlib.load_run(f).get("outputs")):
            n, counts = reading(entry)
            entries[(f.name, where)] = {"n_points": n, "counts": counts,
                                        "n_lengths": len(entry[KEY]) if entry[KEY] is not None
                                        else None}
    rows = []
    for (name, where), e in entries.items():
        src, ref = "void", None
        f3 = None
        if (name, where) in references:
            ref = references[(name, where)]
            f3, src = freq(entries[ref]["counts"]) if ref in entries else None, "reference"
        elif name in groups and where.startswith(groups[name]):
            members = [k for k in entries if k[0] == name and k[1].startswith(groups[name])]
            ref = max(members, key=lambda k: entries[k]["n_points"] or 0)
            f3, src = freq(entries[ref]["counts"]), "group"
        elif freq(e["counts"]) is not None:
            f3, src = freq(e["counts"]), "self"
        row = {"run": name, "where": where, "n_lengths": e["n_lengths"],
               "counts": e["counts"], "freq_source": src,
               "freq_reference": list(ref) if ref else None}
        row.update(Q.min_stretch_tag(e["n_points"], f3, FACTOR))
        rows.append(row)
    tagged = [r for r in rows if r["tag"]]
    return {"rows": rows,
            "n_rows": len(rows),
            "n_assessed": sum(1 for r in rows if r["assessed"]),
            "n_not_assessed": sum(1 for r in rows if not r["assessed"]),
            "n_under_resolved": len(tagged),
            "under_resolved": [{"run": r["run"], "where": r["where"],
                                "n_points": r["n_points"], "min_points": r["min_points"],
                                "n_lengths": r["n_lengths"], "freq_source": r["freq_source"]}
                               for r in tagged]}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    symbols = ["u" if r["tag"] else ("a" if r["assessed"] else "v") for r in out["rows"]]
    comp = emlib.compressibility(
        symbols, "one letter per row: u under-resolved, a assessed and not tagged, "
                 "v not assessed (frequency void)")
    listed = "; ".join(f"{r['run']} {r['where']} ({r['n_points']} points, needs "
                       f"{r['min_points']:.0f}, reads {r['n_lengths']} lengths, f from "
                       f"{r['freq_source']})" for r in out["under_resolved"]) or "none"
    ai = (f"OBSERVED: {out['n_rows']} gap readings; {out['n_assessed']} assessed, "
          f"{out['n_not_assessed']} not assessed (third-gap frequency void); "
          f"{out['n_under_resolved']} tagged under-resolved at factor {FACTOR}: {listed}. "
          "DERIVED: a reading with two lengths and no longer reference cannot be assessed "
          "by this flag, because the frequency it needs is the thing the short stretch may "
          "have missed. AI reading: not interpretable atm beyond the rows.")
    p = emlib.write_run(NAME, seed,
                        {"key": KEY, "factor": FACTOR, "runs_dir": "runs",
                         "references": [[list(k), list(v)] for k, v in REFERENCES.items()],
                         "groups": GROUPS, "skipped_prefixes": list(SKIP)},
                        out, None, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
