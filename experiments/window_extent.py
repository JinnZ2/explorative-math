"""g. window_extent -- one logged window, a longer stretch, and its width
against the integer combinations m + n*slope.

CC0. Co-authored. Stdlib only.

Source: runs/generated_cut_weighted_seed1.json, weight level 0.9. That level
is a hard window centred on the cut, |q| <= h. On the 57 points of that run
it read two gap lengths. Here the same slope and the same h are used on a
stretch of about 2000 points, and the gap-length count is logged at nested
stretch sizes.

The window width is also set beside the numbers m + n*slope, |m|, |n| <= 10.
The perpendicular coordinates of Z^2 are (n - slope*m) / sqrt(1 + slope^2),
so two scalings of the width are logged: as measured in q, and multiplied by
sqrt(1 + slope^2) (the scaling in which the lattice's perpendicular
coordinates are exactly the numbers m + n*slope). The canonical strip width
is logged the same way as a control.

Nothing is interpreted here. Both readings are logged.

  python3 experiments/window_extent.py [seed]
"""
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "window_extent"
SOURCE = "runs/generated_cut_weighted_seed1.json"
LEVEL_INDEX = 0               # level 0.9
TARGETS = (57, 100, 200, 500, 1000, 2000)   # points per nested stretch (approximate)
COMBO_MAX = 10
SEP = 1e-6                    # gap lengths closer than this are one length


def source_window(path=None):
    rec = json.loads((Path(path) if path else emlib.ROOT / SOURCE).read_text())
    spec = rec["params"]["spec"]
    slope = spec["slope"]
    _, perp = C.line_basis_2d(slope)
    lo, hi = C.strip_bounds_1d(perp[0])
    sigma = spec["sigma_over_strip_width"] * (hi - lo)
    level = spec["levels"][LEVEL_INDEX]
    h = sigma * math.sqrt(2.0 * math.log(1.0 / level))
    return {"slope": slope, "level": level, "sigma": sigma, "half_width_in_q": h,
            "canonical_strip_width": hi - lo}


def band_lattice(slope, h, x_max):
    """Points of Z^2 near the line: every (m, n) that can satisfy |q| <= h and
    |x| <= x_max, plus one row of slack either side. cut_and_project applies
    the window itself."""
    nrm = math.hypot(1.0, slope)
    m_max = int(x_max * nrm + h * nrm) + 2
    out = []
    for m in range(-m_max, m_max + 1):
        c = slope * m
        for n in range(math.floor(c - h * nrm) - 1, math.ceil(c + h * nrm) + 2):
            out.append((float(m), float(n)))
    return out


def clusters(gaps, sep=SEP):
    out = []
    for g in sorted(gaps):
        if out and g - out[-1][-1] <= sep:
            out[-1].append(g)
        else:
            out.append([g])
    return [{"length": sum(c) / len(c), "count": len(c)} for c in out]


def nearest_combos(value, slope, k=COMBO_MAX, keep=3):
    rows = [{"m": m, "n": n, "value": m + n * slope, "residual": value - (m + n * slope)}
            for m in range(-k, k + 1) for n in range(-k, k + 1)]
    rows.sort(key=lambda r: abs(r["residual"]))
    return rows[:keep]


def run(seed=0, source=None, targets=TARGETS):
    w = source_window(source)
    slope, h = w["slope"], w["half_width_in_q"]
    width = 2.0 * h
    nrm = math.hypot(1.0, slope)
    par, perp = C.line_basis_2d(slope)
    x_max = max(targets) / (2.0 * width)      # point density along the line = width
    pts = C.cut_and_project(band_lattice(slope, h, x_max), par, perp,
                            lambda q: abs(q[0]) <= h)
    xs = sorted(p[0] for p in pts if abs(p[0]) <= x_max)
    stretches = []
    for t in targets:
        half = t / (2.0 * width)
        sub = [x for x in xs if abs(x) <= half]
        ls = clusters([b - a for a, b in zip(sub, sub[1:])])
        stretches.append({"target_points": t, "half_length": half, "n_points": len(sub),
                          "n_gap_lengths": len(ls), "distinct_gap_lengths": ls,
                          "three_gap_sum": Q.three_gap_sum([c["length"] for c in ls])})
    full = clusters([b - a for a, b in zip(xs, xs[1:])])
    gaps = [b - a for a, b in zip(xs, xs[1:])]
    word = "".join(chr(ord("a") + min(range(len(full)),
                                      key=lambda i: abs(g - full[i]["length"])))
                   for g in gaps)
    first = {chr(ord("a") + i): word.find(chr(ord("a") + i)) for i in range(len(full))}

    def beside(label, value):
        return {"label": label, "value": value, "nearest": nearest_combos(value, slope)}

    cw = w["canonical_strip_width"]
    return {
        "source_window": w, "width_in_q": width, "norm": nrm,
        "width_over_canonical": width / cw,
        "points": xs, "word": word, "first_index_of_each_letter": first,
        "stretches": stretches,
        "combinations": {
            "definition": f"m + n*slope, |m|, |n| <= {COMBO_MAX}; residual = value - combination; "
                          "three nearest kept, nearest first",
            "window_width_in_q": beside("2h", width),
            "window_width_times_norm": beside("2h * sqrt(1 + slope^2)", width * nrm),
            "control_canonical_width_in_q": beside("hi - lo", cw),
            "control_canonical_width_times_norm": beside("(hi - lo) * sqrt(1 + slope^2)",
                                                         cw * nrm),
        },
    }


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    comp = emlib.compressibility(
        list(out["word"]), "full stretch, one letter per gap (a = shortest length)")
    st = "; ".join(f"{s['n_points']} points: {s['n_gap_lengths']}" for s in out["stretches"])
    c = out["combinations"]

    def near(key):
        r = c[key]["nearest"][0]
        return f"{c[key]['value']:.9f} nearest m={r['m']} n={r['n']} residual {r['residual']:.3e}"
    ai = (f"OBSERVED: gap-length count by stretch -- {st}. First index of each gap letter "
          f"in the full word: {out['first_index_of_each_letter']}. Width beside m + n*slope: "
          f"in q, {near('window_width_in_q')}; times norm, {near('window_width_times_norm')}; "
          f"control (canonical width) in q, {near('control_canonical_width_in_q')}; times "
          f"norm, {near('control_canonical_width_times_norm')}. "
          "No reading is given: the run was asked for as a log, not an interpretation.")
    p = emlib.write_run(NAME, seed,
                        {"source_run": SOURCE, "level_index": LEVEL_INDEX,
                         "targets": list(TARGETS), "combo_max": COMBO_MAX,
                         "gap_cluster_sep": SEP},
                        out, None, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
