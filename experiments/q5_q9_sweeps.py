"""Q5 (window width sweep) and Q9 (half-mass dimension seed sweep) for JinnZ2/explorative-math.
CC0. Co-authored. Stdlib only. Imports the repo read-only. Predictions: docs/prereg/prereg_q12_q11_q5_q9_2026-11-18.json (sha256 6ae97bcc...606c).
  python3 experiments/q5_q9_sweeps.py q5 [--small]
  python3 experiments/q5_q9_sweeps.py q9 [--small]
--small: q5 grid step 0.05 (not 0.01); q9 10 seeds (not 100). Own run files (name + _small).
"""
import math
import random
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = None                     # None = the repo runs/ directory
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "experiments"))
from toolkit import loghook  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402
import window_extent as WE  # noqa: E402

SLOPE = 1.0 / P.PHI
SEP = 1e-6
TARGET = 1200


def gap_lengths(h, target=TARGET):
    """distinct gap lengths of the hard window |q| <= h on Z^2 cut at SLOPE, ~target points."""
    width = 2.0 * h
    par, perp = C.line_basis_2d(SLOPE)
    x_max = target / (2.0 * width)
    pts = C.cut_and_project(WE.band_lattice(SLOPE, h, x_max), par, perp, lambda q: abs(q[0]) <= h)
    xs = sorted(p[0] for p in pts if abs(p[0]) <= x_max)
    gaps = [b - a for a, b in zip(xs, xs[1:])]
    cl = WE.clusters(gaps, SEP)
    return len(xs), cl


def q5(small=False):
    par, perp = C.line_basis_2d(SLOPE)
    lo, hi = C.strip_bounds_1d(perp[0])
    canon = hi - lo
    rows = []
    rs = [round(0.05 + (0.05 if small else 0.01) * i, 4) for i in range(0, 40 if small else 196)]   # 0.05 .. 2.00 (small: 0.05 .. 2.00 step 0.05)
    rs += [0.99, 0.995, 0.999, 1.001, 1.005, 1.01]
    for r in sorted(set(rs)):
        h = r * canon / 2.0
        n, cl = gap_lengths(h)
        total = sum(c["count"] for c in cl)
        counts = sorted(c["count"] for c in cl)
        rows.append({"r": r, "width_q": r * canon, "n_points": n, "n_gap_lengths": len(cl),
                     "lengths": [round(c["length"], 9) for c in cl],
                     "gap_counts": counts, "min_gap_frequency": (counts[0] / total) if total else None,
                     "three_gap_sum": Q.three_gap_sum([c["length"] for c in cl]) if len(cl) == 3 else None})
    by_n = {}
    for r in rows:
        by_n.setdefault(r["n_gap_lengths"], []).append(r["r"])
    two_ranges = []
    for r in sorted(by_n.get(2, [])):
        if two_ranges and abs(r - two_ranges[-1][1]) < 0.0115:
            two_ranges[-1][1] = r
        else:
            two_ranges.append([r, r])
    outputs = {"canonical_width_q": canon, "slope": SLOPE, "rows": rows,
               "r_values_by_gap_count": {str(k): v for k, v in sorted(by_n.items())},
               "r_ranges_with_two_lengths_grid_0.01": two_ranges,
               "max_gap_count": max(by_n), "r1_gap_count": next(x["n_gap_lengths"] for x in rows if x["r"] == 1.0),
               "rows_with_min_gap_frequency_below_0.01": [x["r"] for x in rows if x["min_gap_frequency"] is not None and x["min_gap_frequency"] < 0.01]}
    params = {"experiment": "q5_width_sweep", "slope": "1/phi", "target_points": TARGET, "sep": SEP,
              "r_grid": "0.05..2.00 step 0.01 plus 0.99,0.995,0.999,1.001,1.005,1.01",
              "prereg": "docs/prereg/prereg_q12_q11_q5_q9_2026-11-18.json", "repo_commit": "3162c9950136a4a7f14e1e064f8d335db33968c8"}
    params["small"] = small
    return loghook.log_run("q5_width_sweep" + ("_small" if small else ""), 0, params, outputs, out_dir=OUT), outputs


# ------------------------------------------------------------------ Q9
DIMS = list(range(30, 101, 5))
PAIRS = 400
BAND = math.sin(math.radians(5.0))      # |cos| < sin 5 deg  <=>  angle within 90 +/- 5 deg


def share(rng, d):
    k = 0
    for _ in range(PAIRS):
        u, v = P.random_unit(rng, d), P.random_unit(rng, d)
        if abs(P.dot(u, v)) < BAND:
            k += 1
    return k / PAIRS


def fit_cross(ds, ps):
    n = len(ds)
    mx, my = sum(ds) / n, sum(ps) / n
    sxx = sum((x - mx) ** 2 for x in ds)
    b = sum((x - mx) * (y - my) for x, y in zip(ds, ps)) / sxx
    a = my - b * mx
    return (0.5 - a) / b, b


def q9(seeds=100, small=False):
    seeds = 10 if small else seeds
    per_seed = []
    for s in range(seeds):
        rng = random.Random(s)
        ps = [share(rng, d) for d in DIMS]
        cross, slope = fit_cross(DIMS, ps)
        first = next((d for d, p in zip(DIMS, ps) if p >= 0.5), None)
        per_seed.append({"seed": s, "shares": ps, "fit_crossing": cross, "fit_slope_per_dim": slope,
                         "first_grid_dim_at_or_above_0.5": first})
        if s % 10 == 9:
            print("q9 seed", s, flush=True)
    xs = [r["fit_crossing"] for r in per_seed]
    outputs = {"dims": DIMS, "pairs": PAIRS, "per_seed": per_seed,
               "pooled_mean_crossing": statistics.fmean(xs), "between_seed_sd": statistics.stdev(xs),
               "min": min(xs), "max": max(xs), "median": statistics.median(xs),
               "mean_slope_per_dim": statistics.fmean(r["fit_slope_per_dim"] for r in per_seed),
               "exact_curve_crossing_scratch": "between n=61 (0.4994) and n=62 (0.5030), unlogged scratch, Simpson integral of (1-t^2)^((n-3)/2)",
               "prediction_pooled_in_58.3_64.3": 58.3 <= statistics.fmean(xs) <= 64.3,
               "prediction_sd_in_4_10": 4 <= statistics.stdev(xs) <= 10}
    params = {"experiment": "q9_half_mass", "dims": DIMS, "pairs": PAIRS, "seeds": seeds,
              "prereg": "docs/prereg/prereg_q12_q11_q5_q9_2026-11-18.json", "repo_commit": "3162c9950136a4a7f14e1e064f8d335db33968c8"}
    params["small"] = small
    return loghook.log_run("q9_half_mass" + ("_small" if small else ""), 0, params, outputs, out_dir=OUT), outputs


if __name__ == "__main__":
    which = sys.argv[1]
    small = "--small" in sys.argv
    if which == "q5":
        p, o = q5(small)
        print(p)
        print("max gap count:", o["max_gap_count"], "| r=1.0 ->", o["r1_gap_count"])
        print("counts by n:", {k: len(v) for k, v in o["r_values_by_gap_count"].items()})
        print("two-length r ranges:", o["r_ranges_with_two_lengths_grid_0.01"])
        print("low min-gap-frequency r:", o["rows_with_min_gap_frequency_below_0.01"][:20])
        print("r in [0.9,1.1] n:", [(x["r"], x["n_gap_lengths"]) for x in o["rows"] if 0.9 <= x["r"] <= 1.1])
    else:
        p, o = q9(small=small)
        print(p)
        print({k: o[k] for k in ("pooled_mean_crossing", "between_seed_sd", "min", "max", "median", "mean_slope_per_dim",
                                 "prediction_pooled_in_58.3_64.3", "prediction_sd_in_4_10")})
