"""f. weighted_window -- one Z^2 cut, two windows: hard strip and Gaussian weight.

CC0. Co-authored. Stdlib only.

The same lattice (Z^2) and the same cut (a line through the origin, slope
1/phi) are projected twice.

  hard      keep a lattice point when its perpendicular coordinate q lies in
            [lo, hi), the projection of the unit square (strip_window_1d).
            Every kept point counts 1.
  gaussian  keep every lattice point, with weight exp(-q^2 / (2 sigma^2)):
            a function of its distance |q| to the cut. sigma is declared in
            params (half the hard strip's width).

Both are logged. The Gaussian set has no edge, so it is recorded down to a
declared weight floor (the floor is a recording limit, not a window), and
read at several weight levels: the points with weight >= level.

  python3 experiments/weighted_window.py [seed]
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import primitives as P  # noqa: E402

NAME = "weighted_window"
EXTENT = 60                 # lattice points with |m|, |n| <= EXTENT
OFFSET = 1e-7               # shifts the hard window off lattice-point boundaries
WEIGHT_FLOOR = 1e-6         # Gaussian points below this weight are not written
LEVELS = (0.9, math.exp(-0.5), 0.5, 0.1, 0.01)
SEP = 1e-6                  # gap lengths closer than this are one length


def clusters(gaps, sep=SEP):
    out = []
    for g in sorted(gaps):
        if out and g - out[-1][-1] <= sep:
            out[-1].append(g)
        else:
            out.append([g])
    return [{"length": sum(c) / len(c), "count": len(c)} for c in out]


def gaps_of(xs):
    xs = sorted(xs)
    return [b - a for a, b in zip(xs, xs[1:])]


def word_of(gaps, lengths):
    """One letter per gap: a = shortest length, b = next, ..."""
    return "".join(chr(ord("a") + min(range(len(lengths)),
                                      key=lambda i: abs(g - lengths[i]["length"])))
                   for g in gaps)


def run(seed=0, slope=1 / P.PHI, extent=EXTENT):
    par, perp = C.line_basis_2d(slope)
    lattice = P.lattice_zn(2, extent)
    lo, hi = C.strip_bounds_1d(perp[0], OFFSET)
    sigma = (hi - lo) / 2.0
    margin = extent * 0.5 * par[0][0]   # stretch where the lattice box cuts nothing off

    hard = C.cut_and_project(lattice, par, perp, C.strip_window_1d(perp[0], OFFSET))
    hx = sorted(p[0] for p in hard if -margin <= p[0] <= margin)
    hg = gaps_of(hx)
    hl = clusters(hg)

    recs = C.cut_and_project(lattice, par, perp, weight=C.gaussian_weight(sigma))
    inside = [r for r in recs if -margin <= r["point"][0] <= margin]
    total = sum(r["weight"] for r in inside)
    in_strip = sum(r["weight"] for r in inside if lo <= r["perp"][0] < hi)
    kept = sorted((r for r in inside if r["weight"] >= WEIGHT_FLOOR),
                  key=lambda r: r["point"][0])
    levels = []
    for t in LEVELS:
        xs = [r["point"][0] for r in kept if r["weight"] >= t]
        g = gaps_of(xs)
        ls = clusters(g)
        levels.append({"level": t,
                       "half_width_in_q": sigma * math.sqrt(2.0 * math.log(1.0 / t)),
                       "n_points": len(xs),
                       "weight_sum": sum(r["weight"] for r in kept if r["weight"] >= t),
                       "distinct_gap_lengths": ls,
                       "word": word_of(g, ls)})
    return {
        "slope": slope, "parallel": par[0], "perpendicular": perp[0],
        "margin": margin,
        "hard": {"window": [lo, hi], "n_points": len(hx), "points": hx,
                 "distinct_gap_lengths": hl, "word": word_of(hg, hl)},
        "gaussian": {"sigma": sigma, "center": 0.0,
                     "n_lattice_points_in_margin": len(inside),
                     "weight_sum_all": total,
                     "weight_sum_inside_hard_window": in_strip,
                     "weight_share_inside_hard_window": in_strip / total,
                     "n_points_above_floor": len(kept),
                     "points_x_q_weight": [[r["point"][0], r["perp"][0], r["weight"]]
                                           for r in kept],
                     "levels": levels},
    }


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    h, g = out["hard"], out["gaussian"]
    comp = emlib.compressibility(
        list(h["word"]), "hard-window chain, one letter per gap (a = shortest length)")
    per_level = "; ".join(f"level {lv['level']:.4g}: {lv['n_points']} points, "
                          f"{len(lv['distinct_gap_lengths'])} gap lengths"
                          for lv in g["levels"])
    ai = (f"OBSERVED: hard window {h['n_points']} points, "
          f"{len(h['distinct_gap_lengths'])} gap lengths. Gaussian (sigma = half the hard "
          f"width, centred on the cut): share of total weight inside the hard strip "
          f"{g['weight_share_inside_hard_window']:.6f}; {per_level}. "
          "DERIVED: weight >= level is the same as |q| <= sigma*sqrt(2 ln(1/level)), so each "
          "level is a hard window centred on the cut; the Gaussian set is that whole family "
          "at once, each member carrying its weight. AI reading: the count of gap lengths "
          "changes with the level, so 'the gap lengths' of the weighted set is not one "
          "number; which single reading of a weighted set stands beside the hard-window "
          "word is not interpretable atm.")
    p = emlib.write_run(NAME, seed,
                        {"slope": "1/phi", "extent": EXTENT, "window_offset": OFFSET,
                         "windows": ["hard strip (unit-square projection)",
                                     "gaussian weight by distance to cut"],
                         "gaussian_sigma": "half the hard strip width",
                         "weight_floor": WEIGHT_FLOOR, "levels": list(LEVELS),
                         "gap_cluster_sep": SEP},
                        out, None, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
