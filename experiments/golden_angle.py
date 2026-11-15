"""d. golden_angle -- phyllotaxis points: golden angle vs perturbed angles.

CC0. Co-authored. Stdlib only.

Place point k at radius sqrt(k) and angle k * alpha (Vogel's model). For
each divergence angle alpha, measure how well the points pack:
  min_nn         smallest nearest-neighbour distance
  nn_cv          coefficient of variation of nearest-neighbour distances
                 (low = even spacing)
Angles: the golden angle 360(1 - 1/phi), perturbations around it, and
rational fractions of a turn (spokes expected).

  python3 experiments/golden_angle.py [seed]
"""
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
import nature_catalog  # noqa: E402

NAME = "golden_angle"
PHI = (1 + 5 ** 0.5) / 2
GOLDEN = 360.0 * (1.0 - 1.0 / PHI)
N_POINTS = 400
INNER_SKIP = 20  # drop the first points, where the spiral has not formed
ANGLES = {
    "golden": GOLDEN,
    "golden-0.1": GOLDEN - 0.1,
    "golden+0.1": GOLDEN + 0.1,
    "golden-0.5": GOLDEN - 0.5,
    "golden+0.5": GOLDEN + 0.5,
    "137.0": 137.0,
    "138.0": 138.0,
    "rational_2/5_turn": 144.0,
    "rational_3/8_turn": 135.0,
}


def points(alpha_deg, n=N_POINTS):
    a = math.radians(alpha_deg)
    return [(math.sqrt(k) * math.cos(k * a), math.sqrt(k) * math.sin(k * a))
            for k in range(1, n + 1)]


def nn_distances(pts, skip=INNER_SKIP):
    out = []
    for i in range(skip, len(pts)):
        xi, yi = pts[i]
        best = min(math.hypot(xi - xj, yi - yj) for j, (xj, yj) in enumerate(pts) if j != i)
        out.append(best)
    return out


def run(seed=0):
    rows = {}
    for label, alpha in ANGLES.items():
        nn = nn_distances(points(alpha))
        rows[label] = {"alpha_deg": alpha, "min_nn": min(nn),
                       "mean_nn": statistics.fmean(nn),
                       "nn_cv": statistics.pstdev(nn) / statistics.fmean(nn)}
    return rows


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    rows = run(seed)
    order = sorted(rows, key=lambda k: rows[k]["nn_cv"])
    comp = emlib.compressibility(
        [chr(ord("a") + order.index(k)) for k in ANGLES],
        "rank of each angle by nn_cv (a = most even), in ANGLES order")
    ai = (f"OBSERVED: ranking by evenness (nn_cv, low first): {order}. "
          f"Golden min_nn = {rows['golden']['min_nn']:.4f}; 144 deg (2/5 turn) "
          f"min_nn = {rows['rational_2/5_turn']['min_nn']:.4f}. AI reading: not "
          "interpretable atm whether the small perturbations (+/-0.1) differ from "
          "golden by more than the finite-N effect; that needs an N sweep.")
    p = emlib.write_run(NAME, seed,
                        {"n_points": N_POINTS, "inner_skip": INNER_SKIP,
                         "model": "Vogel r=sqrt(k), theta=k*alpha"},
                        {"rows": rows, "rank_by_nn_cv": order},
                        {k: v for k, v in ANGLES.items()}, ai, comp,
                        nature_catalog.match("phyllotaxis_golden_angle"))
    print(p)
    for k in order:
        r = rows[k]
        print(f"{k:18s} alpha={r['alpha_deg']:9.4f} min_nn={r['min_nn']:.4f} nn_cv={r['nn_cv']:.4f}")


if __name__ == "__main__":
    main(sys.argv[1:])
