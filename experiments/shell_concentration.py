"""b. shell_concentration -- fraction of a d-ball's volume within eps of the surface.

CC0. Co-authored. Stdlib only.

Exact: fraction = 1 - (1 - eps)^d. Checked by rejection sampling from the
cube for small d (the acceptance rate falls fast with d, so the Monte Carlo
check stops at MC_MAX_DIM; past that only the exact value is recorded).

  python3 experiments/shell_concentration.py [seed]
"""
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402

NAME = "shell_concentration"
DIMS = (1, 2, 3, 5, 10, 20, 50, 100, 500)
EPS = (0.01, 0.05, 0.1)
MC_MAX_DIM = 6
MC_ACCEPTED = 4000


def exact(d, eps):
    return 1.0 - (1.0 - eps) ** d


def log_ball_to_cube(d):
    """log of (unit d-ball volume) / (volume of the cube [-1,1]^d), via lgamma."""
    return (d / 2) * math.log(math.pi) - math.lgamma(d / 2 + 1) - d * math.log(2)


def monte_carlo(rng, d, eps, accepted=MC_ACCEPTED):
    got, inside, tries = 0, 0, 0
    while got < accepted:
        tries += 1
        x = [rng.uniform(-1.0, 1.0) for _ in range(d)]
        r = emlib.norm(x)
        if r <= 1.0:
            got += 1
            if r >= 1.0 - eps:
                inside += 1
    return inside / accepted, accepted / tries


def run(seed=0):
    rng = random.Random(seed)
    rows = []
    for d in DIMS:
        for e in EPS:
            row = {"d": d, "eps": e, "exact": exact(d, e), "monte_carlo": None,
                   "acceptance_rate": None, "log10_ball_to_cube_volume": log_ball_to_cube(d) / math.log(10)}
            if d <= MC_MAX_DIM:
                mc, acc = monte_carlo(rng, d, e)
                row["monte_carlo"], row["acceptance_rate"] = mc, acc
                ex = row["exact"]
                row["z"] = (mc - ex) / math.sqrt(ex * (1 - ex) / MC_ACCEPTED)
            else:
                row["z"] = None
            rows.append(row)
    return rows


def main(argv):
    seed = int(argv[0]) if argv else 0
    rows = run(seed)
    sym = ["0123456789"[min(9, int(r["exact"] * 10))] for r in rows]
    comp = emlib.compressibility(sym, "exact shell fraction in tenths, rows in (d, eps) order")
    half = {e: next((r["d"] for r in rows if r["eps"] == e and r["exact"] >= 0.5), None) for e in EPS}
    zmax = max(abs(r["z"]) for r in rows if r["z"] is not None)
    ai = (f"OBSERVED: Monte Carlo vs exact where both exist (d<=6): max |z| = {zmax:.2f} "
          f"over {sum(1 for r in rows if r['z'] is not None)} rows at {MC_ACCEPTED} samples; the shell "
          f"fraction passes 0.5 at d={half} for eps={list(EPS)}. Monte Carlo is null above "
          f"d={MC_MAX_DIM} because cube rejection sampling stops working there "
          "(acceptance = ball/cube volume, recorded as log10). AI reading: the acceptance "
          "collapse is the same fact as the shell result seen from outside the ball.")
    p = emlib.write_run(NAME, seed,
                        {"dims": list(DIMS), "eps": list(EPS), "mc_max_dim": MC_MAX_DIM,
                         "mc_accepted": MC_ACCEPTED},
                        {"rows": rows}, None, ai, comp, None)
    print(p)
    for r in rows:
        mc = "  void" if r["monte_carlo"] is None else f"{r['monte_carlo']:.4f}"
        print(f"d={r['d']:4d} eps={r['eps']:.2f} exact={r['exact']:.4f} mc={mc}")


if __name__ == "__main__":
    main(sys.argv[1:])
