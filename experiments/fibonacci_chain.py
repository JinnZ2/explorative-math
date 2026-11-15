"""c. fibonacci_chain -- cut-and-project from Z^2 to a 1D chain.

CC0. Co-authored. Stdlib only.

Take the square lattice Z^2 and a line through the origin with slope s.
Keep the lattice points whose perpendicular coordinate falls inside the
window (the projection of the unit square onto the perpendicular), and
project the kept points onto the line. The gaps between neighbours form the
chain.

At s = 1/phi there are two gap lengths, L and S, with L/S = phi. The L/S
counts tend to phi and the chain has no period. A rational slope (control,
s = 1/2) gives a periodic chain.

  python3 experiments/fibonacci_chain.py [seed]
"""
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
import nature_catalog  # noqa: E402

NAME = "fibonacci_chain"
PHI = (1 + 5 ** 0.5) / 2
EXTENT = 300         # lattice points with |m|, |n| <= EXTENT
OFFSET = 1e-7        # shifts the window off lattice-point boundaries
PREFIX_LENGTHS = (10, 20, 40, 80, 160, 320)


def cut_and_project(slope, extent=EXTENT, offset=OFFSET):
    nrm = math.hypot(1.0, slope)
    par = (1.0 / nrm, slope / nrm)
    perp = (-slope / nrm, 1.0 / nrm)
    corners = [c[0] * perp[0] + c[1] * perp[1] for c in ((0, 0), (1, 0), (0, 1), (1, 1))]
    lo, hi = min(corners) + offset, max(corners) + offset
    pts = []
    for m in range(-extent, extent + 1):
        for n in range(-extent, extent + 1):
            q = m * perp[0] + n * perp[1]
            if lo <= q < hi:
                pts.append(m * par[0] + n * par[1])
    pts.sort()
    # keep the stretch well inside the lattice box, where no point is cut off
    margin = extent * 0.5 * par[0]
    pts = [x for x in pts if -margin <= x <= margin]
    return [b - a for a, b in zip(pts, pts[1:])]


def to_word(gaps, sep=1e-6):
    """Cluster gap lengths (values closer than `sep` are one length). Returns the
    L/S word and the cluster means, or (None, means) when there are not exactly two."""
    clusters = []
    for g in sorted(gaps):
        if clusters and g - clusters[-1][-1] <= sep:
            clusters[-1].append(g)
        else:
            clusters.append([g])
    means = [sum(c) / len(c) for c in clusters]
    if len(means) != 2:
        return None, means
    s, l = means
    return "".join("L" if abs(g - l) < abs(g - s) else "S" for g in gaps), means


def substitution_word(n_iter=20):
    w = "L"
    for _ in range(n_iter):
        w = "".join("LS" if c == "L" else "L" for c in w)
    return w


def smallest_period(w):
    n = len(w)
    for p in range(1, n):
        if all(w[i] == w[i + p] for i in range(n - p)):
            return p
    return n


def run(seed=0):
    out = {}
    for label, slope in (("golden", 1 / PHI), ("control_rational_1_2", 0.5)):
        gaps = cut_and_project(slope)
        word, distinct = to_word(gaps)
        row = {"slope": slope, "n_gaps": len(gaps),
               "distinct_gap_lengths": distinct}
        if word is not None:
            nl, ns = word.count("L"), word.count("S")
            row.update({
                "word": word,
                "L_over_S_length": distinct[1] / distinct[0],
                "count_L": nl, "count_S": ns,
                "count_ratio": nl / ns if ns else None,
                "smallest_period": smallest_period(word),
                "prefix_smallest_periods": {str(k): smallest_period(word[:k])
                                            for k in PREFIX_LENGTHS if k <= len(word)},
            })
        else:
            row.update({"word": None, "L_over_S_length": None, "count_L": None,
                        "count_S": None, "count_ratio": None, "smallest_period": None,
                        "prefix_smallest_periods": None})
        out[label] = row
    sub = substitution_word()
    out["golden"]["factor_of_substitution_word"] = (
        out["golden"]["word"] in sub if out["golden"]["word"] else None)
    return out


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    g = out["golden"]
    comp = emlib.compressibility(list(g["word"]), "golden-slope chain as L/S letters")
    per = g["prefix_smallest_periods"]
    cper = out["control_rational_1_2"]["prefix_smallest_periods"]
    ai = (f"OBSERVED: two gap lengths, L/S = {g['L_over_S_length']:.12f} (phi = {PHI:.12f}); "
          f"counts L={g['count_L']} S={g['count_S']}, ratio {g['count_ratio']:.6f}; "
          f"factor of the substitution word L->LS, S->L: {g['factor_of_substitution_word']}. "
          f"Smallest period of each prefix (length: period) {per}; control slope 1/2 "
          f"{cper}. DERIVED: a finite word always has some period, so aperiodicity shows "
          "only as the smallest period growing with length; it does for the golden chain "
          "and stays fixed for the control. AI reading: which numbers the prefix periods "
          "land on is the part not asked for; not interpretable atm beyond listing them.")
    p = emlib.write_run(NAME, seed,
                        {"slope": "1/phi", "control_slope": str(Fraction(1, 2)),
                         "extent": EXTENT, "window_offset": OFFSET},
                        out, None, ai, comp,
                        nature_catalog.match("fibonacci_chain_quasicrystal"))
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
