"""ai. ten_vertex_volume -- volume of the ten-direction hull against a stated bound.

CC0. Co-authored. Stdlib only. Status: PROPOSED.

Stated bound (REFERENCES/related_studies.md R4, abstract read): a
4-polytope with ten vertices inscribed in the unit sphere S^3 has volume at
most 25/384 * (10 + 2 sqrt 5), about 0.9420, with equality only for the
convex hull of two regular pentagons in orthogonal planes.

Computed here, in exact a + b*phi forms (sqrt 5 = 2 phi - 1):

  A  the hull of a regular 5-cell and its negative (the ten directions of
     sphere_compounds), on the unit sphere
  B  the hull of two regular pentagons in orthogonal planes through the
     origin, on the unit sphere

Each volume is the sum over tetrahedral facets of |det| / 24; the facets
are found in floats (toolkit/hull4.py) and the determinants taken exactly.
Reported: match or miss against the bound. Nothing is interpreted;
observer slots null.

  python3 experiments/ten_vertex_volume.py [seed]
"""
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import complement_fit_v2 as V  # noqa: E402
import emlib  # noqa: E402
from exact import PHI, ZPhi, ZERO  # noqa: E402
from exact import inflation as XI  # noqa: E402
from exact import volume as XV  # noqa: E402
from toolkit import complement as CP  # noqa: E402
from toolkit import hull4 as H4  # noqa: E402
from toolkit import primitives as P  # noqa: E402

NAME = "ten_vertex_volume"
PHI_F = (1 + 5 ** 0.5) / 2
BOUND = ZPhi(Fraction(25, 48), Fraction(25, 96))      # 25/384 (10 + 2 sqrt 5), sqrt 5 = 2 phi - 1
STATED = "25/384 * (10 + 2 sqrt 5)"


def as_float(x):
    return float(x.a) + float(x.b) * PHI_F


def entry(name, exact_points, float_points, scale=None):
    """scale: an exact factor the determinant must be multiplied by, when the
    exact coordinates are coefficients of a fixed irrational (see pentagons)."""
    fs = H4.facets(float_points)
    vol = XV.volume_from_facets(exact_points, fs)
    if scale is not None:
        vol = vol * scale
    fvol = sum(abs(_det_f([float_points[i] for i in f])) for f in fs) / 24.0
    diff = vol - BOUND
    return {"name": name, "n_points": len(float_points),
            "on_unit_sphere_max_deviation": max(abs(P.norm(p) - 1.0) for p in float_points),
            "n_facets": len(fs), "facets_all_tetrahedra": all(len(f) == 4 for f in fs),
            "facets_close_up": H4.closed(fs),
            "volume_exact": vol.pair(), "volume_exact_text": str(vol),
            "volume_float_from_exact": as_float(vol), "volume_float_direct": fvol,
            "bound_exact": BOUND.pair(), "bound_float": as_float(BOUND),
            "volume_minus_bound_exact": diff.pair(),
            "volume_minus_bound_float": as_float(diff),
            "volume_over_bound": as_float(vol) / as_float(BOUND),
            "equals_the_bound_exactly": vol == BOUND,
            "below_the_bound": vol < BOUND, "above_the_bound": vol > BOUND,
            "verdict": "MATCH" if vol == BOUND else "MISS",
            "distinct_pairwise_angles_deg":
                CP.distinct_values(V.signature(CP.unit_dirs(float_points)), tol=0.05)}


def _det_f(rows):
    n = len(rows)
    if n == 1:
        return rows[0][0]
    return sum((-1) ** j * rows[0][j] * _det_f([r[:j] + r[j + 1:] for r in rows[1:]])
               for j in range(n))


def run(seed=0):
    five = XV.five_cell_unit()
    ten = five + XV.negated(five)
    ten_f = [[as_float(x) for x in v] for v in ten]
    a = entry("A: 5-cell and its negative", ten, ten_f)
    a["same_set_as_the_sphere_compounds_entry"] = V.same_signature(
        V.signature(ten_f), V.signature(V.catalogue()["5-cell + its negative (as built)"]),
        tol=1e-6)
    a["five_cell_is_regular_and_unit_exactly"] = (
        all(XV.dot(v, v) == ZPhi(1, 0) for v in five)
        and {XV.dot(five[i], five[j]) for i in range(5) for j in range(i + 1, 5)}
        == {ZPhi(Fraction(-1, 4), 0)})
    # two pentagons in orthogonal planes: (zeta^k, 0) and (0, zeta^k). A plane
    # point u + v*w has coordinates (u, v*s), s = sqrt(2 + phi); the exact rows
    # hold (u, v) and the determinant picks up s^2 = 2 + phi (one s per plane).
    pent = XI.pentagon()
    rows = [(z.u, z.v, ZERO, ZERO) for z in pent] + [(ZERO, ZERO, z.u, z.v) for z in pent]
    s = (2 + PHI_F) ** 0.5
    rows_f = [[as_float(r[0]), as_float(r[1]) * s, as_float(r[2]), as_float(r[3]) * s]
              for r in rows]
    b = entry("B: two regular pentagons in orthogonal planes", rows, rows_f,
              scale=ZPhi(2, 0) + PHI)
    b["exact_note"] = ("coordinates are (u, v*s) per plane with s = sqrt(2 + phi); the exact "
                       "determinant of the (u, v) rows is multiplied by s^2 = 2 + phi")
    return {"stated_bound": STATED, "bound_exact_form": "25/48 + 25/96 * phi",
            "source": "REFERENCES/related_studies.md R4 (abstract read)",
            "bound": {"seed": seed, "arithmetic": "determinants in exact a + b*phi forms; "
                      "facets found in floats (tol 1e-7) and checked to close up",
                      "only_two_point_sets": "this run does not search ten-point sets; it "
                                             "reads A and B"},
            "A": a, "B": b,
            "A_over_B": a["volume_float_from_exact"] / b["volume_float_from_exact"]}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    comp = emlib.compressibility([out["A"]["verdict"][1], out["B"]["verdict"][1]],
                                 "second letter of the verdict for A then B (A match, I miss)")
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME, seed, {"stated_bound": STATED, "status_before_run": "PROPOSED",
                                     "point_sets": ["A: 5-cell and its negative",
                                                    "B: two pentagons in orthogonal planes"]},
                        out, None, None, comp, None)
    print(p)
    for k in ("A", "B"):
        e = out[k]
        print(f"{e['name']}: volume {e['volume_exact_text']} = {e['volume_float_from_exact']:.10f}"
              f" (direct float {e['volume_float_direct']:.10f}); bound {e['bound_float']:.10f}; "
              f"{e['verdict']}; volume/bound {e['volume_over_bound']:.6f}; facets "
              f"{e['n_facets']}, close up {e['facets_close_up']}")
    print("A over B", out["A_over_B"])


if __name__ == "__main__":
    main(sys.argv[1:])
