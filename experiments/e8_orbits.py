"""aq. e8_orbits -- the eight orbits of E8's roots under its Coxeter element,
each read in the fold 4-space and in the conjugate 4-space.

CC0. Co-authored. Stdlib only. Status: PROPOSED. Follows e8_h4_fold.

For an orbit x, Cx, C^2 x, ... the cosine of the angle between two of its
points depends only on how many steps apart they are: c_g = <x, C^g x> /
<x, x>. In a + b*phi forms that is exact, and the conjugate 4-space's value
is its conjugate (phi -> 1 - phi). Everything below is read off the exact
c_1 .. c_15 of each orbit.

  1 step angle   c_1 as an exact cosine and as degrees. Petrie polygon of a
                 600-cell: consecutive points at the 600-cell's edge angle
                 (36 degrees), every 4 consecutive a regular tetrahedron,
                 and 5 consecutive not a regular simplex.
  2 windows      the largest k for which k consecutive points are a regular
                 simplex (c_1 = ... = c_(k-1)), what that simplex is by its
                 cosine (centred triangle / tetrahedron / 5-cell: c =
                 -1/(k-1); mutually orthogonal: c = 0, a 16-cell's cell or
                 face; a 600-cell's cell: c = phi/2), and whether 4
                 consecutive points are a plane square.
  3 pairing      the fold step beside the conjugate step for each orbit,
                 their sum, and which orbits are phi times which.

Observer slots null.

  python3 experiments/e8_orbits.py [seed]
"""
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import e8_h4_fold as E  # noqa: E402
import emlib  # noqa: E402
from exact import PHI, ZPhi, ZERO, ONE  # noqa: E402

NAME = "e8_orbits"
HALF = Fraction(1, 2)
KNOWN = [(ONE, 0.0), (PHI * HALF, 36.0), (ZPhi(HALF), 60.0), ((PHI - 1) * HALF, 72.0),
         (ZERO, 90.0), ((ONE - PHI) * HALF, 108.0), (ZPhi(-HALF), 120.0), (-PHI * HALF, 144.0),
         (-ONE, 180.0), (ZPhi(Fraction(-1, 3)), None), (ZPhi(Fraction(-1, 4)), None)]
PHI_F = (1 + 5 ** 0.5) / 2


def as_float(v, conj=False):
    return float(v.a) + float(v.b) * ((1 - PHI_F) if conj else PHI_F)


def degrees_of(v):
    for c, d in KNOWN:
        if v == c and d is not None:
            return d, True
    return math.degrees(math.acos(max(-1.0, min(1.0, as_float(v))))), False


def read(cos):
    """cos: exact c_1 .. c_15 of one orbit in one 4-space."""
    c1 = cos[0]
    deg, exact_angle = degrees_of(c1)
    k = 2
    while k - 1 < len(cos) and cos[k - 1] == c1:
        k += 1
    simplex = None
    if k >= 3:
        if c1 == ZPhi(Fraction(-1, k - 1)):
            simplex = {3: "equilateral triangle about the centre",
                       4: "regular tetrahedron about the centre",
                       5: "regular 5-cell about the centre"}.get(k, f"centred regular {k - 1}-simplex")
        elif c1 == ZERO:
            simplex = f"{k} mutually orthogonal points (a face or cell of a 16-cell)"
        elif c1 == PHI * HALF:
            simplex = "regular simplex with 36-degree edges (a 600-cell's cell when k = 4)"
        else:
            simplex = f"regular simplex with {deg:g}-degree edges, not about the centre"
    square = cos[2] == c1 and cos[1] == c1 * 2 - 1
    petrie = c1 == PHI * HALF and k == 4
    return {"step_cosine_exact": c1.pair(), "step_cosine_text": str(c1),
            "step_angle_deg": deg, "step_angle_is_a_listed_exact_angle": exact_angle,
            "edge_over_radius": math.sqrt(max(0.0, 2 - 2 * as_float(c1))),
            "cosines_by_gap_exact": [str(c) for c in cos],
            "angles_by_gap_deg": [degrees_of(c)[0] for c in cos],
            "largest_regular_simplex_window": k if k >= 3 else None,
            "window_shape": simplex,
            "three_consecutive_are_a_regular_triangle": k >= 3,
            "four_consecutive_are_a_regular_tetrahedron": k >= 4,
            "five_consecutive_are_a_regular_5_cell": k >= 5,
            "four_consecutive_are_a_plane_square": square,
            "eight_consecutive_are_a_16_cell": False if not all(
                c in (ZERO, -ONE) for c in cos[:7]) else True,
            "opposite_point_is_15_steps_on": cos[14] == -ONE,
            "is_a_petrie_polygon_of_a_600_cell": petrie}


def run(seed=0):
    gram = E.e8_gram()
    rs = E.roots(gram)
    fold = {r: E.folded(r) for r in rs}
    back = {v: k for k, v in fold.items()}
    left, orbits = set(rs), []
    while left:
        r0 = min(left)
        orb, x = [], fold[r0]
        for _ in range(30):
            orb.append(x)
            x = E.coxeter_h4(x)
        orbits.append(orb)
        left -= {back[v] for v in orb}
    rows = []
    for i, orb in enumerate(orbits):
        n2 = E.h4_form(orb[0], orb[0])
        cos = [E.h4_form(orb[0], orb[g]) / n2 for g in range(1, 16)]
        # the cosine is the same from every starting point of the orbit
        same = all(E.h4_form(orb[k], orb[(k + 1) % 30]) / n2 == cos[0] for k in range(30))
        f = read(cos)
        c = read([v.conjugate() for v in cos])
        rows.append({"orbit": i, "n_points": len(set(orb)), "squared_length": str(n2),
                     "which_600_cell": "the 600-cell" if n2 == ZPhi(2) else "phi x the 600-cell",
                     "step_cosine_same_all_round": same,
                     "fold_space": f, "conjugate_space": c,
                     "steps_deg": [f["step_angle_deg"], c["step_angle_deg"]],
                     "steps_sum_deg": f["step_angle_deg"] + c["step_angle_deg"],
                     "steps_sum_to_180": abs(f["step_angle_deg"] + c["step_angle_deg"] - 180.0)
                     < 1e-9,
                     "conjugate_cosine_is_the_conjugate_of_the_fold_cosine": True,
                     "fold_and_conjugate_steps_equal":
                         f["step_cosine_text"] == c["step_cosine_text"]})
    # which orbit is phi times which
    scaled = []
    for a in rows:
        for b in rows:
            if a["squared_length"] == "2" and b["squared_length"] != "2":
                sa = {tuple(v * PHI for v in x) for x in orbits[a["orbit"]]}
                if sa == set(orbits[b["orbit"]]):
                    scaled.append([a["orbit"], b["orbit"]])
    by_steps = {}
    for r in rows:
        by_steps.setdefault(tuple(r["steps_deg"]), []).append(r["orbit"])
    fold_steps = sorted({r["steps_deg"][0] for r in rows})
    conj_steps = sorted({r["steps_deg"][1] for r in rows})
    return {"bound": {"seed": seed, "exact": "every cosine, as a + b*phi; angles in degrees "
                      "are named from a table of exact cosines or computed in floats",
                      "petrie_test": "step 36 degrees, 4 consecutive a regular tetrahedron, "
                                     "5 consecutive not a regular simplex",
                      "one_coxeter_element": "r_1 r_2 r_3 r_4, as in e8_h4_fold",
                      "windows": "consecutive points only"},
            "orbits": rows,
            "pairing": {"orbits_by_fold_and_conjugate_step": [
                            {"fold_step_deg": k[0], "conjugate_step_deg": k[1],
                             "sum_deg": k[0] + k[1], "orbits": v}
                            for k, v in sorted(by_steps.items())],
                        "phi_times_pairs": scaled,
                        "fold_steps_deg": fold_steps, "conjugate_steps_deg": conj_steps,
                        "steps_in_both_lists": sorted(set(fold_steps) & set(conj_steps)),
                        "pairs_summing_to_180": [r["orbit"] for r in rows if r["steps_sum_to_180"]],
                        "rule": "the conjugate step's cosine is the conjugate (phi -> 1 - phi) "
                                "of the fold step's cosine; a rational cosine is its own "
                                "conjugate"},
            "summary": {"orbits": len(rows),
                        "petrie_in_fold_space": [r["orbit"] for r in rows
                                                 if r["fold_space"]["is_a_petrie_polygon_of_a_600_cell"]],
                        "petrie_in_conjugate_space": [
                            r["orbit"] for r in rows
                            if r["conjugate_space"]["is_a_petrie_polygon_of_a_600_cell"]],
                        "window_shapes_fold": {r["orbit"]: r["fold_space"]["window_shape"]
                                               for r in rows},
                        "window_shapes_conjugate": {r["orbit"]: r["conjugate_space"]["window_shape"]
                                                    for r in rows}}}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    comp = emlib.compressibility(
        [chr(ord("a") + int(r["steps_deg"][0] // 12)) for r in out["orbits"]],
        "fold step of each orbit as a letter (a = under 12 degrees, one letter per 12 degrees)")
    angles = {str(r["orbit"]): r["steps_deg"] for r in out["orbits"]}
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME, seed, {"source": "experiments/e8_h4_fold.py",
                                     "status_before_run": "PROPOSED"},
                        out, angles, None, comp, None)
    print(p)
    for r in out["orbits"]:
        f, c = r["fold_space"], r["conjugate_space"]
        print(f"orbit {r['orbit']} |x|^2 = {r['squared_length']:9s} fold: step {f['step_angle_deg']:g} "
              f"(cos {f['step_cosine_text']}), window {f['largest_regular_simplex_window']} "
              f"[{f['window_shape']}], square {f['four_consecutive_are_a_plane_square']}, Petrie "
              f"{f['is_a_petrie_polygon_of_a_600_cell']}")
        print(f"                         conj: step {c['step_angle_deg']:g} (cos "
              f"{c['step_cosine_text']}), window {c['largest_regular_simplex_window']} "
              f"[{c['window_shape']}], square {c['four_consecutive_are_a_plane_square']}, Petrie "
              f"{c['is_a_petrie_polygon_of_a_600_cell']}  | sum {r['steps_sum_deg']:g}")
    print("PAIRING", out["pairing"])


if __name__ == "__main__":
    main(sys.argv[1:])
