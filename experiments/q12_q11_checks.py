"""Q12 (extended) and Q11 (aligned re-checks) for JinnZ2/explorative-math.
CC0. Co-authored. Stdlib only. Imports the repo read-only; patches nothing on disk (one in-memory
patch of explorative_math.geometry._rng, stated in the Q11 orthogonality record).
Predictions: docs/prereg/prereg_q12_q11_q5_q9_2026-11-18.json (sha256 6ae97bcc...606c), written before this ran.
Deterministic, runs in about a second.
  python3 experiments/q12_q11_checks.py
"""
import math
import random
import sys
from decimal import Decimal, getcontext
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = None                     # None = the repo runs/ directory
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "experiments"))
import emlib  # noqa: E402
from toolkit import loghook  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402
from explorative_math import geometry as A_geo  # noqa: E402
from explorative_math import projections as A_proj  # noqa: E402
import shell_concentration as B_shell  # noqa: E402

getcontext().prec = 80
PHI = P.PHI


def divisors_upto(k, kmax):
    return [d for d in range(2, kmax + 1) if k % d == 0]


# ---------------------------------------------------------------- Q12 extended
def q12():
    rows, union_rows = [], []
    ks = (5, 8, 13, 16, 24)
    kmaxes = (12, 24, 32)
    refuted = []
    for k in ks:
        poly = P.regular_polygon(k)
        for kmax in kmaxes:
            got = Q.rotational_symmetry(poly, kmax, 1e-6)
            exp = divisors_upto(k, kmax)
            rows.append({"k": k, "k_max": kmax, "got": got, "divisors_le_kmax": exp,
                         "equals_divisor_set": got == exp,
                         "bound": Q.rotational_symmetry_bound(poly, kmax, 1e-6)})
            if got != exp:
                refuted.append(("kgon_not_divisor_set", k, kmax, got, exp))
        for label, deg in (("rot_7deg", 7.0), ("rot_180_over_k", 180.0 / k)):
            th = math.radians(deg)
            c, s = math.cos(th), math.sin(th)
            scaled = [[PHI * (c * x - s * y), PHI * (s * x + c * y)] for x, y in poly]
            union = poly + scaled
            for kmax in kmaxes:
                got = Q.rotational_symmetry(union, kmax, 1e-6)
                base = Q.rotational_symmetry(poly, kmax, 1e-6)
                union_rows.append({"k": k, "rotation": label, "rotation_deg": deg, "k_max": kmax,
                                   "union_points": len(union), "union_list": got,
                                   "bare_kgon_list": base, "same_list": got == base,
                                   "k_max_ge_k": kmax >= k})
                if kmax >= k and got != base:
                    refuted.append(("union_differs", k, label, kmax, got, base))
    outputs = {"kgons": rows, "unions": union_rows,
               "prediction_refuted_by": [list(map(str, r)) for r in refuted],
               "prediction_status": "NONCONFORMS" if refuted else "CONFORMS"}
    params = {"experiment": "q12_extended", "ks": list(ks), "k_maxes": list(kmaxes),
              "rel_tol": 1e-6, "prereg": "docs/prereg/prereg_q12_q11_q5_q9_2026-11-18.json",
              "repo_commit": "3162c9950136a4a7f14e1e064f8d335db33968c8"}
    return loghook.log_run("q12_extended", 0, params, outputs, out_dir=OUT), outputs


# ---------------------------------------------------------------- Q11 aligned
def ulp_err(x, truth):
    return float((Decimal(x) - truth) / Decimal(math.ulp(x)))


def q11_constants():
    sqrt5 = Decimal(5).sqrt()
    phi = (1 + sqrt5) / 2
    truth = Decimal(360) / (phi * phi)
    a = A_proj.GOLDEN_ANGLE_DEGREES
    b = 360.0 * (1.0 - 1.0 / P.PHI)
    return {"truth_60_digits": str(truth)[:62],
            "A_value": a, "A_error_ulp": ulp_err(a, truth),
            "B_value": b, "B_error_ulp": ulp_err(b, truth),
            "A_minus_B": a - b,
            "aligned_diff_when_B_uses_A_formula": a - 360.0 / (P.PHI ** 2),
            "aligned_diff_when_A_uses_B_formula": 360.0 * (1.0 - 1.0 / A_proj.PHI) - b,
            "within_2_ulp_both": abs(ulp_err(a, truth)) <= 2 and abs(ulp_err(b, truth)) <= 2}


def q11_shell():
    cells, a_closer, b_closer, tie, refuted = [], 0, 0, 0, []
    for d in (1, 2, 3, 5, 10, 20, 50, 100, 500):
        for t in (0.01, 0.05, 0.1, 0.5):
            a = A_geo.shell_fraction(d, t)
            b = B_shell.exact(d, t)
            truth = 1 - (1 - Decimal(t)) ** d
            ea, eb = abs(Decimal(a) - truth), abs(Decimal(b) - truth)
            differ = a != b
            who = "tie" if ea == eb else ("A" if ea < eb else "B")
            if differ:
                a_closer += who == "A"
                b_closer += who == "B"
                tie += who == "tie"
                if who == "B":
                    refuted.append((d, t))
            cells.append({"d": d, "t": t, "A": a, "B": b, "differ": differ, "closer_to_truth": who,
                          "A_err": float(ea), "B_err": float(eb)})
    n_differ = sum(c["differ"] for c in cells)
    return {"cells": cells, "n_differ": n_differ, "A_closer": a_closer, "B_closer": b_closer,
            "tie": tie, "cells_where_B_closer": refuted,
            "aligned_diff_when_both_use_A_formula": 0.0}


def q11_polytopes():
    out = []
    shapes = {"simplex": ("simplex", P.simplex), "hypercube": ("hypercube", P.hypercube),
              "crosspolytope": ("crosspolytope", P.cross_polytope)}
    for name, (a_name, b_fn) in shapes.items():
        for n in (3, 4, 5):
            va, vb = A_proj._vertices(a_name, n), b_fn(n)
            ra = max(P.norm(v) for v in va)
            rb = max(P.norm(v) for v in vb)
            same_dim = len(va[0]) == len(vb[0])
            rec = {"shape": name, "n": n, "A_ambient": len(va[0]), "B_ambient": len(vb[0]),
                   "A_circumradius": ra, "B_circumradius": rb, "same_ambient": same_dim}
            if same_dim and len(va) == len(vb):
                na = sorted(tuple(round(x / ra, 12) for x in v) for v in va)
                nb = sorted(tuple(round(x / rb, 12) for x in v) for v in vb)
                rec["identical_after_unit_circumradius_and_sort"] = na == nb
                rec["max_abs_coord_diff"] = max(abs(x - y) for u, w in zip(na, nb) for x, y in zip(u, w))
            else:
                rec["identical_after_unit_circumradius_and_sort"] = None
                rec["note"] = "ambient dimension or vertex count differs; coordinates not comparable"
            out.append(rec)
    return out


def match_error(a_pts, b_pts):
    """max over a-points of the distance to the nearest b-point, and the reverse."""
    def one_way(x, y):
        return max(min(math.dist(p, q) for q in y) for p in x)
    return max(one_way(a_pts, b_pts), one_way(b_pts, a_pts))


def q11_phi():
    golden = A_proj.GOLDEN_ANGLE_DEGREES
    rows = []
    for name in ("simplex", "crosspolytope"):
        run = next(A_proj.iter_projection_runs(
            {"dimension": 4, "shape": name, "frames": 4, "seed": 0, "composition": "phi"}))
        for frame in [o for o in run["outputs"] if o["kind"] == "projection_frame"]:
            pts = frame["projected_points"]
            a_union = frame["composed_points"]
            rec = {"shape": name, "frame": frame["frame_index"], "A_points": len(a_union)}
            for sign in (+1, -1):
                th = math.radians(sign * golden)
                c, s = math.cos(th), math.sin(th)
                scaled = [[PHI * (c * x - s * y), PHI * (s * x + c * y)] for x, y in pts]
                b_union = [list(p) for p in pts] + scaled
                key = "plus_golden" if sign > 0 else "minus_golden"
                rec[key] = {"B_points": len(b_union),
                            "max_nearest_point_distance": match_error(a_union, b_union),
                            "orderwise_max_abs_diff": (max(abs(x - y) for p, q in zip(a_union, b_union)
                                                          for x, y in zip(p, q))
                                                      if len(a_union) == len(b_union) else None)}
            rec["matches_within_1e-9"] = (rec["plus_golden"]["max_nearest_point_distance"] < 1e-9 or
                                          rec["minus_golden"]["max_nearest_point_distance"] < 1e-9)
            rows.append(rec)
    return rows


def q11_orthogonality():
    """Align the stream choice: patch A's _rng (in memory) to Random(0) and feed B the same stream."""
    import hashlib  # noqa: F401
    orig = A_geo._rng
    rows = []
    try:
        for d in (2, 20, 50, 100, 200):
            A_geo._rng = lambda params, _s=0: random.Random(_s)
            cfg = A_geo.validate_geometry_config({"dimensions": [d], "samples": 400, "seed": 0,
                                                  "experiment": "orthogonality"})
            run_a = A_geo._orthogonality(d, cfg)
            pairs_a = [o for o in run_a["outputs"] if o["kind"] == "pair"]
            summ = [o for o in run_a["outputs"] if o["kind"] == "summary"][0]
            rng = random.Random(0)
            pairs_b = [(P.random_unit(rng, d), P.random_unit(rng, d)) for _ in range(400)]
            same_vectors = sum(1 for pa, (u, v) in zip(pairs_a, pairs_b)
                               if pa["left"] == u and pa["right"] == v)
            ang_b = [P.angle_deg(u, v) for u, v in pairs_b]
            rows.append({"d": d, "pairs": 400, "pairs_with_identical_vectors": same_vectors,
                         "A_mean_angle": summ["mean_angle_degrees"],
                         "B_mean_angle": math.fsum(ang_b) / 400,
                         "diff": summ["mean_angle_degrees"] - math.fsum(ang_b) / 400,
                         "A_mean_squared_dot": summ["mean_squared_dot"],
                         "B_mean_squared_dot_from_same_pairs": math.fsum(P.dot(u, v) ** 2 for u, v in pairs_b) / 400})
    finally:
        A_geo._rng = orig
    return rows


def q11():
    consts, shell, poly, phi, orth = (q11_constants(), q11_shell(), q11_polytopes(),
                                      q11_phi(), q11_orthogonality())
    checks = {
        "constants_within_2ulp_and_aligned_zero": consts["within_2_ulp_both"] and
            consts["aligned_diff_when_B_uses_A_formula"] == 0.0,
        "shell_A_at_least_as_close_in_all_DIFFER_cells": shell["B_closer"] == 0,
        "polytopes_hypercube_cross_identical_after_alignment": all(
            r["identical_after_unit_circumradius_and_sort"] is True for r in poly if r["shape"] != "simplex"),
        "polytopes_simplex_remaining_difference": sorted({(r["A_ambient"], r["B_ambient"]) for r in poly
                                                           if r["shape"] == "simplex"}),
        "phi_all_8_frames_match_with_one_sign": all(r["matches_within_1e-9"] for r in phi),
    }
    outputs = {"constants": consts, "shell": shell, "polytopes": poly, "phi_composition": phi,
               "orthogonality_aligned_stream": orth, "prediction_checks": checks,
               "orthogonality_note": "A's _rng replaced in memory by Random(0) for this call only; "
                                     "no repo file changed. Compare with compare/runs/orthogonality_stats.jsonl (5 DIFFER)."}
    params = {"experiment": "q11_aligned", "groups": ["constants", "shell", "polytopes", "phi", "orthogonality"],
              "prereg": "docs/prereg/prereg_q12_q11_q5_q9_2026-11-18.json",
              "repo_commit": "3162c9950136a4a7f14e1e064f8d335db33968c8"}
    return loghook.log_run("q11_aligned", 0, params, outputs, out_dir=OUT), outputs


if __name__ == "__main__":
    p12, o12 = q12()
    print("Q12", p12, o12["prediction_status"], "refuted_by:", o12["prediction_refuted_by"][:5])
    for r in o12["kgons"]:
        print("  kgon", r["k"], "kmax", r["k_max"], r["got"])
    print("  unions same_list when k_max>=k:", all(u["same_list"] for u in o12["unions"] if u["k_max_ge_k"]),
          "| unions that differ from bare k-gon anywhere:",
          [(u["k"], u["rotation"], u["k_max"], u["union_list"], u["bare_kgon_list"]) for u in o12["unions"] if not u["same_list"]])
    p11, o11 = q11()
    print("Q11", p11)
    print(" checks:", o11["prediction_checks"])
    c = o11["constants"]
    print(" constants: A err %.3f ulp, B err %.3f ulp, A-B=%.3g" % (c["A_error_ulp"], c["B_error_ulp"], c["A_minus_B"]))
    s = o11["shell"]
    print(" shell: differ", s["n_differ"], "A closer", s["A_closer"], "B closer", s["B_closer"], "tie", s["tie"], s["cells_where_B_closer"])
    for r in o11["polytopes"]:
        print("  poly", r["shape"], r["n"], r["A_ambient"], r["B_ambient"], r["identical_after_unit_circumradius_and_sort"])
    for r in o11["phi_composition"]:
        print("  phi", r["shape"], r["frame"], "+%.2e" % r["plus_golden"]["max_nearest_point_distance"],
              "-%.2e" % r["minus_golden"]["max_nearest_point_distance"], r["matches_within_1e-9"])
    for r in o11["orthogonality_aligned_stream"]:
        print("  orth d", r["d"], "identical pairs", r["pairs_with_identical_vectors"], "/400", "angle diff %.3g" % r["diff"])
