"""k. symmetry_resolution -- the rotational-symmetry probe at four bounds.

CC0. Co-authored. Stdlib only.

open_questions.md Q12: the probe returns the same list for two different
unions. Is that a real equality, or the probe's bound?

  known    two shapes whose rotational orders are known and different: the
           regular 8-gon (order 8) and the regular 16-gon (order 16). The
           probe is read at k_max 12 and 24, relative tolerance 1e-6 and
           1e-9, and each bound logs whether the two outputs separate.
  unions   the two unions Q12 is about, remade the way compare/cross.py
           makes them (Build A's phi_scale_rotate, Build B's phi_scale, on
           Build A's projected frames), read at the same four bounds.

The second part imports Build A (explorative_math) read-only, as compare/
does. Neither build is changed.

  python3 experiments/symmetry_resolution.py [seed]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from explorative_math import projections as A_proj  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "symmetry_resolution"
K_MAXES = (12, 24)
REL_TOLS = (1e-6, 1e-9)
KNOWN = (("regular_8_gon", 8), ("regular_16_gon", 16))


def grid(a, b):
    """The probe on two point sets at every (k_max, rel_tol)."""
    rows = []
    for k_max in K_MAXES:
        for tol in REL_TOLS:
            ra = Q.rotational_symmetry(a, k_max, tol)
            rb = Q.rotational_symmetry(b, k_max, tol)
            rows.append({"k_max": k_max, "rel_tol": tol, "a": ra, "b": rb,
                         "separate": ra != rb,
                         "a_bound": Q.rotational_symmetry_bound(a, k_max, tol),
                         "b_bound": Q.rotational_symmetry_bound(b, k_max, tol)})
    return rows


def known():
    (na, ka), (nb, kb) = KNOWN
    a, b = P.regular_polygon(ka), P.regular_polygon(kb)
    return {"a": {"name": na, "known_order": ka, "points": a},
            "b": {"name": nb, "known_order": kb, "points": b},
            "rows": grid(a, b)}


def unions():
    out = []
    for name in ("simplex", "crosspolytope"):
        run = next(A_proj.iter_projection_runs(
            {"dimension": 4, "shape": name, "frames": 4, "seed": 0, "composition": "phi"}))
        for frame in [o for o in run["outputs"] if o["kind"] == "projection_frame"]:
            a_union = frame["composed_points"]
            b_union = [r["point"] for r in C.phi_stack(frame["projected_points"], (0, 1))]
            dirs = [len({round(x, 9) for x in Q.polar_angles(u)}) for u in (a_union, b_union)]
            out.append({"shape": name, "dimension": 4, "frame": frame["frame_index"],
                        "a_operation": "phi_scale_rotate (Build A)",
                        "b_operation": "phi_scale (Build B)",
                        "n_points": [len(a_union), len(b_union)],
                        "distinct_directions": dirs,
                        "rows": grid(a_union, b_union)})
    return out


def run(seed=0):
    k, u = known(), unions()
    return {"known": k, "unions": u,
            "known_separate_by_bound": [[r["k_max"], r["rel_tol"], r["separate"]]
                                        for r in k["rows"]],
            "union_frames_separating_by_bound": [
                [k_max, tol, sum(1 for f in u for r in f["rows"]
                                 if r["k_max"] == k_max and r["rel_tol"] == tol
                                 and r["separate"]), len(u)]
                for k_max in K_MAXES for tol in REL_TOLS]}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    symbols = ["s" if r["separate"] else "e" for r in out["known"]["rows"]]
    symbols += ["s" if r["separate"] else "e" for f in out["unions"] for r in f["rows"]]
    comp = emlib.compressibility(
        symbols, "one letter per (pair, bound): s outputs separate, e outputs equal; "
                 "known pair first, then union frames, bounds in grid order")
    kn = "; ".join(f"k_max {r['k_max']} tol {r['rel_tol']:g}: {r['a']} vs {r['b']} -> "
                   f"{'separate' if r['separate'] else 'equal'}" for r in out["known"]["rows"])
    un = "; ".join(f"k_max {k} tol {t:g}: {n} of {tot} frames separate"
                   for k, t, n, tot in out["union_frames_separating_by_bound"])
    ai = (f"OBSERVED, known shapes (8-gon vs 16-gon, {len(out['known']['a']['points'])} and "
          f"{len(out['known']['b']['points'])} points): {kn}. OBSERVED, the Q12 unions "
          f"(8 frames, 10 or 16 points each): {un}. DERIVED: for the known pair the outputs "
          "are equal exactly where k_max is below the larger order, so there the equal "
          "reading is the probe's bound; the tolerance change moves nothing in this run. "
          "A set of n points cannot have rotational order above n, and every union has at "
          "most 16 points, so k_max 24 covers every order the unions could have. "
          "AI reading: not interpretable atm beyond the table.")
    p = emlib.write_run(NAME, seed,
                        {"k_maxes": list(K_MAXES), "rel_tols": list(REL_TOLS),
                         "known_shapes": [list(k) for k in KNOWN],
                         "unions_source": "explorative_math.projections.iter_projection_runs "
                                          "(dimension 4, frames 4, seed 0, composition phi), "
                                          "as in compare/cross.py",
                         "about": "centroid"},
                        out, None, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
