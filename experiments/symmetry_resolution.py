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

  unions_mirror   mirror_symmetry on the same two unions, at 360 and 720
           axes and both tolerances, with the bound of each reading.

The union parts import Build A (explorative_math) read-only, as compare/
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
N_AXES = (360, 720)
PREDICTION = ("DERIVED before the run, by the operator: on the two unions the mirror "
              "lists separate; the rotation lists are equal by construction.")
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


def mirror_grid(a, b):
    """mirror_symmetry on two point sets at every (n_axes, rel_tol)."""
    rows = []
    for n_axes in N_AXES:
        for tol in REL_TOLS:
            ra = Q.mirror_symmetry(a, n_axes, tol)
            rb = Q.mirror_symmetry(b, n_axes, tol)
            rows.append({"n_axes": n_axes, "rel_tol": tol, "a": ra, "b": rb,
                         "separate": ra != rb,
                         "a_bound": Q.mirror_symmetry_bound(a, n_axes, tol),
                         "b_bound": Q.mirror_symmetry_bound(b, n_axes, tol)})
    return rows


def union_frames():
    """(shape, frame index, Build A union, Build B union), as compare/cross.py makes them."""
    for name in ("simplex", "crosspolytope"):
        run = next(A_proj.iter_projection_runs(
            {"dimension": 4, "shape": name, "frames": 4, "seed": 0, "composition": "phi"}))
        for frame in [o for o in run["outputs"] if o["kind"] == "projection_frame"]:
            yield (name, frame["frame_index"], frame["composed_points"],
                   [r["point"] for r in C.phi_stack(frame["projected_points"], (0, 1))])


def unions_mirror():
    return [{"shape": name, "dimension": 4, "frame": i,
             "a_operation": "phi_scale_rotate (Build A)",
             "b_operation": "phi_scale (Build B)",
             "n_points": [len(a), len(b)],
             "rows": mirror_grid(a, b),
             "from_pairs": [{"rel_tol": tol,
                             "a": Q.mirror_symmetry_from_pairs(a, tol),
                             "b": Q.mirror_symmetry_from_pairs(b, tol),
                             "bound": {"axes": "every pair of points; no angle grid",
                                       "rel_tol": tol, "about": "centroid"}}
                            for tol in REL_TOLS]}
            for name, i, a, b in union_frames()]


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
    k, u, m = known(), unions(), unions_mirror()
    rot_equal = all(not r["separate"] for f in u for r in f["rows"])
    mir_frames = sum(1 for f in m if any(r["separate"] for r in f["rows"]))
    return {"known": k, "unions": u, "unions_mirror": m,
            "mirror_frames_separating_by_bound": [
                [n_axes, tol, sum(1 for f in m for r in f["rows"]
                                  if r["n_axes"] == n_axes and r["rel_tol"] == tol
                                  and r["separate"]), len(m)]
                for n_axes in N_AXES for tol in REL_TOLS],
            "mirror_frames_separating_from_pairs": [
                [tol, sum(1 for f in m for r in f["from_pairs"]
                          if r["rel_tol"] == tol and r["a"] != r["b"]), len(m)]
                for tol in REL_TOLS],
            "prediction": {"text": PREDICTION,
                           "frames_where_mirror_lists_separate_from_pairs":
                               sum(1 for f in m if any(r["a"] != r["b"]
                                                       for r in f["from_pairs"])),
                           "rotation_lists_equal_on_all_frames_and_bounds": rot_equal,
                           "frames_where_mirror_lists_separate_at_some_bound": mir_frames,
                           "frames": len(m)},
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
    mi = "; ".join(f"{n} axes tol {t:g}: {c} of {tot} frames separate"
                   for n, t, c, tot in out["mirror_frames_separating_by_bound"])
    per = "; ".join(f"{f['shape']} frame {f['frame']}: A {len(f['rows'][0]['a'])} axes, "
                    f"B {len(f['rows'][0]['b'])} axes" for f in out["unions_mirror"])
    fp = "; ".join(f"tol {t:g}: {c} of {tot} frames separate"
                   for t, c, tot in out["mirror_frames_separating_from_pairs"])
    fpper = "; ".join(f"{f['shape']} frame {f['frame']}: A {len(f['from_pairs'][0]['a'])}, "
                      f"B {len(f['from_pairs'][0]['b'])}" for f in out["unions_mirror"])
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
          f"OBSERVED, mirror_symmetry on the same unions: {mi}. Axis counts at 360 axes, "
          f"tolerance 1e-6: {per}. Prediction logged before the reading: {PREDICTION} "
          "DERIVED: the probe tests only axes on its grid (0.5 or 0.25 degrees apart) "
          "through the centroid, so a mirror line between grid angles reads as no mirror; "
          "an empty list there is bounded by that step. "
          f"OBSERVED, mirror axes with no grid (candidates from every pair of points): {fp}. "
          f"Axis counts, tolerance 1e-6: {fpper}. "
          "AI reading: not interpretable atm beyond the table.")
    p = emlib.write_run(NAME, seed,
                        {"k_maxes": list(K_MAXES), "rel_tols": list(REL_TOLS),
                         "known_shapes": [list(k) for k in KNOWN],
                         "unions_source": "explorative_math.projections.iter_projection_runs "
                                          "(dimension 4, frames 4, seed 0, composition phi), "
                                          "as in compare/cross.py",
                         "n_axes": list(N_AXES), "prediction": PREDICTION,
                         "about": "centroid"},
                        out, None, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
