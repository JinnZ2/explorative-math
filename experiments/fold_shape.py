"""z. fold_shape -- anchor, fold the SHAPE, rerun on the folded points.

CC0. Co-authored. Stdlib only.

Carrier's method, untested by the carrier; status PROPOSED: fold the shape,
not the search space. (experiments/fold_loop.py folds the search space; it
is labelled fold_search_space and does not test this method.)

  fold      the mirrors found so far reflect every projected point into
            their fundamental wedge (angle 180 / n for n mirrors). Points
            that land together are one folded point; its multiplicity is
            how many landed there.
  rerun     pursuit on the FOLDED points only: the six rotation angles are
            stepped, every found mirror held exact, to lower the mirror
            residual of the folded set (mirror lines through the centre).
            A near mirror of the folded set is the candidate next anchor.
            It is then tried on the whole shape: added to the held mirrors
            and pulled back. If every held mirror is exact, it is the next
            anchor. If not, the folded set had a mirror the shape does not
            have; that is logged.
  iterate   fold again by all mirrors now present; repeat. Per level:
            folded point count, multiplicities, wedge angle, axes found.
            Stop at no new structure, or when the folded set is a single
            point.

  controls, same starts
    c  wrong fold: reflect across ONE random axis instead of the found ones
    d  no fold: the same loop on the unfolded points
  compared with the anchored generate (logged constrained_generate and
  anchor_selection runs, rule m = h).

  exact     where the shape allows it the fold is carried with no rounding
            (exact/fold.py): the Coxeter-plane projection of the 5-cell is
            a regular pentagon, whose points are exact in a + b*phi forms.
            The 16-cell and the 24-cell need sqrt(2) and sqrt(3), which
            exact/ does not hold; their folds are in floats, with the
            merging tolerance logged.

Every endpoint is read with probes.degeneracy. Observer slots null.

  python3 experiments/fold_shape.py [seed]
  python3 experiments/fold_shape.py --small [seed]   # phone bound; own run file
"""
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import constrained_generate as K  # noqa: E402
import emlib  # noqa: E402
import fold_loop as FL  # noqa: E402
import mirror_pursuit as U  # noqa: E402
from exact import Cyc5, ZETA  # noqa: E402
from exact import fold as XF  # noqa: E402
from exact import inflation as XI  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "fold_shape"
SHAPES = FL.SHAPES
EXTRA_SHAPES = FL.EXTRA_SHAPES
STARTS = 8
EXTRA_STARTS = 2
LEVEL_TRIALS = 400
MAX_LEVELS = 6
STEP0 = 30.0
STEP_MIN = 1e-5
SNAP_TRY = 0.02             # a folded-set mirror with residual below this is tried on the shape
ANCHOR_TOL = K.ANCHOR_TOL
SEP_DEG = 3.0               # a candidate this close to a found axis is not a new axis
MERGE_TOL = 1e-6            # folded points closer than this (x radius) are one folded point
COUNT_TOL = 1e-6
VARIANTS = ("fold_shape", "wrong_fold", "no_fold")
SMALL = {"starts": 2, "level_trials": 120, "anchor_max_evals": 1000, "extra_starts": 0}
METHOD = ("PROPOSED, carrier's method, untested by the carrier: reflect the shape into the "
          "fundamental wedge of the mirrors found; run the pursuit on the folded points; new "
          "structure there is the next anchor; fold again.")


# ---- the fold, in floats -------------------------------------------------------

def merge(points, scale):
    """Folded points closer than MERGE_TOL * scale are one. Returns
    (points, multiplicities)."""
    reps, mult = [], []
    for p in points:
        for k, r in enumerate(reps):
            if math.dist(p, r) < MERGE_TOL * scale:
                mult[k] += 1
                break
        else:
            reps.append(list(p))
            mult.append(1)
    return reps, mult


def fold_by_mirrors(pts, axes_deg):
    """Reflect every point into the wedge [a0, a0 + 180/n] of n equally spaced
    mirrors through the origin (a0 the first of them)."""
    n = len(axes_deg)
    if n == 0:
        return [list(p) for p in pts]
    a0 = math.radians(sorted(axes_deg)[0])
    period, wedge = 2 * math.pi / n, math.pi / n
    out = []
    for x, y in pts:
        r = math.hypot(x, y)
        t = (math.atan2(y, x) - a0) % period
        if t > wedge:
            t = period - t
        out.append([r * math.cos(a0 + t), r * math.sin(a0 + t)])
    return out


def fold_by_one_axis(pts, axis_deg):
    """Reflect every point on the clockwise side of one axis across it."""
    a = math.radians(axis_deg)
    ux, uy = math.cos(a), math.sin(a)
    c2, s2 = math.cos(2 * a), math.sin(2 * a)
    out = []
    for x, y in pts:
        if ux * y - uy * x < 0:
            x, y = c2 * x + s2 * y, s2 * x - c2 * y
        out.append([x, y])
    return out


def spacing_defect(axes_deg):
    """How far the found axes are from equally spaced (degrees); 0 for one axis."""
    n = len(axes_deg)
    if n < 2:
        return 0.0
    s = sorted(a % 180.0 for a in axes_deg)
    gaps = [(b - a) for a, b in zip(s, s[1:])] + [s[0] + 180.0 - s[-1]]
    return max(abs(g - 180.0 / n) for g in gaps)


def folded(pts, found, variant, wrong_axis):
    scale = max((math.hypot(x, y) for x, y in pts), default=0.0) or 1.0
    if variant == "no_fold":
        f = [list(p) for p in pts]
    elif variant == "wrong_fold":
        f = fold_by_one_axis(pts, wrong_axis)
    else:
        f = fold_by_mirrors(pts, found)
    reps, mult = merge(f, scale)
    return reps, mult, scale


def folded_score(reps, scale, found):
    """Smallest mirror residual of the folded set, over mirror lines through
    the ORIGIN whose direction is the mean polar angle of a pair of folded
    points and is not within SEP_DEG of a found axis. Returns (residual,
    axis, candidates)."""
    n = len(reps)
    th = [math.atan2(y, x) for x, y in reps if math.hypot(x, y) > 1e-12 * scale]
    uniq = {}
    for i, a in enumerate(th):
        for b in th[i:]:
            deg = math.degrees((a + b) / 2.0) % 180.0
            if all(FL.circ(deg, f) > SEP_DEG for f in found):
                uniq.setdefault(round(deg, 9), deg)
    best, best_deg = None, None
    for deg in uniq.values():
        a = math.radians(deg)
        c2, s2 = math.cos(2 * a), math.sin(2 * a)
        total = 0.0
        for x, y in reps:
            rx, ry = c2 * x + s2 * y, s2 * x - c2 * y
            total += math.sqrt(min((rx - u) ** 2 + (ry - v) ** 2 for u, v in reps))
            if best is not None and total >= best:
                break
        else:
            best, best_deg = total, deg
    return ((best / n / scale) if best is not None else None), best_deg, len(uniq)


# ---- one level, and the loop ---------------------------------------------------

def level(verts, planes, angles, held, found, variant, wrong_axis, level_trials,
          score_fn=None):
    """score_fn(points, found) -> (residual or None, axis, candidates, folded points):
    another way to score a position; None uses the folded-set mirror residual."""
    ang = list(angles)
    c = {"proj": 0, "trials": 0, "cand": 0, "scores": 0, "points": 0,
         "snaps_tried": 0, "snaps_refused": 0}

    def score(a):
        if score_fn is not None:
            r, deg, k, npts = score_fn(K.project(verts, planes, a), found)
        else:
            reps, _, scale = folded(K.project(verts, planes, a), found, variant, wrong_axis)
            r, deg, k = folded_score(reps, scale, found)
            npts = len(reps)
        c["cand"] += k
        c["scores"] += 1
        c["points"] += npts
        return (float("inf") if r is None else r), deg

    def try_snap(a, deg):
        c["snaps_tried"] += 1
        pair = K.pairing(K.project(verts, planes, a), deg)
        trial_held = held + [(deg, pair)]
        t, e = FL.pull_back(verts, planes, a, trial_held)
        c["proj"] += e
        if FL.held_ok(verts, planes, t, trial_held):
            return t, trial_held
        c["snaps_refused"] += 1
        return None

    best, bdeg = score(ang)
    start_score = best
    step, stop, snapped = STEP0, "step below step_min", None
    refused_at = []
    if best < SNAP_TRY:
        snapped = try_snap(ang, bdeg)
        if snapped is None:
            refused_at.append({"score": best, "axis_deg": bdeg})
    while snapped is None and step >= STEP_MIN:
        if c["trials"] >= level_trials:
            stop = "level_trials reached"
            break
        moved = False
        for i in range(len(ang)):
            for d in (step, -step):
                t = list(ang)
                t[i] += d
                t, e = FL.pull_back(verts, planes, t, held)
                c["proj"] += e
                c["trials"] += 1
                if not FL.held_ok(verts, planes, t, held):
                    continue
                r, deg = score(t)
                if r < best:
                    ang, best, bdeg, moved = t, r, deg, True
                    if best < SNAP_TRY:
                        snapped = try_snap(ang, bdeg)
                        if snapped is None and len(refused_at) < 20:
                            refused_at.append({"score": best, "axis_deg": bdeg})
                    break
            if snapped is not None:
                break
        if not moved:
            step /= 2.0
    new_axis = None
    if snapped is not None:
        ang, held = snapped
        new_axis, stop = held[-1][0], "new axis: folded-set mirror confirmed on the shape"
    return ang, held, {
        "start_score": start_score, "final_score": best, "new_axis_deg": new_axis,
        "trials": c["trials"], "projection_evals": c["proj"], "score_evals": c["scores"],
        "mean_candidate_axes_per_score": c["cand"] / c["scores"] if c["scores"] else None,
        "mean_folded_points_per_score": c["points"] / c["scores"] if c["scores"] else None,
        "candidate_axes_total": c["cand"],
        "folded_mirrors_tried_on_the_shape": c["snaps_tried"],
        "folded_mirrors_the_shape_does_not_have": c["snaps_refused"],
        "refused_examples": refused_at[:5], "stop": stop}


def loop(verts, planes, anc, h, variant, wrong_axis, level_trials, max_levels=MAX_LEVELS,
         score_fn=None):
    a = anc["final_reading"]["axis_deg"]
    ang = list(anc["final_angles_deg"])
    held = [(a, K.pairing(K.project(verts, planes, ang), a))]
    levels, stop = [], "max_levels reached"
    for _ in range(max_levels):
        pts = K.project(verts, planes, ang)
        found = Q.mirror_symmetry_from_pairs(pts, COUNT_TOL)
        reps, mult, _ = folded(pts, found, variant, wrong_axis)
        head = {"axes_found_before": len(found),
                "wedge_angle_deg": (180.0 / len(found) if variant == "fold_shape" and found
                                    else 180.0 if variant == "wrong_fold" else 360.0),
                "axes_spacing_defect_deg": spacing_defect(found),
                "folded_point_count": len(reps), "multiplicities": sorted(mult, reverse=True),
                "n_points": len(pts)}
        if len(reps) <= 1:
            levels.append(dict(head, stop="single point", trials=0, projection_evals=0,
                               score_evals=0, candidate_axes_total=0, new_axis_deg=None,
                               axes_found_after=len(found),
                               degeneracy=Q.degeneracy(pts)))
            stop = "single point"
            break
        ang, held, log = level(verts, planes, ang, held, found, variant, wrong_axis,
                               level_trials, score_fn)
        pts = K.project(verts, planes, ang)
        after = Q.mirror_symmetry_from_pairs(pts, COUNT_TOL)
        reps2, mult2, _ = folded(pts, after, variant, wrong_axis)
        log.update(head)
        log.update({"axes_found_after": len(after), "axes_after_deg": after,
                    "folded_point_count_after": len(reps2),
                    "multiplicities_after": sorted(mult2, reverse=True),
                    "degeneracy": Q.degeneracy(pts), "held_mirrors": len(held),
                    "max_held_residual": max(Q.mirror_residual_at(pts, [x for x, _ in held]))})
        levels.append(log)
        if log["new_axis_deg"] is None or len(after) <= log["axes_found_before"]:
            stop = "no new structure"
            break
    pts = K.project(verts, planes, ang)
    axes = Q.mirror_symmetry_from_pairs(pts, COUNT_TOL)
    reps, mult, _ = folded(pts, axes, variant, wrong_axis)
    deg = Q.degeneracy(pts)
    return {"variant": variant,
            "wrong_fold_axis_deg": wrong_axis if variant == "wrong_fold" else None,
            "levels": levels, "n_levels": len(levels), "stop": stop,
            "final_angles_deg": ang, "final_projected": pts,
            "final_mirror_axes_deg": axes, "final_count": len(axes),
            "final_count_at_1e-4": len(Q.mirror_symmetry_from_pairs(pts, 1e-4)),
            "final_folded_point_count": len(reps),
            "final_multiplicities": sorted(mult, reverse=True),
            "reaches_h": len(axes) == h,
            "degeneracy": deg, "degenerate": deg["tags"]["1e-06"]["degenerate"],
            "trials_total": sum(lv["trials"] for lv in levels),
            "projection_evals_total": sum(lv["projection_evals"] for lv in levels),
            "score_evals_total": sum(lv["score_evals"] for lv in levels),
            "candidate_axes_total": sum(lv["candidate_axes_total"] for lv in levels),
            "folded_mirrors_the_shape_does_not_have":
                sum(lv.get("folded_mirrors_the_shape_does_not_have", 0) for lv in levels)}


def run_shape(name, make, h, starts, rng, level_trials, anchor_max_evals, ref, seed):
    verts = make()
    planes = P.all_planes(len(verts[0]))
    rows = []
    for i in range(starts):
        start = [rng.uniform(0.0, 360.0) for _ in planes]
        anc = U.descend(verts, start, max_evals=anchor_max_evals)
        wrong = random.Random(f"{seed}/{name}/{i}/wrong-fold-axis").uniform(0.0, 180.0)
        rows.append({"start_index": i, "start_angles_deg": start,
                     "anchor": {"final_angles_deg": anc["final_angles_deg"],
                                "final_residual": anc["final_residual"],
                                "axis_deg": anc["final_reading"]["axis_deg"],
                                "evals": anc["evals"], "stop": anc["stop"],
                                "n_mirror_axes": len(Q.mirror_symmetry_from_pairs(
                                    anc["final_projected"], COUNT_TOL)),
                                "degeneracy": Q.degeneracy(anc["final_projected"])},
                     "variants": {v: loop(verts, planes, anc, h, v, wrong, level_trials)
                                  for v in VARIANTS}})

    def summ(v):
        gs = [r["variants"][v] for r in rows]
        return {"final_counts": [g["final_count"] for g in gs],
                "final_counts_at_1e-4": [g["final_count_at_1e-4"] for g in gs],
                "degenerate": [g["degenerate"] for g in gs],
                "n_degenerate_endpoints": sum(1 for g in gs if g["degenerate"]),
                "starts_reaching_h": sum(1 for g in gs if g["reaches_h"]),
                "starts_reaching_h_excluding_degenerate":
                    sum(1 for g in gs if g["reaches_h"] and not g["degenerate"]),
                "stops": [g["stop"] for g in gs],
                "n_levels": [g["n_levels"] for g in gs],
                "axes_per_level": [[lv["axes_found_after"] for lv in g["levels"]] for g in gs],
                "wedge_angle_per_level": [[lv["wedge_angle_deg"] for lv in g["levels"]]
                                          for g in gs],
                "folded_point_count_per_level": [[lv["folded_point_count"] for lv in g["levels"]]
                                                 for g in gs],
                "multiplicities_per_level": [[lv["multiplicities"] for lv in g["levels"]]
                                             for g in gs],
                "final_folded_point_counts": [g["final_folded_point_count"] for g in gs],
                "final_multiplicities": [g["final_multiplicities"] for g in gs],
                "folded_mirrors_the_shape_does_not_have":
                    [g["folded_mirrors_the_shape_does_not_have"] for g in gs],
                "trials_total": [g["trials_total"] for g in gs],
                "projection_evals_total": [g["projection_evals_total"] for g in gs],
                "score_evals_total": [g["score_evals_total"] for g in gs],
                "candidate_axes_total": [g["candidate_axes_total"] for g in gs]}
    same = {}
    for key, per in ref.items():
        r = (per or {}).get(name)
        same[key] = (r["starts"][:starts] == [x["start_angles_deg"] for x in rows]) if r else None
    return {"shape": name, "h": h, "n_vertices": len(verts),
            "compare": {"target_h": h, "count_rel_tol": COUNT_TOL,
                        "starts_equal_reference_starts": same,
                        "anchor_evals": [r["anchor"]["evals"] for r in rows],
                        "b_fold_shape": summ("fold_shape"), "c_wrong_fold": summ("wrong_fold"),
                        "d_no_fold": summ("no_fold"),
                        "a_anchored_generate":
                            {key: ({k: v for k, v in per[name].items() if k != "starts"}
                                   if per and name in per else None)
                             for key, per in ref.items()}},
            "starts": rows}


# ---- the fold with no rounding, where the shape allows it -----------------------

def exact_pentagon():
    """The 5-cell's Coxeter-plane projection is a regular pentagon. Its fold is
    carried in exact a + b*phi forms: one mirror, then all five."""
    pts = XI.pentagon()
    one = Cyc5(1, 0)
    z10 = -(ZETA * ZETA * ZETA)                       # direction at 36 degrees
    dirs = {"0": one, "36": z10, "72": ZETA, "108": z10 * ZETA, "144": ZETA * ZETA}
    mirrors = [d for d, u in dirs.items() if XF.is_mirror(pts, u)]
    f1 = XF.half_plane_fold(pts, one)
    f1_pts = [f for f, _ in f1]
    f1_mirrors = [d for d, u in dirs.items() if d != "0" and XF.is_mirror(f1_pts, u)]
    f5 = XF.fold(pts, one, z10)
    return {"shape": "regular pentagon (5-cell on its Coxeter plane)",
            "arithmetic": "exact a + b*phi forms; equality is ==, no tolerance",
            "n_points": len(pts),
            "mirror_directions_deg": mirrors,
            "level_1": {"fold": "one mirror (direction 0)", "wedge_angle_deg": 180,
                        "folded_point_count": len(f1),
                        "multiplicities": sorted((m for _, m in f1), reverse=True),
                        "folded_points": [f.pair() for f, _ in f1],
                        "mirrors_of_the_folded_set_deg": f1_mirrors,
                        "each_is_a_mirror_of_the_shape":
                            [XF.is_mirror(pts, dirs[d]) for d in f1_mirrors]},
            "level_2": {"fold": "all five mirrors (wedge between 0 and 36)",
                        "wedge_angle_deg": 36,
                        "folded_point_count": len(f5),
                        "multiplicities": sorted((m for _, m in f5), reverse=True),
                        "folded_points": [f.pair() for f, _ in f5],
                        "stop": "single point" if len(f5) == 1 else None},
            "not_exact_here": {"16-cell": "its Coxeter projection needs sqrt(2)",
                               "24-cell": "its Coxeter projection needs sqrt(3)"}}


def run(seed=0, starts=STARTS, level_trials=LEVEL_TRIALS, anchor_max_evals=U.MAX_EVALS,
        extra_starts=EXTRA_STARTS):
    rng = random.Random(seed)
    ref = FL.references()
    shapes = [run_shape(n, m, h, starts, rng, level_trials, anchor_max_evals, ref, seed)
              for n, m, h in SHAPES]
    extra_rng = random.Random(f"{seed}/extra-shapes")
    if extra_starts:
        shapes += [run_shape(n, m, h, extra_starts, extra_rng, level_trials, anchor_max_evals,
                             ref, seed) for n, m, h in EXTRA_SHAPES]
    return {"method": METHOD, "label": "fold_shape",
            "bound": {"starts_per_shape": starts, "extra_starts": extra_starts,
                      "extra_shapes": [n for n, _, _ in EXTRA_SHAPES] if extra_starts else [],
                      "extra_shapes_not_run": [] if extra_starts else
                      [n for n, _, _ in EXTRA_SHAPES],
                      "seed": seed, "level_trials": level_trials, "max_levels": MAX_LEVELS,
                      "step0_deg": STEP0, "step_min_deg": STEP_MIN, "snap_try": SNAP_TRY,
                      "anchor_tol": ANCHOR_TOL, "sep_deg": SEP_DEG, "merge_tol": MERGE_TOL,
                      "count_rel_tol": COUNT_TOL, "anchor_max_evals": anchor_max_evals,
                      "fold": "every projected point reflected into the wedge of the mirrors "
                              "found (assumed equally spaced; the defect is logged per level)",
                      "folded_mirror": "a line through the centre of the shape, direction from "
                                       "pairs of folded points, nearest-point matching, "
                                       "multiplicities not weighed",
                      "confirmation": "a folded-set mirror counts only if the whole shape has "
                                      "it: added to the held mirrors, pulled back, all exact",
                      "wrong_fold": "one random axis per start, the same at every level",
                      "no_fold": "the unfolded points",
                      "h_used_for": "scoring only",
                      "float_fold_tolerance": MERGE_TOL,
                      "references_not_like_for_like": "the anchored generate uses the rule "
                          "m = h and 2500 trials; this loop uses no rule and level_trials per "
                          "level",
                      "not_built": "B, a quotient of orientation space"},
            "exact_fold": exact_pentagon(),
            "shapes": shapes}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["starts"], SMALL["level_trials"], SMALL["anchor_max_evals"],
               SMALL["extra_starts"]) if small else run(seed))
    symbols = [str(min(c, 9)) for s in out["shapes"]
               for v in ("b_fold_shape", "c_wrong_fold", "d_no_fold")
               for c in s["compare"][v]["final_counts"]]
    comp = emlib.compressibility(
        symbols, "final mirror-axis count per start (digit, capped at 9): per shape, fold "
                 "shape, wrong fold, no fold")
    angles = {s["shape"]: {v: [r["variants"][v]["final_mirror_axes_deg"] for r in s["starts"]]
                           for v in VARIANTS} for s in out["shapes"]}
    b = out["bound"]
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "label": "fold_shape",
                         "shapes": [[s["shape"], s["h"], len(s["starts"])] for s in out["shapes"]],
                         "level_trials": b["level_trials"], "max_levels": MAX_LEVELS,
                         "snap_try": SNAP_TRY, "sep_deg": SEP_DEG, "merge_tol": MERGE_TOL,
                         "variants": list(VARIANTS), "references": FL.REFERENCES,
                         "method": METHOD, "status_before_run": "PROPOSED"},
                        out, angles, None, comp, None, gz=True)
    print(p)
    e = out["exact_fold"]
    print("exact pentagon: level 1", e["level_1"]["folded_point_count"], "points, multiplicities",
          e["level_1"]["multiplicities"], "folded-set mirrors", e["level_1"]["mirrors_of_the_folded_set_deg"],
          "each a mirror of the shape", e["level_1"]["each_is_a_mirror_of_the_shape"],
          "| level 2", e["level_2"]["folded_point_count"], "point, multiplicities",
          e["level_2"]["multiplicities"])
    for s in out["shapes"]:
        c = s["compare"]
        print(s["shape"], "h =", s["h"], "starts", len(s["starts"]))
        for key, lab in (("b_fold_shape", "b fold shape"), ("c_wrong_fold", "c wrong fold"),
                         ("d_no_fold", "d no fold   ")):
            g = c[key]
            print(f"  {lab} counts {g['final_counts']} reach h {g['starts_reaching_h']} "
                  f"(non-degenerate {g['starts_reaching_h_excluding_degenerate']}; degenerate "
                  f"endpoints {g['n_degenerate_endpoints']})  trials {g['trials_total']}  "
                  f"candidate axes {g['candidate_axes_total']}  refused on the shape "
                  f"{g['folded_mirrors_the_shape_does_not_have']}")
            print(f"       folded points per level {g['folded_point_count_per_level']}  final "
                  f"folded {g['final_folded_point_counts']}")
        for key, r in c["a_anchored_generate"].items():
            if r:
                print(f"  a {key}: counts {r['final_counts']} reach h {r['starts_reaching_h']}")


if __name__ == "__main__":
    main(sys.argv[1:])
