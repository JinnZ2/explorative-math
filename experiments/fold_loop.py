"""x. fold_loop -- anchor, fold, rerun: search only inside the fundamental domain.

CC0. Co-authored. Stdlib only.

Label: fold_search_space. This run folds the SEARCH SPACE (the directions a
candidate axis can have). The carrier's method folds the SHAPE; that is
experiments/fold_shape.py. This run is kept as run, under this label, and
does not test the carrier's method.

Status: PROPOSED. The choices that make this reading concrete are listed
under "What fold means here" and logged in the run's bound.

  anchor   stage-1 pursuit (mirror_pursuit): one exact mirror, axis a, with
           its pairing kept. It is held from then on.
  fold     the mirrors found so far act on the directions an axis can
           have. A direction and its images under those mirrors carry the
           same residual, so only one of each family need be looked at:
           the fundamental domain, a wedge starting at a.
             n mirrors found, n odd    wedge width  90 / n degrees
             n mirrors found, n even   wedge width 180 / n degrees
  rerun    pursuit from the current position, every found mirror held
           exact (Gauss-Newton pull-back on the stacked constraints),
           lowering the smallest residual among candidate axes INSIDE the
           wedge. When a candidate comes close it is snapped: added to the
           held set and pulled back; if every held mirror is then exact,
           it is the next anchor. The mirrors now present are read, the
           wedge is folded again, and the loop repeats until no new axis
           is found or the count reaches h.

  compared on the same starts and seed
    a  full-space anchored generate: the logged constrained_generate and
       anchor_selection runs (rule m = h, 2500 trials)
    b  fold loop
    c  control, WRONG fold: a wedge of the same width placed at a random
       direction instead of at the found axis
    d  the same loop with no fold: the whole half-turn of directions

What fold means here
  - it folds the DIRECTIONS of candidate axes in the projection plane, not
    the six rotation angles. The rotation is restricted only by the held
    mirrors.
  - the loop uses no rule a + k * 180 / h: it does not know h, except to
    stop and to be scored.
  - a candidate within SEP_DEG of a found axis is not a new axis.

Every endpoint is read with probes.degeneracy. Observer slots null.

  python3 experiments/fold_loop.py [seed]
  python3 experiments/fold_loop.py --small [seed]   # phone bound; own run file
"""
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import constrained_generate as K  # noqa: E402
import emlib  # noqa: E402
import mirror_pursuit as U  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "fold_loop"
SHAPES = (("5-cell", lambda: P.simplex(4), 5), ("16-cell", lambda: P.cross_polytope(4), 8))
EXTRA_SHAPES = (("24-cell", P.cell24, 12),)
STARTS = 8
EXTRA_STARTS = 2            # starts for the extra shape
LEVEL_TRIALS = 400          # trial moves per level
MAX_LEVELS = 6
STEP0 = 30.0
STEP_MIN = 1e-5
SNAP_TRY = 0.02             # a candidate with residual below this is tried as the next anchor
ANCHOR_TOL = K.ANCHOR_TOL   # every held mirror stays below this
SEP_DEG = 3.0               # a candidate this close to a found axis is not a new axis
COUNT_TOL = 1e-6
VARIANTS = ("fold", "wrong_fold", "no_fold")
SMALL = {"starts": 2, "level_trials": 120, "anchor_max_evals": 1000, "extra_starts": 0}
REFERENCES = {"constrained_generate": "runs/constrained_generate_seed0.json",
              "anchor_selection": "runs/anchor_selection_seed0.json"}
LABEL = "fold_search_space"     # the carrier's method folds the shape: see fold_shape.py
METHOD = ("PROPOSED, carried, untested by the carrier: after an anchor is found, quotient the "
          "search by the symmetry found so far and search only the fundamental domain; the "
          "next axis becomes the next anchor; fold again; repeat.")


def wedge_width(n_axes):
    return 90.0 / n_axes if n_axes % 2 else 180.0 / n_axes


def circ(a, b):
    """Distance between two axis directions (period 180 degrees)."""
    d = abs(a - b) % 180.0
    return min(d, 180.0 - d)


def in_wedge(deg, lo, width):
    return (deg - lo) % 180.0 <= width + 1e-9


def wedge_score(pts, lo, width, found):
    """Smallest mirror residual among pair-candidate axes inside the wedge and
    not within SEP_DEG of a found axis. Returns (residual, axis, candidates
    in the wedge, per-axis residual sums actually carried out)."""
    n = len(pts)
    c = P.centroid(pts)
    q = [P.sub(p, c) for p in pts]
    scale = max((math.hypot(x, y) for x, y in q), default=0.0) or 1.0
    th = [math.atan2(y, x) for x, y in q if math.hypot(x, y) > 1e-12 * scale]
    uniq = {}
    for i, a in enumerate(th):
        for b in th[i:]:
            deg = math.degrees((a + b) / 2.0) % 180.0
            if in_wedge(deg, lo, width) and all(circ(deg, f) > SEP_DEG for f in found):
                uniq.setdefault(round(deg, 9), deg)
    best, best_deg = None, None
    for deg in uniq.values():
        a = math.radians(deg)
        c2, s2 = math.cos(2 * a), math.sin(2 * a)
        total = 0.0
        for x, y in q:
            rx, ry = c2 * x + s2 * y, s2 * x - c2 * y
            total += math.sqrt(min((rx - u) ** 2 + (ry - v) ** 2 for u, v in q))
            if best is not None and total >= best:
                break
        else:
            best, best_deg = total, deg
    return ((best / n / scale) if best is not None else None), best_deg, len(uniq)


def pull_back(verts, planes, angles, held):
    """Gauss-Newton onto every held mirror at once. held: [(axis_deg, pairing)]."""
    ang = list(angles)
    evals = 0

    def c_of(a):
        pts = K.project(verts, planes, a)
        out = []
        for axis, pair in held:
            out.extend(K.constraint(pts, axis, pair))
        return out
    for _ in range(K.GN_ITERS):
        c0 = c_of(ang)
        evals += 1
        if max(abs(v) for v in c0) < K.GN_STOP:
            break
        cols = []
        for i in range(len(ang)):
            t = list(ang)
            t[i] += K.FD_DEG
            ci = c_of(t)
            evals += 1
            cols.append([(b - a) / K.FD_DEG for a, b in zip(c0, ci)])
        n = len(ang)
        jtj = [[sum(u * v for u, v in zip(cols[i], cols[j])) for j in range(n)] for i in range(n)]
        lam = K.DAMPING * (sum(jtj[i][i] for i in range(n)) or 1.0)
        for i in range(n):
            jtj[i][i] += lam
        step = K.solve(jtj, [-sum(u * v for u, v in zip(cols[i], c0)) for i in range(n)])
        if step is None:
            break
        ang = [a + d for a, d in zip(ang, step)]
    return ang, evals


def held_ok(verts, planes, angles, held):
    pts = K.project(verts, planes, angles)
    return max(Q.mirror_residual_at(pts, [a for a, _ in held])) < ANCHOR_TOL


def level(verts, planes, angles, held, found, lo, width, level_trials):
    """One rerun: pursue inside the wedge until a candidate snaps, or stop."""
    ang = list(angles)
    proj = trials = cand_sum = scores = snaps_tried = 0

    def score(a):
        nonlocal cand_sum, scores
        r, deg, k = wedge_score(K.project(verts, planes, a), lo, width, found)
        cand_sum += k
        scores += 1
        return (float("inf") if r is None else r), deg

    def try_snap(a, deg):
        nonlocal proj, snaps_tried
        snaps_tried += 1
        pair = K.pairing(K.project(verts, planes, a), deg)
        trial_held = held + [(deg, pair)]
        t, e = pull_back(verts, planes, a, trial_held)
        proj += e
        return (t, trial_held) if held_ok(verts, planes, t, trial_held) else None

    best, bdeg = score(ang)
    start_score = best
    step, stop, snapped = STEP0, "step below step_min", None
    if best < SNAP_TRY:
        snapped = try_snap(ang, bdeg)
    while snapped is None and step >= STEP_MIN:
        if trials >= level_trials:
            stop = "level_trials reached"
            break
        moved = False
        for i in range(len(ang)):
            for d in (step, -step):
                t = list(ang)
                t[i] += d
                t, e = pull_back(verts, planes, t, held)
                proj += e
                trials += 1
                if not held_ok(verts, planes, t, held):
                    continue
                r, deg = score(t)
                if r < best:
                    ang, best, bdeg, moved = t, r, deg, True
                    if best < SNAP_TRY:
                        snapped = try_snap(ang, bdeg)
                    break
            if snapped is not None:
                break
        if not moved:
            step /= 2.0
    new_axis = None
    if snapped is not None:
        ang, held = snapped
        new_axis, stop = held[-1][0], "new axis snapped"
    return ang, held, {"wedge_start_deg": lo % 180.0, "domain_size_deg": width,
                       "axes_found_before": len(found), "start_score": start_score,
                       "final_score": best, "new_axis_deg": new_axis, "trials": trials,
                       "projection_evals": proj, "score_evals": scores,
                       "mean_candidate_axes_per_score": (cand_sum / scores) if scores else None,
                       "snaps_tried": snaps_tried, "stop": stop}


def loop(verts, planes, anc, h, variant, rng, level_trials, max_levels=MAX_LEVELS):
    a = anc["final_reading"]["axis_deg"]
    ang = list(anc["final_angles_deg"])
    held = [(a, K.pairing(K.project(verts, planes, ang), a))]
    offset = rng.uniform(SEP_DEG, 180.0 - SEP_DEG)     # drawn for every variant; used by one
    levels = []
    for _ in range(max_levels):
        pts = K.project(verts, planes, ang)
        found = Q.mirror_symmetry_from_pairs(pts, COUNT_TOL)
        if len(found) >= h:
            break
        n = max(len(found), 1)
        if variant == "no_fold":
            lo, width = a, 180.0
        elif variant == "wrong_fold":
            lo, width = a + offset, wedge_width(n)
        else:
            lo, width = a, wedge_width(n)
        ang, held, log = level(verts, planes, ang, held, found, lo, width, level_trials)
        pts = K.project(verts, planes, ang)
        after = Q.mirror_symmetry_from_pairs(pts, COUNT_TOL)
        log.update({"axes_found_after": len(after), "axes_after_deg": after,
                    "degeneracy": Q.degeneracy(pts),
                    "held_mirrors": len(held),
                    "max_held_residual": max(Q.mirror_residual_at(pts, [x for x, _ in held]))})
        levels.append(log)
        if log["new_axis_deg"] is None or len(after) <= log["axes_found_before"]:
            break
    pts = K.project(verts, planes, ang)
    axes = Q.mirror_symmetry_from_pairs(pts, COUNT_TOL)
    deg = Q.degeneracy(pts)
    return {"variant": variant, "wrong_fold_offset_deg": offset if variant == "wrong_fold" else None,
            "levels": levels, "n_levels": len(levels),
            "final_angles_deg": ang, "final_projected": pts,
            "final_mirror_axes_deg": axes, "final_count": len(axes),
            "final_count_at_1e-4": len(Q.mirror_symmetry_from_pairs(pts, 1e-4)),
            "reaches_h": len(axes) == h,
            "degeneracy": deg, "degenerate": deg["tags"]["1e-06"]["degenerate"],
            "held_pairings": [{"axis_deg": x, "pairing": p, "points_on_axis": K.points_on_axis(p)}
                              for x, p in held],
            "trials_total": sum(lv["trials"] for lv in levels),
            "projection_evals_total": sum(lv["projection_evals"] for lv in levels),
            "score_evals_total": sum(lv["score_evals"] for lv in levels),
            "candidate_axes_total": sum((lv["mean_candidate_axes_per_score"] or 0.0)
                                        * lv["score_evals"] for lv in levels),
            "stop": (levels[-1]["stop"] if levels else "anchor already has h axes")}


def references():
    out = {}
    for key, path in REFERENCES.items():
        try:
            rec = emlib.load_run(emlib.ROOT / path)
        except FileNotFoundError:
            out[key] = None
            continue
        per = {}
        for s in rec["outputs"]["shapes"]:
            c = s["compare"]
            if key == "constrained_generate":
                g = c["constrained_bound_check"]["h"]
                per[s["shape"]] = {"starts": [r["start_angles_deg"] for r in s["starts"]],
                                   "final_counts": g["final_counts"],
                                   "degenerate": g["degenerate"],
                                   "starts_reaching_h": g["starts_reaching_h"],
                                   "starts_reaching_h_excluding_degenerate":
                                       g["starts_reaching_h_excluding_degenerate"],
                                   "trials": g["trials"],
                                   "projection_evals": g["projection_evals"],
                                   "rule": "m = h, known in advance", "max_trials": 2500}
            else:
                g = c["selected"]["2500"]["h"]
                per[s["shape"]] = {"starts": [r["start_angles_deg"] for r in s["starts"]],
                                   "final_counts": g["final_counts"],
                                   "degenerate": g["degenerate"],
                                   "starts_reaching_h": g["starts_reaching_h"],
                                   "starts_reaching_h_excluding_degenerate":
                                       g["starts_reaching_h_excluding_degenerate"],
                                   "rule": "m = h, known in advance; anchors re-drawn by type",
                                   "max_trials": 2500}
        out[key] = per
    return out


def run_shape(name, make, h, starts, rng, level_trials, anchor_max_evals, ref, seed):
    verts = make()
    planes = P.all_planes(len(verts[0]))
    rows = []
    for i in range(starts):
        start = [rng.uniform(0.0, 360.0) for _ in planes]
        anc = U.descend(verts, start, max_evals=anchor_max_evals)
        variants = {}
        for v in VARIANTS:
            own = random.Random(f"{seed}/{name}/{i}/wrong-fold")     # same offset every time
            variants[v] = loop(verts, planes, anc, h, v, own, level_trials)
        rows.append({"start_index": i, "start_angles_deg": start,
                     "anchor": {"final_angles_deg": anc["final_angles_deg"],
                                "final_residual": anc["final_residual"],
                                "axis_deg": anc["final_reading"]["axis_deg"],
                                "evals": anc["evals"], "stop": anc["stop"],
                                "n_mirror_axes": len(Q.mirror_symmetry_from_pairs(
                                    anc["final_projected"], COUNT_TOL)),
                                "degeneracy": Q.degeneracy(anc["final_projected"])},
                     "variants": variants})

    def summ(v):
        gs = [r["variants"][v] for r in rows]
        return {"final_counts": [g["final_count"] for g in gs],
                "final_counts_at_1e-4": [g["final_count_at_1e-4"] for g in gs],
                "degenerate": [g["degenerate"] for g in gs],
                "n_degenerate_endpoints": sum(1 for g in gs if g["degenerate"]),
                "starts_reaching_h": sum(1 for g in gs if g["reaches_h"]),
                "starts_reaching_h_excluding_degenerate":
                    sum(1 for g in gs if g["reaches_h"] and not g["degenerate"]),
                "n_levels": [g["n_levels"] for g in gs],
                "axes_per_level": [[lv["axes_found_after"] for lv in g["levels"]] for g in gs],
                "domain_size_per_level": [[lv["domain_size_deg"] for lv in g["levels"]]
                                          for g in gs],
                "trials_total": [g["trials_total"] for g in gs],
                "projection_evals_total": [g["projection_evals_total"] for g in gs],
                "score_evals_total": [g["score_evals_total"] for g in gs],
                "candidate_axes_total": [g["candidate_axes_total"] for g in gs],
                "stops": [g["stop"] for g in gs]}
    same = {}
    for key, per in ref.items():
        r = (per or {}).get(name)
        same[key] = (r["starts"][:starts] == [x["start_angles_deg"] for x in rows]) if r else None
    return {"shape": name, "h": h, "n_vertices": len(verts),
            "compare": {"target_h": h, "count_rel_tol": COUNT_TOL,
                        "starts_equal_reference_starts": same,
                        "anchor_evals": [r["anchor"]["evals"] for r in rows],
                        "b_fold": summ("fold"), "c_wrong_fold": summ("wrong_fold"),
                        "d_no_fold": summ("no_fold"),
                        "a_full_space_anchored_generate":
                            {key: ({k: v for k, v in per[name].items() if k != "starts"}
                                   if per and name in per else None)
                             for key, per in ref.items()}},
            "starts": rows}


def run(seed=0, starts=STARTS, level_trials=LEVEL_TRIALS, anchor_max_evals=U.MAX_EVALS,
        extra_starts=EXTRA_STARTS):
    rng = random.Random(seed)
    ref = references()
    shapes = [run_shape(n, m, h, starts, rng, level_trials, anchor_max_evals, ref, seed)
              for n, m, h in SHAPES]
    extra_rng = random.Random(f"{seed}/extra-shapes")
    if extra_starts:
        shapes += [run_shape(n, m, h, extra_starts, extra_rng, level_trials, anchor_max_evals,
                             ref, seed) for n, m, h in EXTRA_SHAPES]
    return {"method": METHOD, "label": LABEL,
            "tests_the_carriers_method": False,
            "bound": {"starts_per_shape": starts, "extra_shapes": [n for n, _, _ in EXTRA_SHAPES]
                      if extra_starts else [], "extra_starts": extra_starts,
                      "extra_shapes_not_run": [] if extra_starts else
                      [n for n, _, _ in EXTRA_SHAPES],
                      "seed": seed, "level_trials": level_trials, "max_levels": MAX_LEVELS,
                      "step0_deg": STEP0, "step_min_deg": STEP_MIN, "snap_try": SNAP_TRY,
                      "anchor_tol": ANCHOR_TOL, "sep_deg": SEP_DEG, "count_rel_tol": COUNT_TOL,
                      "anchor_max_evals": anchor_max_evals,
                      "fold_acts_on": "directions of candidate axes in the projection plane; "
                                      "the rotation is restricted only by the held mirrors",
                      "wedge_width": "90 / n for n odd, 180 / n for n even; n = mirrors present",
                      "wrong_fold": "a wedge of that width at a random direction, one offset "
                                    "per start",
                      "no_fold": "all directions (180 degrees)",
                      "h_used_for": "stopping and scoring only; no rule a + k*180/h",
                      "pairing": "nearest-point, read when an axis is anchored, then kept",
                      "references_not_like_for_like": "a uses the rule m = h and 2500 trials; "
                          "b, c, d use no rule and level_trials per level",
                      "paths": "not kept; one record per level"},
            "shapes": shapes}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["starts"], SMALL["level_trials"], SMALL["anchor_max_evals"],
               SMALL["extra_starts"]) if small else run(seed))
    symbols = [str(min(c, 9)) for s in out["shapes"]
               for v in ("b_fold", "c_wrong_fold", "d_no_fold")
               for c in s["compare"][v]["final_counts"]]
    comp = emlib.compressibility(
        symbols, "final mirror-axis count per start (digit, capped at 9): per shape, fold, "
                 "wrong fold, no fold")
    angles = {s["shape"]: {v: [r["variants"][v]["final_mirror_axes_deg"] for r in s["starts"]]
                           for v in VARIANTS} for s in out["shapes"]}
    b = out["bound"]
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "label": LABEL,
                         "shapes": [[s["shape"], s["h"], len(s["starts"])] for s in out["shapes"]],
                         "level_trials": b["level_trials"], "max_levels": MAX_LEVELS,
                         "snap_try": SNAP_TRY, "sep_deg": SEP_DEG,
                         "variants": list(VARIANTS), "references": REFERENCES,
                         "method": METHOD, "status_before_run": "PROPOSED"},
                        out, angles, None, comp, None, gz=True)
    print(p)
    for s in out["shapes"]:
        c = s["compare"]
        print(s["shape"], "h =", s["h"], "starts", len(s["starts"]))
        for key, lab in (("b_fold", "b fold      "), ("c_wrong_fold", "c wrong fold"),
                         ("d_no_fold", "d no fold   ")):
            g = c[key]
            print(f"  {lab} counts {g['final_counts']} reach h {g['starts_reaching_h']} "
                  f"(non-degenerate {g['starts_reaching_h_excluding_degenerate']}; degenerate "
                  f"endpoints {g['n_degenerate_endpoints']})  trials {g['trials_total']}  "
                  f"candidate axes {[round(x) for x in g['candidate_axes_total']]}  axes per "
                  f"level {g['axes_per_level']}")
        for key, r in c["a_full_space_anchored_generate"].items():
            if r:
                print(f"  a {key}: counts {r['final_counts']} reach h "
                      f"{r['starts_reaching_h']}")


if __name__ == "__main__":
    main(sys.argv[1:])
