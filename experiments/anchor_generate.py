"""p. anchor_generate -- find one mirror, then generate the other axes by rule.

CC0. Co-authored. Stdlib only.

Hypothesis, PROPOSED by the operator: find one anchor, then generate the
rest by rule, instead of searching for each mirror.

  anchor     the stage-1 pursuit of mirror_pursuit (descent on
             probes.mirror_residual), from the same starts. Its end point
             has one mirror; a is that mirror's axis angle.
  generate   axes a + k * 180/m degrees, k = 0..m-1, held fixed in the
             projection plane. Score = sum of probes.mirror_residual_at on
             those axes only. The score is descended from the anchor's
             angles, same moves as the other pursuits.
               m = h        the dihedral rule (Coxeter-plane pattern)
               m = h - 1    control, wrong rule
               m = h + 1    control, wrong rule
  compare    beside runs/layered_pursuit_seed0 (layered, and alt k = h):
             starts reaching h mirrors, evaluations used, final residual
             per axis.

Shapes: 5-cell (h = 5), 16-cell (h = 8); 8 starts each, seed 0, the starts
of mirror_pursuit and layered_pursuit. Nothing is interpreted here; both
observer slots are written null.

  python3 experiments/anchor_generate.py [seed]
"""
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
import mirror_pursuit as U  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "anchor_generate"
SHAPES = (("5-cell", lambda: P.simplex(4), 5), ("16-cell", lambda: P.cross_polytope(4), 8))
STARTS = 8
STEP0 = 30.0                # first step of the generate descent, degrees
STEP_MIN = 1e-6             # stop when the step is below this, degrees
SCORE_STOP = 1e-11          # stop when the score is below this
MAX_EVALS = 700             # per generate descent (the layered_pursuit limit)
BOUND_CHECK_EVALS = 4000    # the same descent again at a higher limit, logged beside it
COUNT_TOLS = (1e-6, 1e-4)
LAYERED = "runs/layered_pursuit_seed0.json"
HYPOTHESIS = ("PROPOSED by the operator: find one anchor, then generate the rest by rule, "
              "instead of searching for each mirror. Controls (m = h - 1, h + 1) are "
              "expected to fail.")


def project(verts, planes, angles_deg):
    return C.rotate_project(verts, [(pl, math.radians(a)) for pl, a in zip(planes, angles_deg)])


def counts(pts):
    return [{"rel_tol": t, "mirror_axes_deg": Q.mirror_symmetry_from_pairs(pts, t),
             "n_mirror_axes": len(Q.mirror_symmetry_from_pairs(pts, t))} for t in COUNT_TOLS]


def generated_axes(anchor_deg, m):
    return [(anchor_deg + k * 180.0 / m) % 180.0 for k in range(m)]


def generate(verts, planes, anchor_angles, axes, max_evals=MAX_EVALS):
    ang = list(anchor_angles)
    evals = 0

    def f(a):
        nonlocal evals
        evals += 1
        return sum(Q.mirror_residual_at(project(verts, planes, a), axes))

    best = f(ang)
    step, n_pass = STEP0, 0
    path = [{"pass": 0, "step_deg": step, "angles_deg": list(ang), "score": best}]
    stop = "step below step_min"
    while step >= STEP_MIN:
        if best < SCORE_STOP:
            stop = "score below score_stop"
            break
        if evals >= max_evals:
            stop = "max_evals reached"
            break
        moved = False
        for i in range(len(ang)):
            for d in (step, -step):
                trial = list(ang)
                trial[i] += d
                r = f(trial)
                if r < best:
                    ang, best, moved = trial, r, True
                    break
        n_pass += 1
        path.append({"pass": n_pass, "step_deg": step, "angles_deg": list(ang), "score": best})
        if not moved:
            step /= 2.0
    pts = project(verts, planes, ang)
    return {"m": len(axes), "generated_axes_deg": list(axes),
            "start_score": path[0]["score"], "final_score": best,
            "final_residual_per_axis": Q.mirror_residual_at(pts, axes),
            "final_angles_deg": ang, "passes": n_pass, "evals": evals, "stop": stop,
            "counts": counts(pts),
            "n_distinct_projected_points": len({(round(x, 6), round(y, 6)) for x, y in pts}),
            "final_projected": pts, "path": path}


def layered_reference(path=None):
    try:
        rec = emlib.load_run(Path(path) if path else emlib.ROOT / LAYERED)
    except FileNotFoundError:
        return {}
    out = {}
    for s in rec["outputs"]["shapes"]:
        out[s["shape"]] = {
            "starts": [r["start_angles_deg"] for r in s["starts"]],
            "layered_final_counts": s["layered_final_counts"],
            "layered_evals": [sum(st["evals"] for st in r["layered_stages"])
                              for r in s["starts"]],
            "layered_final_residual_per_axis": [
                [a["residual"] for a in r["layered_stages"][-1]["axis_residuals"][:s["h"]]]
                for r in s["starts"]],
            "alt_final_counts": s["alt_final_counts"],
            "alt_evals": [r["alt"]["evals"] for r in s["starts"]],
            "alt_final_residual_per_axis": [
                [a["residual"] for a in r["alt"]["axis_residuals"][:s["h"]]]
                for r in s["starts"]]}
    return out


def run(seed=0, starts=STARTS, max_evals=MAX_EVALS, anchor_max_evals=U.MAX_EVALS,
        layered_path=None, bound_check_evals=BOUND_CHECK_EVALS):
    rng = random.Random(seed)
    ref = layered_reference(layered_path)
    shapes = []
    for name, make, h in SHAPES:
        verts = make()
        planes = P.all_planes(len(verts[0]))
        rows = []
        for i in range(starts):
            start = [rng.uniform(0.0, 360.0) for _ in planes]
            anc = U.descend(verts, start, max_evals=anchor_max_evals)
            a = anc["final_reading"]["axis_deg"]
            anchor = {"final_angles_deg": anc["final_angles_deg"],
                      "final_residual": anc["final_residual"], "axis_deg": a,
                      "evals": anc["evals"], "stop": anc["stop"],
                      "counts": counts(anc["final_projected"])}
            rules = {}
            for label, m in (("h_minus_1", h - 1), ("h", h), ("h_plus_1", h + 1)):
                rules[label] = generate(verts, planes, anc["final_angles_deg"],
                                        generated_axes(a, m), max_evals)
                if bound_check_evals:
                    again = generate(verts, planes, anc["final_angles_deg"],
                                     generated_axes(a, m), bound_check_evals)
                    rules[label]["bound_check"] = {
                        "max_evals": bound_check_evals, "evals": again["evals"],
                        "stop": again["stop"], "final_score": again["final_score"],
                        "final_residual_per_axis": again["final_residual_per_axis"],
                        "final_angles_deg": again["final_angles_deg"],
                        "counts": again["counts"]}
            rows.append({"start_index": i, "start_angles_deg": start, "anchor": anchor,
                         "rules": rules})
        r = ref.get(name)
        same = (r["starts"][:starts] == [x["start_angles_deg"] for x in rows]) if r else None

        def col(label, key):
            return [x["rules"][label][key] for x in rows]

        def reach(label):
            return sum(1 for x in rows if x["rules"][label]["counts"][0]["n_mirror_axes"] == h)
        compare = {
            "target_h": h, "count_rel_tol": COUNT_TOLS[0],
            "starts_equal_layered_pursuit_starts": same,
            "anchor_generate": {
                label: {"m": rows[0]["rules"][label]["m"],
                        "final_counts": [c[0]["n_mirror_axes"] for c in col(label, "counts")],
                        "starts_reaching_h": reach(label),
                        "generate_evals": col(label, "evals"),
                        "anchor_evals": [x["anchor"]["evals"] for x in rows],
                        "total_evals": [x["anchor"]["evals"] + x["rules"][label]["evals"]
                                        for x in rows],
                        "final_scores": col(label, "final_score"),
                        "final_residual_per_axis": col(label, "final_residual_per_axis"),
                        "stops": col(label, "stop"),
                        "bound_check": ({
                            "max_evals": bound_check_evals,
                            "final_counts": [x["rules"][label]["bound_check"]["counts"][0]
                                             ["n_mirror_axes"] for x in rows],
                            "starts_reaching_h": sum(
                                1 for x in rows if x["rules"][label]["bound_check"]["counts"][0]
                                ["n_mirror_axes"] == h),
                            "generate_evals": [x["rules"][label]["bound_check"]["evals"]
                                               for x in rows],
                            "final_scores": [x["rules"][label]["bound_check"]["final_score"]
                                             for x in rows],
                            "stops": [x["rules"][label]["bound_check"]["stop"] for x in rows]}
                            if bound_check_evals else None)}
                for label in ("h_minus_1", "h", "h_plus_1")},
            "layered": ({"final_counts": r["layered_final_counts"],
                         "starts_reaching_h": sum(1 for c in r["layered_final_counts"] if c == h),
                         "total_evals": r["layered_evals"],
                         "final_residual_per_axis": r["layered_final_residual_per_axis"]}
                        if r else None),
            "alt": ({"final_counts": r["alt_final_counts"],
                     "starts_reaching_h": sum(1 for c in r["alt_final_counts"] if c == h),
                     "total_evals": r["alt_evals"],
                     "final_residual_per_axis": r["alt_final_residual_per_axis"]}
                    if r else None)}
        shapes.append({"shape": name, "h": h, "n_vertices": len(verts),
                       "planes": [list(p) for p in planes],
                       "compare": compare, "starts": rows})
    return {"hypothesis": HYPOTHESIS,
            "bound": {"starts_per_shape": starts, "seed": seed,
                      "anchor": {"index": "mirror_residual", "step0_deg": U.STEP0,
                                 "step_min_deg": U.STEP_MIN, "residual_stop": U.RESIDUAL_STOP,
                                 "max_evals": anchor_max_evals},
                      "generate": {"step0_deg": STEP0, "step_min_deg": STEP_MIN,
                                   "score_stop": SCORE_STOP, "max_evals": max_evals,
                                   "bound_check_max_evals": bound_check_evals,
                                   "axes": "fixed in the projection plane at a + k*180/m",
                                   "matching": "nearest point, not one-to-one"},
                      "count_rel_tols": list(COUNT_TOLS),
                      "moves": "one angle at a time, +-step, first improvement kept",
                      "evals_not_like_for_like": "an anchor evaluation reads every pair "
                          "candidate; a generate evaluation reads m axes; a layered or alt "
                          "evaluation reads every pair candidate"},
            "shapes": shapes}


def main(argv):
    seed = int(argv[0]) if argv else 0
    out = run(seed)
    symbols = [str(min(c, 9)) for s in out["shapes"]
               for label in ("h_minus_1", "h", "h_plus_1")
               for c in s["compare"]["anchor_generate"][label]["final_counts"]]
    comp = emlib.compressibility(
        symbols, "final mirror-axis count per start (digit, capped at 9): per shape, "
                 "rule h-1 then h then h+1")
    angles = {s["shape"]: {label: [r["rules"][label]["final_angles_deg"] for r in s["starts"]]
                           for label in ("h_minus_1", "h", "h_plus_1")}
              for s in out["shapes"]}
    # observer slots null by instruction: ai_interest is void
    p = emlib.write_run(NAME, seed,
                        {"shapes": [[n, h] for n, _, h in SHAPES], "starts": STARTS,
                         "rules": ["h - 1", "h", "h + 1"], "step0_deg": STEP0,
                         "step_min_deg": STEP_MIN, "score_stop": SCORE_STOP,
                         "max_evals_generate": MAX_EVALS, "max_evals_anchor": U.MAX_EVALS,
                         "bound_check_max_evals": BOUND_CHECK_EVALS,
                         "count_rel_tols": list(COUNT_TOLS),
                         "layered_reference": LAYERED, "hypothesis": HYPOTHESIS},
                        out, angles, None, comp, None, gz=True)
    print(p)
    for s in out["shapes"]:
        c = s["compare"]
        print(s["shape"], "h =", s["h"], "same starts:", c["starts_equal_layered_pursuit_starts"])
        for label in ("h_minus_1", "h", "h_plus_1"):
            g = c["anchor_generate"][label]
            print(f"  rule {label:9s} m={g['m']}  counts {g['final_counts']}  reach h "
                  f"{g['starts_reaching_h']}  evals {g['total_evals']}  scores "
                  f"{[float(f'{x:.1e}') for x in g['final_scores']]}")
            b = g["bound_check"]
            print(f"       at {b['max_evals']} evals: counts {b['final_counts']}  reach h "
                  f"{b['starts_reaching_h']}  stops on limit "
                  f"{sum(1 for x in b['stops'] if x == 'max_evals reached')}  scores "
                  f"{[float(f'{x:.1e}') for x in b['final_scores']]}")
        for k in ("layered", "alt"):
            if c[k]:
                print(f"  {k:14s} counts {c[k]['final_counts']}  reach h "
                      f"{c[k]['starts_reaching_h']}  evals {c[k]['total_evals']}")


if __name__ == "__main__":
    main(sys.argv[1:])
