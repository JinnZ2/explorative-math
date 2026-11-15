"""o. layered_pursuit -- descend toward more than one mirror, two ways.

CC0. Co-authored. Stdlib only.

runs/mirror_pursuit_seed0 descended the smallest axis residual and stopped
on one mirror (5-cell) or two (16-cell): that index is 0 as soon as one
mirror exists. Here the index asks for more axes. The score is
probes.mirror_k_score: the sum of the k smallest residuals over distinct
candidate axes (axes at least MIN_SEP_DEG apart), 0 when k mirrors exist.

  layered   stage s = 1, 2, ..., h descends the score with k = s, each
            stage starting where the one before stopped. The mirror-axis
            count is logged after every stage.
  alt       one descent of the score with k = h, from the same start.

"Layered" is read here as that staging. The operator's own definition was
not in the relay; if another staging was meant, this is a different run.

Same shapes, same random starts as mirror_pursuit (seed 0, 8 per shape),
same moves: one rotation angle at a time, +-step, first improvement kept,
step halved after a pass with no move. h is the Coxeter number (5, 8),
and the final counts are set beside the Coxeter-plane counts.

  python3 experiments/layered_pursuit.py [seed]
"""
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "layered_pursuit"
SHAPES = (("5-cell", lambda: P.simplex(4), 5), ("16-cell", lambda: P.cross_polytope(4), 8))
STARTS = 8
STEP0 = 30.0                # first step of every descent, degrees
STEP_MIN = 1e-6             # stop when the step is below this, degrees
SCORE_STOP = 1e-11          # stop when the score is below this
MAX_EVALS = 700             # per descent (per stage, and for alt)
MIN_SEP_DEG = 5.0           # two candidate axes closer than this are one axis
COUNT_TOLS = (1e-6, 1e-4)
PURSUIT = "runs/mirror_pursuit_seed0.json"
COXETER = "runs/coxeter_plane_seed0.json"


def project(verts, planes, angles_deg):
    return C.rotate_project(verts, [(pl, math.radians(a)) for pl, a in zip(planes, angles_deg)])


def read(pts):
    """Mirror-axis count at each tolerance, with what bounds it."""
    return [{"rel_tol": t, "n_mirror_axes": len(Q.mirror_symmetry_from_pairs(pts, t)),
             "mirror_axes_deg": Q.mirror_symmetry_from_pairs(pts, t)} for t in COUNT_TOLS]


def descend(verts, planes, start, k, max_evals=MAX_EVALS):
    ang = list(start)
    evals = 0

    def f(a):
        nonlocal evals
        evals += 1
        s = Q.mirror_k_score(project(verts, planes, a), k, MIN_SEP_DEG)
        return float("inf") if s is None else s

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
    return {"k": k, "start_score": path[0]["score"], "final_score": best,
            "final_angles_deg": ang, "passes": n_pass, "evals": evals, "stop": stop,
            "counts": read(pts),
            "n_distinct_projected_points": len({(round(x, 6), round(y, 6)) for x, y in pts}),
            "axis_residuals": Q.mirror_axis_residuals(pts, MIN_SEP_DEG)[:k + 2],
            "final_projected": pts, "path": path}


def layered(verts, planes, start, h, max_evals):
    stages, ang = [], list(start)
    for k in range(1, h + 1):
        st = descend(verts, planes, ang, k, max_evals)
        stages.append(st)
        ang = st["final_angles_deg"]
    return stages


def reference(path, pick):
    try:
        return pick(emlib.load_run(Path(path)))
    except FileNotFoundError:
        return None


def run(seed=0, starts=STARTS, max_evals=MAX_EVALS, pursuit_path=None, coxeter_path=None):
    rng = random.Random(seed)
    cox = reference(coxeter_path or emlib.ROOT / COXETER,
                    lambda r: {c["name"]: c["readings"][0]["n_mirror_axes"]
                               for c in r["outputs"]["cases"]}) or {}
    prev = reference(pursuit_path or emlib.ROOT / PURSUIT,
                     lambda r: {s["shape"]: [x["start_angles_deg"] for x in s["starts"]]
                                for s in r["outputs"]["shapes"]}) or {}
    out = []
    for name, make, h in SHAPES:
        verts = make()
        planes = P.all_planes(len(verts[0]))
        rows = []
        for i in range(starts):
            start = [rng.uniform(0.0, 360.0) for _ in planes]
            stages = layered(verts, planes, start, h, max_evals)
            alt = descend(verts, planes, start, h, max_evals)
            rows.append({"start_index": i, "start_angles_deg": start,
                         "start_counts": read(project(verts, planes, start)),
                         "layered_stages": stages,
                         "layered_axis_count_per_stage":
                             [s["counts"][0]["n_mirror_axes"] for s in stages],
                         "layered_final_count": stages[-1]["counts"][0]["n_mirror_axes"],
                         "alt": alt, "alt_final_count": alt["counts"][0]["n_mirror_axes"]})
        target = cox.get(name)
        same = (prev.get(name, [])[:starts] == [r["start_angles_deg"] for r in rows]
                if name in prev else None)
        out.append({"shape": name, "h": h, "n_vertices": len(verts),
                    "planes": [list(p) for p in planes],
                    "coxeter_plane_mirror_count": target,
                    "starts_equal_mirror_pursuit_starts": same,
                    "layered_final_counts": [r["layered_final_count"] for r in rows],
                    "layered_max_count_any_stage":
                        [max(r["layered_axis_count_per_stage"]) for r in rows],
                    "alt_final_counts": [r["alt_final_count"] for r in rows],
                    "layered_starts_reaching_coxeter_count":
                        (sum(1 for r in rows if r["layered_final_count"] == target)
                         if target is not None else None),
                    "alt_starts_reaching_coxeter_count":
                        (sum(1 for r in rows if r["alt_final_count"] == target)
                         if target is not None else None),
                    "starts": rows})
    return {"bound": {"starts_per_shape": starts, "step0_deg": STEP0, "step_min_deg": STEP_MIN,
                      "score_stop": SCORE_STOP, "max_evals_per_descent": max_evals,
                      "min_sep_deg": MIN_SEP_DEG, "count_rel_tols": list(COUNT_TOLS),
                      "moves": "one angle at a time, +-step, first improvement kept",
                      "axes": "mean polar angle of every pair of points",
                      "matching": "nearest point, not one-to-one"},
            "shapes": out}


def main(argv):
    seed = int(argv[0]) if argv else 0
    out = run(seed)
    symbols = [str(min(c, 9)) for s in out["shapes"] for r in s["starts"]
               for c in r["layered_axis_count_per_stage"]]
    comp = emlib.compressibility(
        symbols, "layered pursuit: mirror-axis count after each stage (digit, capped at 9), "
                 "start by start, 5-cell then 16-cell")
    rows = "; ".join(
        f"{s['shape']} (h = {s['h']}, Coxeter-plane count {s['coxeter_plane_mirror_count']}): "
        f"layered axis count per stage, one list per start "
        f"{[r['layered_axis_count_per_stage'] for r in s['starts']]}; layered final counts "
        f"{s['layered_final_counts']}, reaching the Coxeter count "
        f"{s['layered_starts_reaching_coxeter_count']} of {len(s['starts'])}; alt (k = h) "
        f"final counts {s['alt_final_counts']}, reaching it "
        f"{s['alt_starts_reaching_coxeter_count']} of {len(s['starts'])}; final scores "
        f"layered {[float(f'{r['layered_stages'][-1]['final_score']:.1e}') for r in s['starts']]}"
        f", alt {[float(f'{r['alt']['final_score']:.1e}') for r in s['starts']]}; distinct "
        f"projected points at the end, layered "
        f"{[r['layered_stages'][-1]['n_distinct_projected_points'] for r in s['starts']]}, alt "
        f"{[r['alt']['n_distinct_projected_points'] for r in s['starts']]} of {s['n_vertices']}; "
        f"same starts as mirror_pursuit: {s['starts_equal_mirror_pursuit_starts']}"
        for s in out["shapes"])
    stops = sorted({st["stop"] for s in out["shapes"] for r in s["starts"]
                    for st in r["layered_stages"] + [r["alt"]]})
    ai = (f"OBSERVED, counts at tolerance {COUNT_TOLS[0]:g}: {rows}. Stop reasons seen: "
          f"{stops}. DERIVED: a count is read at a tolerance, a score is not; a descent "
          "that stops on its evaluation limit with a small score can read fewer axes than "
          "the score suggests. Axes closer than the declared separation are one axis to "
          "the score. AI reading: not interpretable atm beyond the rows.")
    angles = {s["shape"]: {"layered": [r["layered_stages"][-1]["final_angles_deg"]
                                       for r in s["starts"]],
                           "alt": [r["alt"]["final_angles_deg"] for r in s["starts"]]}
              for s in out["shapes"]}
    p = emlib.write_run(NAME, seed,
                        {"shapes": [[n, h] for n, _, h in SHAPES], "starts": STARTS,
                         "step0_deg": STEP0, "step_min_deg": STEP_MIN,
                         "score_stop": SCORE_STOP, "max_evals_per_descent": MAX_EVALS,
                         "min_sep_deg": MIN_SEP_DEG, "count_rel_tols": list(COUNT_TOLS),
                         "pursuit_reference": PURSUIT, "coxeter_reference": COXETER,
                         "layered": "stage s descends the sum of the s smallest distinct-axis "
                                    "residuals, s = 1..h, from the previous stage's end",
                         "alt": "one descent of the sum of the h smallest, same start"},
                        out, angles, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
