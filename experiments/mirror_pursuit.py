"""n. mirror_pursuit -- descend the mirror residual from random starts.

CC0. Co-authored. Stdlib only.

A 4-polytope is rotated by one angle in each of the six coordinate planes
of R^4 and projected onto the first two coordinates. From a random start
the six angles are stepped, one at a time, to lower
probes.mirror_residual of the projected points (coordinate descent: try
+step and -step on each angle, keep a move that lowers the residual, halve
the step after a pass with no move). Every pass is logged.

At the end the projection is read for mirror axes
(mirror_symmetry_from_pairs) and set beside the Coxeter-plane count of
runs/coxeter_plane_seed0.json.

Prediction, PROPOSED before the run by the operator: asked as a question,
not a claim -- does the descent reach the Coxeter-plane mirror counts?

Related name, kept as a pointer: projection pursuit. The index pursued here
is the mirror residual.

  python3 experiments/mirror_pursuit.py [seed]
"""
import json
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "mirror_pursuit"
SHAPES = (("5-cell", lambda: P.simplex(4)), ("16-cell", lambda: P.cross_polytope(4)))
STARTS = 8                  # random starts per shape
STEP0 = 30.0                # first step, degrees
STEP_MIN = 1e-7             # stop when the step is below this, degrees
RESIDUAL_STOP = 1e-12       # stop when the residual is below this
MAX_EVALS = 4000            # stop after this many residual evaluations
COUNT_TOLS = (1e-6, 1e-4)   # tolerances the final mirror count is read at
COXETER = "runs/coxeter_plane_seed0.json"
QUESTION = ("PROPOSED before the run, by the operator, as a question: does the descent "
            "reach the Coxeter-plane mirror counts?")


def project(verts, planes, angles_deg):
    return C.rotate_project(verts, [(pl, math.radians(a)) for pl, a in zip(planes, angles_deg)])


def descend(verts, start_deg, step0=STEP0, step_min=STEP_MIN, max_evals=MAX_EVALS):
    planes = P.all_planes(len(verts[0]))
    ang = list(start_deg)
    evals = 0

    def f(a):
        nonlocal evals
        evals += 1
        return Q.mirror_residual(project(verts, planes, a))["residual"]

    best = f(ang)
    step = step0
    path = [{"pass": 0, "step_deg": step, "angles_deg": list(ang), "residual": best,
             "evals": evals}]
    n_pass = 0
    stop = "step below step_min"
    while step >= step_min:
        if best < RESIDUAL_STOP:
            stop = "residual below residual_stop"
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
        path.append({"pass": n_pass, "step_deg": step, "angles_deg": list(ang),
                     "residual": best, "evals": evals})
        if not moved:
            step /= 2.0
    pts = project(verts, planes, ang)
    counts = [{"rel_tol": t, "mirror_axes_deg": Q.mirror_symmetry_from_pairs(pts, t),
               "rotational_symmetry": Q.rotational_symmetry(pts, 24, t),
               "rotational_symmetry_bound": Q.rotational_symmetry_bound(pts, 24, t)}
              for t in COUNT_TOLS]
    for c in counts:
        c["n_mirror_axes"] = len(c["mirror_axes_deg"])
    return {"start_angles_deg": list(start_deg), "final_angles_deg": ang,
            "start_residual": path[0]["residual"], "final_residual": best,
            "final_reading": Q.mirror_residual(pts),
            "passes": n_pass, "evals": evals, "stop": stop,
            "final_projected": pts,
            "n_distinct_projected_points": len({(round(x, 6), round(y, 6)) for x, y in pts}),
            "final_counts": counts, "path": path}


def coxeter_counts(path=None):
    f = Path(path) if path else emlib.ROOT / COXETER
    try:
        rec = emlib.load_run(f)
    except FileNotFoundError:
        return {}
    return {c["name"]: c["readings"][0]["n_mirror_axes"] for c in rec["outputs"]["cases"]}


def run(seed=0, starts=STARTS, max_evals=MAX_EVALS, coxeter_path=None):
    rng = random.Random(seed)
    ref = coxeter_counts(coxeter_path)
    out = []
    for name, make in SHAPES:
        verts = make()
        planes = P.all_planes(len(verts[0]))
        runs = []
        for k in range(starts):
            start = [rng.uniform(0.0, 360.0) for _ in planes]
            runs.append(dict(descend(verts, start, max_evals=max_evals), start_index=k))
        target = ref.get(name)
        finals = [r["final_counts"][0]["n_mirror_axes"] for r in runs]
        out.append({"shape": name, "n_vertices": len(verts),
                    "planes": [list(p) for p in planes],
                    "coxeter_plane_mirror_count": target,
                    "final_mirror_counts": finals,
                    "final_residuals": [r["final_residual"] for r in runs],
                    "starts_reaching_coxeter_count":
                        (sum(1 for c in finals if c == target) if target is not None else None),
                    "starts_reaching_any_mirror": sum(1 for c in finals if c > 0),
                    "starts": runs})
    return {"question": QUESTION,
            "bound": {"starts_per_shape": starts, "step0_deg": STEP0, "step_min_deg": STEP_MIN,
                      "residual_stop": RESIDUAL_STOP, "max_evals": max_evals,
                      "count_rel_tols": list(COUNT_TOLS), "k_max": 24,
                      "moves": "one angle at a time, +-step, first improvement kept",
                      "shapes_not_run": ["tesseract", "24-cell"],
                      "shapes_not_run_reason": "residual cost per evaluation; not attempted"},
            "shapes": out}


def main(argv):
    seed = int(argv[0]) if argv else 0
    out = run(seed)
    symbols = [str(min(c, 9)) for s in out["shapes"] for c in s["final_mirror_counts"]]
    comp = emlib.compressibility(
        symbols, "final mirror-axis count per start (digit, capped at 9), 5-cell starts "
                 "then 16-cell starts")
    rows = "; ".join(
        f"{s['shape']}: final mirror counts at tolerance {COUNT_TOLS[0]:g} "
        f"{s['final_mirror_counts']} (at {COUNT_TOLS[1]:g}: "
        f"{[r['final_counts'][1]['n_mirror_axes'] for r in s['starts']]}), final residuals "
        f"{[float(f'{r:.2e}') for r in s['final_residuals']]}, distinct projected points "
        f"{[r['n_distinct_projected_points'] for r in s['starts']]} of {s['n_vertices']}, "
        f"Coxeter-plane count {s['coxeter_plane_mirror_count']}, starts reaching it "
        f"{s['starts_reaching_coxeter_count']} of {len(s['starts'])}"
        for s in out["shapes"])
    stops = sorted({r["stop"] for s in out["shapes"] for r in s["starts"]})
    ai = (f"OBSERVED, {STARTS} random starts per shape, seed {seed}: {rows}. Stop reasons: "
          f"{stops}. Question logged before the reading: {QUESTION} "
          "DERIVED: the residual is 0 at any projection with one mirror, so a descent "
          "that reaches 0 has found a mirror, not the largest number of mirrors; nothing "
          "in the index prefers more axes. Coordinate descent stops at the first point no "
          "single-angle move improves, which can be a local minimum above 0. "
          "AI reading: not interpretable atm beyond the rows.")
    angles = {s["shape"]: [r["final_angles_deg"] for r in s["starts"]] for s in out["shapes"]}
    p = emlib.write_run(NAME, seed,
                        {"shapes": [n for n, _ in SHAPES], "starts": STARTS,
                         "step0_deg": STEP0, "step_min_deg": STEP_MIN,
                         "residual_stop": RESIDUAL_STOP, "max_evals": MAX_EVALS,
                         "count_rel_tols": list(COUNT_TOLS),
                         "coxeter_reference": COXETER, "question": QUESTION},
                        out, angles, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
