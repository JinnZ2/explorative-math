"""y. axis_split -- the 16-cell starts whose anchor has 2 points on the axis:
what is logged about the ones that reach 8 and the ones that do not.

CC0. Co-authored. Stdlib only.

In runs/anchor_selection_seed0 every 16-cell anchor has 2 points on the
axis; 5 starts reach 8 mirrors and 3 do not. This run remakes those eight
selected anchors and their descents (rule m = h, pairing kept, 2500
trials), this time keeping the path, and logs side by side for every
start:

  which points sit on the axis   vertex labels (+e0 .. -e3), and whether the
                                 two are an antipodal pair
  the pairing                    every swapped pair, with whether it swaps
                                 a vertex with its antipode
  the four generators            the 16-cell's points are +-p0 .. +-p3 in
                                 the plane; each p's radius and its angle
                                 from the anchor axis, at the anchor and
                                 at the end
  residual path                  score and anchor residual after every pass
  degeneracy                     at the anchor and at the end

The Coxeter-plane projection is logged the same way as a reference.
Nothing is interpreted: no reading is given and the observer slots are null.

  python3 experiments/axis_split.py [seed]
  python3 experiments/axis_split.py --small [seed]   # phone bound; own run file
"""
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anchor_selection as S  # noqa: E402
import constrained_generate as K  # noqa: E402
import emlib  # noqa: E402
import mirror_pursuit as U  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "axis_split"
SHAPE, H = "16-cell", 8
TRIALS = K.BOUND_CHECK_TRIALS
SMALL = {"starts": 2, "trials": 120, "anchor_max_evals": 1000, "max_attempts": 6}
SELECTION = "runs/anchor_selection_seed{seed}.json"


def label(i):
    return ("+" if i % 2 == 0 else "-") + f"e{i // 2}"


def antipode(i):
    return i ^ 1


def describe_pairing(pair):
    on = [i for i, j in enumerate(pair) if i == j]
    swaps = sorted({tuple(sorted((i, j))) for i, j in enumerate(pair) if i != j})
    return {"points_on_axis": [label(i) for i in on],
            "on_axis_indices": on,
            "on_axis_is_antipodal_pair": len(on) == 2 and antipode(on[0]) == on[1],
            "on_axis_generators": sorted({i // 2 for i in on}),
            "swapped_pairs": [[label(i), label(j)] for i, j in swaps],
            "swaps_with_own_antipode": sum(1 for i, j in swaps if antipode(i) == j),
            "swaps_between_generators": sum(1 for i, j in swaps if antipode(i) != j),
            "generator_permutation": sorted({tuple(sorted((i // 2, j // 2))) for i, j in swaps
                                             if i // 2 != j // 2})}


def generators(pts, axis_deg):
    """The four points p_i = image of +e_i: radius, and angle from the axis
    (degrees, in [0, 180))."""
    out = []
    for i in range(4):
        x, y = pts[2 * i]
        out.append({"generator": f"e{i}", "radius": math.hypot(x, y),
                    "angle_from_axis_deg": (math.degrees(math.atan2(y, x)) - axis_deg) % 180.0})
    radii = sorted(g["radius"] for g in out)
    return {"points": out, "radii_sorted": radii,
            "radius_ratio_min_over_max": radii[0] / radii[-1] if radii[-1] else None}


def coxeter_reference():
    verts = P.cross_polytope(4)
    basis = P.coxeter_plane([[1.0, -1.0, 0, 0], [0, 1.0, -1.0, 0], [0, 0, 1.0, -1.0],
                             [0, 0, 0, 1.0]], H)
    pts = P.project(verts, basis=basis)
    axes = Q.mirror_symmetry_from_pairs(pts, 1e-6)
    return {"n_mirror_axes": len(axes),
            "mirrors": [dict(describe_pairing(K.pairing(pts, a)), axis_deg=a) for a in axes],
            "generators": generators(pts, axes[0]),
            "degeneracy": Q.degeneracy(pts)}


def run(seed=0, starts=K.STARTS, trials=TRIALS, anchor_max_evals=U.MAX_EVALS,
        max_attempts=S.MAX_ATTEMPTS, selection_path=None):
    rng = random.Random(seed)
    # the 16-cell's starts come after the 5-cell's in every run's random stream
    for name, make, _ in K.SHAPES:
        if name == SHAPE:
            break
        for _ in range(starts):
            for _ in P.all_planes(len(make()[0])):
                rng.uniform(0.0, 360.0)
    verts = P.cross_polytope(4)
    planes = P.all_planes(4)
    want = S.PREFERRED[SHAPE]
    try:
        logged = emlib.load_run(Path(selection_path) if selection_path
                                else emlib.ROOT / SELECTION.format(seed=seed))
        lc = next(s["compare"] for s in logged["outputs"]["shapes"] if s["shape"] == SHAPE)
        key = str(trials)
        logged_counts = (lc["selected"][key]["h"]["final_counts"]
                         if key in lc["selected"] else None)
    except FileNotFoundError:
        logged_counts = None
    rows = []
    for i in range(starts):
        start = [rng.uniform(0.0, 360.0) for _ in planes]
        extra = random.Random(f"{seed}/{SHAPE}/{i}/re-anchor")
        attempts = [S.anchor_of(verts, planes, start, anchor_max_evals)]
        while attempts[-1]["points_on_axis"] != want and len(attempts) < max_attempts:
            fresh = [extra.uniform(0.0, 360.0) for _ in planes]
            attempts.append(S.anchor_of(verts, planes, fresh, anchor_max_evals))
        anc = attempts[-1]
        a = anc["axis_deg"]
        g = K.constrained_descent(verts, planes, anc["final_angles_deg"], a, H, trials)
        pts0 = K.project(verts, planes, anc["final_angles_deg"])
        pts1 = g["final_projected"]
        count = g["counts"][0]["n_mirror_axes"]
        rows.append({
            "start_index": i, "re_anchored": len(attempts) > 1, "attempts": len(attempts),
            "reaches_h": count == H, "final_count": count,
            "final_count_at_1e-4": g["counts"][1]["n_mirror_axes"],
            "final_mirror_axes_deg": g["counts"][0]["mirror_axes_deg"],
            "anchor_axis_deg": a,
            "anchor_pairing": describe_pairing(anc["pairing"]),
            "final_pairing": describe_pairing(g["final_pairing"]),
            "anchor_mirror_count": anc["counts"][0]["n_mirror_axes"],
            "generators_at_anchor": generators(pts0, a),
            "generators_at_end": generators(pts1, a),
            "degeneracy_at_anchor": Q.degeneracy(pts0),
            "degeneracy_at_end": g["degeneracy"],
            "anchor_residual": anc["final_residual"],
            "anchor_stop": anc["stop"],
            "start_score": g["start_score"], "final_score": g["final_score"],
            "final_residual_per_axis": g["final_residual_per_axis"],
            "generated_axes_deg": g["generated_axes_deg"],
            "trials": g["trials"], "moves_accepted": g["moves_accepted"],
            "trials_refused_by_constraint": g["trials_refused_by_constraint"],
            "stop": g["stop"],
            "residual_path": [{"pass": p["pass"], "step_deg": p["step_deg"], "score": p["score"],
                               "anchor_residual": p["anchor_residual"]} for p in g["path"]],
            "angles_path": [p["angles_deg"] for p in g["path"]]})
    reached = [r for r in rows if r["reaches_h"]]
    failed = [r for r in rows if not r["reaches_h"]]

    def side(rs):
        return {"starts": [r["start_index"] for r in rs],
                "final_counts": [r["final_count"] for r in rs],
                "on_axis_points": [r["anchor_pairing"]["points_on_axis"] for r in rs],
                "on_axis_is_antipodal_pair":
                    [r["anchor_pairing"]["on_axis_is_antipodal_pair"] for r in rs],
                "swaps_with_own_antipode":
                    [r["anchor_pairing"]["swaps_with_own_antipode"] for r in rs],
                "swaps_between_generators":
                    [r["anchor_pairing"]["swaps_between_generators"] for r in rs],
                "anchor_mirror_count": [r["anchor_mirror_count"] for r in rs],
                "anchor_radius_ratio_min_over_max":
                    [r["generators_at_anchor"]["radius_ratio_min_over_max"] for r in rs],
                "end_radius_ratio_min_over_max":
                    [r["generators_at_end"]["radius_ratio_min_over_max"] for r in rs],
                "anchor_min_pair_distance":
                    [r["degeneracy_at_anchor"]["min_pair_distance"] for r in rs],
                "end_min_pair_distance":
                    [r["degeneracy_at_end"]["min_pair_distance"] for r in rs],
                "degenerate_at_end":
                    [r["degeneracy_at_end"]["tags"]["1e-06"]["degenerate"] for r in rs],
                "start_scores": [r["start_score"] for r in rs],
                "final_scores": [r["final_score"] for r in rs],
                "passes": [len(r["residual_path"]) - 1 for r in rs],
                "moves_accepted": [r["moves_accepted"] for r in rs],
                "stops": [r["stop"] for r in rs],
                "re_anchored": [r["re_anchored"] for r in rs]}
    return {"shape": SHAPE, "h": H,
            "bound": {"seed": seed, "starts": starts, "trials": trials,
                      "preferred_points_on_axis": want, "max_attempts": max_attempts,
                      "anchor_max_evals": anchor_max_evals, "count_rel_tol": K.COUNT_TOLS[0],
                      "pairing": "kept from the anchor", "rule": "m = h",
                      "vertex_labels": [label(i) for i in range(8)]},
            "counts_equal_logged_selection_run":
                (logged_counts[:starts] == [r["final_count"] for r in rows])
                if logged_counts is not None else None,
            "logged_selection_counts": logged_counts,
            "all_anchors_have_preferred_type":
                all(len(r["anchor_pairing"]["on_axis_indices"]) == want for r in rows),
            "reached": side(reached), "not_reached": side(failed),
            "coxeter_plane_reference": coxeter_reference(),
            "starts": rows}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["starts"], SMALL["trials"], SMALL["anchor_max_evals"],
               SMALL["max_attempts"]) if small else run(seed))
    symbols = ["r" if r["reaches_h"] else "n" for r in out["starts"]]
    comp = emlib.compressibility(symbols, "one letter per start: r reaches h, n does not")
    angles = {str(r["start_index"]): r["final_mirror_axes_deg"] for r in out["starts"]}
    # observer slots null, and no reading: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "shape": SHAPE, "h": H,
                         "starts": out["bound"]["starts"], "trials": out["bound"]["trials"],
                         "selection_reference": SELECTION.format(seed=seed),
                         "asked_for": "what separates the 2-on-axis starts that reach 8 from "
                                      "those that do not; logged, not interpreted"},
                        out, angles, None, comp, None, gz=True)
    print(p)
    print("counts equal the logged selection run:", out["counts_equal_logged_selection_run"])
    for key in ("reached", "not_reached"):
        s = out[key]
        print(key, "starts", s["starts"], "counts", s["final_counts"])
        print("   on axis", s["on_axis_points"], "antipodal", s["on_axis_is_antipodal_pair"])
        print("   swaps with own antipode", s["swaps_with_own_antipode"], "between generators",
              s["swaps_between_generators"], "anchor mirrors", s["anchor_mirror_count"])
        print("   radius ratio at anchor", [round(x, 4) for x in s["anchor_radius_ratio_min_over_max"]],
              "at end", [round(x, 4) for x in s["end_radius_ratio_min_over_max"]])
        print("   min pair distance at anchor", [float(f"{x:.2e}") for x in s["anchor_min_pair_distance"]],
              "at end", [float(f"{x:.2e}") for x in s["end_min_pair_distance"]],
              "degenerate at end", s["degenerate_at_end"])
        print("   start score", [round(x, 3) for x in s["start_scores"]], "final score",
              [float(f"{x:.1e}") for x in s["final_scores"]], "passes", s["passes"],
              "stops", s["stops"])
    c = out["coxeter_plane_reference"]
    print("coxeter reference: on-axis counts", [len(m["on_axis_indices"]) for m in c["mirrors"]],
          "antipodal", [m["on_axis_is_antipodal_pair"] for m in c["mirrors"]],
          "swaps with own antipode", [m["swaps_with_own_antipode"] for m in c["mirrors"]],
          "radius ratio", round(c["generators"]["radius_ratio_min_over_max"], 6))


if __name__ == "__main__":
    main(sys.argv[1:])
