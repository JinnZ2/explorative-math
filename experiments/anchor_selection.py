"""w. anchor_selection -- choose the anchor by its pairing type before generating.

CC0. Co-authored. Stdlib only.

An anchor is a position (the six rotation angles) AND a pairing (which
point is the mirror image of which). The pairing has a type: the number
of points on the mirror axis. In runs/constrained_generate_seed0 the
starts that reached the full count were those whose anchor had the type of
a Coxeter-plane mirror.

Hypothesis, PROPOSED by the operator: classify the anchor by points on the
axis before generating; prefer the type that matched the full symmetry;
re-anchor if the anchor is of another type.

  random anchor     the logged constrained_generate run: the anchor each
                    start happened to find (pairing kept).
  selected anchor   the same starts. An anchor of the preferred type is
                    used as found. Otherwise the stage-1 pursuit is run
                    again from fresh random angles until it lands on the
                    preferred type, or MAX_ATTEMPTS is used up (then the
                    last anchor is used and the start is marked).
  preferred type    16-cell: 2 points on the axis. 5-cell: 1 point on the
                    axis. Both are the types that reached h in the logged
                    run, read off these same starts.

Then, as constrained_generate: rule m = h and the controls h - 1, h + 1,
pairing kept from the anchor, at both trial limits. Observer slots null.

  python3 experiments/anchor_selection.py [seed]
  python3 experiments/anchor_selection.py --small [seed]   # phone bound; own run file
"""
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anchor_generate as A  # noqa: E402
import constrained_generate as K  # noqa: E402
import emlib  # noqa: E402
import mirror_pursuit as U  # noqa: E402
import pairing_reread as R  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "anchor_selection"
PREFERRED = {"5-cell": 1, "16-cell": 2}       # points on the axis; fixed from seed 0
PREFERRED_FROM_SEED = 0                       # a run at another seed is a holdout for them
MAX_ATTEMPTS = 12                             # stage-1 pursuits per start, the first included
LIMITS = (K.MAX_TRIALS, K.BOUND_CHECK_TRIALS)
SMALL = {"starts": 2, "limits": (80,), "anchor_max_evals": 1000, "max_attempts": 6}
RULES = R.RULES
HYPOTHESIS = ("PROPOSED by the operator: classify the anchor by points on the axis, prefer the "
              "type that matched the full symmetry, re-anchor otherwise; more starts reach h "
              "than with the anchor as found.")
KEEP = tuple(k for k in R.KEEP if k not in ("pull_backs_with_a_different_pairing",))


def anchor_of(verts, planes, start, anchor_max_evals):
    anc = U.descend(verts, start, max_evals=anchor_max_evals)
    a = anc["final_reading"]["axis_deg"]
    pair = K.pairing(anc["final_projected"], a)
    return {"start_angles_deg": list(start), "final_angles_deg": anc["final_angles_deg"],
            "final_residual": anc["final_residual"], "axis_deg": a, "pairing": pair,
            "points_on_axis": K.points_on_axis(pair), "evals": anc["evals"], "stop": anc["stop"],
            "counts": A.counts(anc["final_projected"]),
            "degeneracy": Q.degeneracy(anc["final_projected"])}


def run(seed=0, starts=K.STARTS, limits=LIMITS, anchor_max_evals=U.MAX_EVALS,
        max_attempts=MAX_ATTEMPTS, reference_path=None):
    rng = random.Random(seed)
    ref = R.kept_reference(reference_path, seed)
    shapes = []
    for name, make, h in K.SHAPES:
        verts = make()
        planes = P.all_planes(len(verts[0]))
        want = PREFERRED[name]
        rows = []
        for i in range(starts):
            start = [rng.uniform(0.0, 360.0) for _ in planes]
            # fresh angles for re-anchoring come from their own stream, so the
            # eight starts stay the ones every other run used
            extra = random.Random(f"{seed}/{name}/{i}/re-anchor")
            attempts = [anchor_of(verts, planes, start, anchor_max_evals)]
            while attempts[-1]["points_on_axis"] != want and len(attempts) < max_attempts:
                fresh = [extra.uniform(0.0, 360.0) for _ in planes]
                attempts.append(anchor_of(verts, planes, fresh, anchor_max_evals))
            anc = attempts[-1]
            rules = {}
            for label, m in (("h_minus_1", h - 1), ("h", h), ("h_plus_1", h + 1)):
                rules[label] = {}
                for limit in limits:
                    g = K.constrained_descent(verts, planes, anc["final_angles_deg"],
                                              anc["axis_deg"], m, limit)
                    rules[label][str(limit)] = {k: g[k] for k in KEEP}
            rows.append({"start_index": i, "start_angles_deg": start,
                         "first_anchor_points_on_axis": attempts[0]["points_on_axis"],
                         "re_anchored": len(attempts) > 1,
                         "attempts": len(attempts),
                         "attempt_types": [a["points_on_axis"] for a in attempts],
                         "anchor_evals_total": sum(a["evals"] for a in attempts),
                         "selected_anchor_is_preferred_type": anc["points_on_axis"] == want,
                         "attempt_log": [{k: a[k] for k in ("start_angles_deg", "points_on_axis",
                                                            "final_residual", "evals", "stop")}
                                         for a in attempts],
                         "anchor": anc, "rules": rules})
        r = ref.get(name)
        re_idx = [x["start_index"] for x in rows if x["re_anchored"]]
        compare = {"target_h": h, "count_rel_tol": K.COUNT_TOLS[0],
                   "preferred_points_on_axis": want,
                   "starts_equal_constrained_generate_starts":
                       (r["starts"][:starts] == [x["start_angles_deg"] for x in rows])
                       if r else None,
                   "first_anchor_points_on_axis":
                       [x["first_anchor_points_on_axis"] for x in rows],
                   "selected_anchor_points_on_axis": [x["anchor"]["points_on_axis"] for x in rows],
                   "re_anchored_starts": re_idx,
                   "attempts": [x["attempts"] for x in rows],
                   "all_attempt_types": sorted(t for x in rows for t in x["attempt_types"]),
                   "starts_left_on_another_type":
                       [x["start_index"] for x in rows
                        if not x["selected_anchor_is_preferred_type"]],
                   "anchor_evals_total": [x["anchor_evals_total"] for x in rows],
                   "selected": {}, "random": (r["kept"] if r else None)}
        for limit in limits:
            per = {}
            for label in RULES:
                gs = [x["rules"][label][str(limit)] for x in rows]
                cnt = [g["counts"][0]["n_mirror_axes"] for g in gs]
                deg = [g["degeneracy"]["tags"]["1e-06"]["degenerate"] for g in gs]
                per[label] = {
                    "final_counts": cnt,
                    "final_counts_at_1e-4": [g["counts"][1]["n_mirror_axes"] for g in gs],
                    "degenerate": deg,
                    "starts_reaching_h": sum(1 for c in cnt if c == h),
                    "starts_reaching_h_excluding_degenerate":
                        sum(1 for c, d in zip(cnt, deg) if c == h and not d),
                    "re_anchored_starts_reaching_h":
                        sum(1 for x, c in zip(rows, cnt) if x["re_anchored"] and c == h),
                    "starts_holding_anchor": sum(1 for g in gs if g["anchor_held"]),
                    "final_scores": [g["final_score"] for g in gs],
                    "stops": [g["stop"] for g in gs]}
            compare["selected"][str(limit)] = per
        shapes.append({"shape": name, "h": h, "n_vertices": len(verts),
                       "compare": compare, "starts": rows})
    return {"hypothesis": HYPOTHESIS,
            "bound": {"starts_per_shape": starts, "seed": seed, "trial_limits": list(limits),
                      "preferred_points_on_axis": PREFERRED, "max_attempts": max_attempts,
                      "anchor_max_evals": anchor_max_evals, "anchor_tol": K.ANCHOR_TOL,
                      "count_rel_tols": list(K.COUNT_TOLS),
                      "re_anchor": "stage-1 pursuit again from fresh uniform random angles",
                      "preferred_from_seed": PREFERRED_FROM_SEED,
                      "holdout": seed != PREFERRED_FROM_SEED,
                      "preference_source": (
                          "read off the logged constrained_generate run on these same starts; "
                          "starts that are not re-anchored repeat that run's descents"
                          if seed == PREFERRED_FROM_SEED else
                          f"fixed from seed {PREFERRED_FROM_SEED}; these starts were not used "
                          "to choose them, and nothing was retuned"),
                      "selection_costs_more": "a re-anchored start uses several stage-1 "
                                              "pursuits; anchor_evals_total logs them",
                      "pairing": "kept from the selected anchor",
                      "paths": "not kept for these descents"},
            "shapes": shapes}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["starts"], SMALL["limits"], SMALL["anchor_max_evals"],
               SMALL["max_attempts"]) if small else run(seed))
    limits = out["bound"]["trial_limits"]
    symbols = [str(min(c, 9)) for s in out["shapes"] for label in RULES
               for c in s["compare"]["selected"][str(limits[-1])][label]["final_counts"]]
    comp = emlib.compressibility(
        symbols, "final mirror-axis count per start at the largest trial limit (digit, capped "
                 "at 9): per shape, rule h-1 then h then h+1")
    angles = {s["shape"]: {label: [r["rules"][label][str(limits[-1])]["final_angles_deg"]
                                   for r in s["starts"]] for label in RULES}
              for s in out["shapes"]}
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "shapes": [[n, h] for n, _, h in K.SHAPES],
                         "starts": out["bound"]["starts_per_shape"], "trial_limits": limits,
                         "preferred_points_on_axis": PREFERRED,
                         "max_attempts": out["bound"]["max_attempts"],
                         "rules": ["h - 1", "h", "h + 1"],
                         "reference": R.REFERENCE.format(seed=seed),
                         "preferred_from_seed": PREFERRED_FROM_SEED,
                         "holdout": seed != PREFERRED_FROM_SEED,
                         "hypothesis": HYPOTHESIS, "status_before_run": "PROPOSED"},
                        out, angles, None, comp, None, gz=True)
    print(p)
    for s in out["shapes"]:
        c = s["compare"]
        print(s["shape"], "h =", s["h"], "preferred", c["preferred_points_on_axis"],
              "| first types", c["first_anchor_points_on_axis"], "-> selected",
              c["selected_anchor_points_on_axis"], "| attempts", c["attempts"],
              "| all attempt types", c["all_attempt_types"])
        for limit in limits:
            for label in RULES:
                g = c["selected"][str(limit)][label]
                k = ((c["random"] or {}).get(str(limit)) or {}).get(label)
                print(f"  {limit:5d} {label:9s} selected counts {g['final_counts']} (1e-4: "
                      f"{g['final_counts_at_1e-4']}) reach h {g['starts_reaching_h']} "
                      f"(non-degenerate {g['starts_reaching_h_excluding_degenerate']}; "
                      f"re-anchored reaching {g['re_anchored_starts_reaching_h']})"
                      + (f"   | random counts {k['final_counts']} reach h "
                         f"{k['starts_reaching_h']}" if k else ""))


if __name__ == "__main__":
    main(sys.argv[1:])
