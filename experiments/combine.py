"""aa. combine -- fold the shape + anchor by type + rule m = h + weighed multiplicities.

CC0. Co-authored. Stdlib only.

Each piece has been run alone. Here they are run together and one at a
time through ONE loop (the fold_shape loop: fold, pursue, try the candidate
on the whole shape, fold again), on the same starts, so reach and cost can
be set side by side.

  fold      the projected points reflected into the wedge of the mirrors
            found (fold_shape), or not folded
  anchor    as found, or chosen by pairing type (anchor_selection: stage-1
            pursuit repeated until the anchor has the preferred number of
            points on its axis; types fixed from seed 0)
  axes      candidate next axes from pairs of folded points, or only the
            rule's axes a + k * 180 / h
  weighing  a mirror of the folded set must map each folded point to one of
            the SAME multiplicity, each point counting by its multiplicity;
            or multiplicities not weighed

  variants
    fold_only           fold                                  (the fold_shape run again)
    fold_weighted       fold + weighing
    fold_rule           fold + rule
    fold_type           fold + anchor by type
    type_rule_no_fold   anchor by type + rule, no fold
    combined_bare       fold + anchor by type + rule
    combined            fold + anchor by type + rule + weighing

The score at a position is the smallest residual among the candidate axes;
a candidate below SNAP_TRY is tried on the whole shape, as in fold_shape.
The logged rule-alone and type-and-rule runs (constrained_generate,
anchor_selection) are set beside these as references; they descend a sum
over all rule axes for up to 2500 trials and are not like for like.

Every endpoint is read with probes.degeneracy. Observer slots null.

  python3 experiments/combine.py [seed]
  python3 experiments/combine.py --small [seed]   # phone bound; own run file
"""
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anchor_selection as S  # noqa: E402
import constrained_generate as K  # noqa: E402
import emlib  # noqa: E402
import fold_loop as FL  # noqa: E402
import fold_shape as FS  # noqa: E402
import mirror_pursuit as U  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "combine"
SHAPES = FS.SHAPES                      # 5-cell, 16-cell (the shapes with a preferred type)
STARTS = 8
LEVEL_TRIALS = FS.LEVEL_TRIALS
MAX_ATTEMPTS = S.MAX_ATTEMPTS
SMALL = {"starts": 2, "level_trials": 120, "anchor_max_evals": 1000, "max_attempts": 6}
#            name                 anchor       axes     fold   weighed
VARIANTS = (("fold_only",         "as_found", "pairs", True,  False),
            ("fold_weighted",     "as_found", "pairs", True,  True),
            ("fold_rule",         "as_found", "rule",  True,  False),
            ("fold_type",         "by_type",  "pairs", True,  False),
            ("type_rule_no_fold", "by_type",  "rule",  False, False),
            ("combined_bare",     "by_type",  "rule",  True,  False),
            ("combined",          "by_type",  "rule",  True,  True))
REFERENCES = {"fold_shape": "runs/fold_shape_seed{seed}.json",
              "constrained_generate": "runs/constrained_generate_seed{seed}.json",
              "anchor_selection": "runs/anchor_selection_seed{seed}.json"}
ASKED = ("PROPOSED: combine fold shape + anchor by type + rule m = h, and weigh "
         "multiplicities in the folded-set mirror test; compare reach and cost with each alone.")


def residual_on_axis(reps, mult, scale, deg, weighed):
    """Mirror residual of the folded set on one axis through the centre.
    weighed: each point is matched only to points of its own multiplicity
    and counts by its multiplicity."""
    a = math.radians(deg)
    c2, s2 = math.cos(2 * a), math.sin(2 * a)
    total, weight = 0.0, 0.0
    for (x, y), m in zip(reps, mult):
        rx, ry = c2 * x + s2 * y, s2 * x - c2 * y
        if weighed:
            d2 = min((rx - u) ** 2 + (ry - v) ** 2
                     for (u, v), k in zip(reps, mult) if k == m)
            total += m * math.sqrt(d2)
            weight += m
        else:
            total += math.sqrt(min((rx - u) ** 2 + (ry - v) ** 2 for u, v in reps))
            weight += 1
    return total / weight / scale


def make_score(fold, axes_mode, weighed, anchor_axis, h):
    def score(pts, found):
        reps, mult, scale = FS.folded(pts, found, "fold_shape" if fold else "no_fold", None)
        if axes_mode == "rule":
            cands = [(anchor_axis + k * 180.0 / h) % 180.0 for k in range(h)]
            cands = [d for d in cands if all(FL.circ(d, f) > FS.SEP_DEG for f in found)]
            if fold and found:
                a0, width = sorted(found)[0], 180.0 / len(found)
                cands = [d for d in cands if FL.in_wedge(d, a0, width)]
        else:
            th = [math.atan2(y, x) for x, y in reps if math.hypot(x, y) > 1e-12 * scale]
            uniq = {}
            for i, a in enumerate(th):
                for b in th[i:]:
                    d = math.degrees((a + b) / 2.0) % 180.0
                    if all(FL.circ(d, f) > FS.SEP_DEG for f in found):
                        uniq.setdefault(round(d, 9), d)
            cands = list(uniq.values())
        best, best_deg = None, None
        for d in cands:
            r = residual_on_axis(reps, mult, scale, d, weighed)
            if best is None or r < best:
                best, best_deg = r, d
        return best, best_deg, len(cands), len(reps)
    return score


def exact_weighing_check():
    """The regular pentagon folded by one mirror, in exact forms: is the folded
    set's own mirror one that keeps multiplicities?"""
    from exact import Cyc5, ZETA
    from exact import fold as XF
    from exact import inflation as XI
    pts = XI.pentagon()
    f = XF.half_plane_fold(pts, Cyc5(1, 0))
    fp = [z for z, _ in f]
    return {"shape": "regular pentagon folded by one mirror (direction 0)",
            "arithmetic": "exact a + b*phi forms; no tolerance",
            "folded_multiplicities": sorted((m for _, m in f), reverse=True),
            "folded_set_has_a_mirror_at_72_deg": XF.is_mirror(fp, ZETA),
            "that_mirror_keeps_multiplicities": XF.is_mirror_weighed(f, ZETA),
            "the_pentagon_has_that_mirror": XF.is_mirror(pts, ZETA)}


def as_loop_anchor(a):
    """An anchor_selection anchor in the form the fold_shape loop reads."""
    return {"final_reading": {"axis_deg": a["axis_deg"]},
            "final_angles_deg": a["final_angles_deg"]}


def references(seed):
    out = {}
    for key, path in REFERENCES.items():
        try:
            rec = emlib.load_run(emlib.ROOT / path.format(seed=seed))
        except FileNotFoundError:
            out[key] = None
            continue
        per = {}
        for s in rec["outputs"]["shapes"]:
            c = s["compare"]
            if key == "fold_shape":
                g = c["b_fold_shape"]
                per[s["shape"]] = {"final_counts": g["final_counts"],
                                   "starts_reaching_h": g["starts_reaching_h"],
                                   "trials_total": g["trials_total"],
                                   "what": "fold alone (the fold_only variant should repeat it)"}
            elif key == "constrained_generate":
                g = c["constrained_bound_check"]["h"]
                per[s["shape"]] = {"final_counts": g["final_counts"],
                                   "starts_reaching_h": g["starts_reaching_h"],
                                   "trials": g["trials"],
                                   "what": "rule m = h alone, sum over all rule axes, 2500 trials"}
            else:
                g = c["selected"]["2500"]["h"]
                per[s["shape"]] = {"final_counts": g["final_counts"],
                                   "starts_reaching_h": g["starts_reaching_h"],
                                   "anchor_evals_total": c["anchor_evals_total"],
                                   "what": "anchor by type + rule m = h, sum over all rule "
                                           "axes, 2500 trials"}
        out[key] = per
    return out


def run(seed=0, starts=STARTS, level_trials=LEVEL_TRIALS, anchor_max_evals=U.MAX_EVALS,
        max_attempts=MAX_ATTEMPTS):
    rng = random.Random(seed)
    ref = references(seed)
    shapes = []
    for name, make, h in SHAPES:
        verts = make()
        planes = P.all_planes(len(verts[0]))
        want = S.PREFERRED[name]
        rows = []
        for i in range(starts):
            start = [rng.uniform(0.0, 360.0) for _ in planes]
            extra = random.Random(f"{seed}/{name}/{i}/re-anchor")
            attempts = [S.anchor_of(verts, planes, start, anchor_max_evals)]
            while attempts[-1]["points_on_axis"] != want and len(attempts) < max_attempts:
                fresh = [extra.uniform(0.0, 360.0) for _ in planes]
                attempts.append(S.anchor_of(verts, planes, fresh, anchor_max_evals))
            anchors = {"as_found": attempts[0], "by_type": attempts[-1]}
            cost = {"as_found": attempts[0]["evals"], "by_type": sum(a["evals"] for a in attempts)}
            variants = {}
            for vname, amode, axes_mode, fold, weighed in VARIANTS:
                anc = anchors[amode]
                g = FS.loop(verts, planes, as_loop_anchor(anc), h,
                            "fold_shape" if fold else "no_fold", None, level_trials,
                            score_fn=make_score(fold, axes_mode, weighed, anc["axis_deg"], h))
                g.update({"variant": vname, "anchor": amode, "axes": axes_mode, "fold": fold,
                          "multiplicities_weighed": weighed,
                          "anchor_evals": cost[amode],
                          "anchor_points_on_axis": anc["points_on_axis"]})
                variants[vname] = g
            rows.append({"start_index": i, "start_angles_deg": start,
                         "first_anchor_points_on_axis": attempts[0]["points_on_axis"],
                         "attempts": len(attempts), "re_anchored": len(attempts) > 1,
                         "variants": variants})

        def summ(vname):
            gs = [r["variants"][vname] for r in rows]
            trials = [g["trials_total"] for g in gs]
            return {"anchor": gs[0]["anchor"], "axes": gs[0]["axes"], "fold": gs[0]["fold"],
                    "multiplicities_weighed": gs[0]["multiplicities_weighed"],
                    "final_counts": [g["final_count"] for g in gs],
                    "final_counts_at_1e-4": [g["final_count_at_1e-4"] for g in gs],
                    "degenerate": [g["degenerate"] for g in gs],
                    "central": [g["degeneracy"]["tags"]["1e-06"]["central"] for g in gs],
                    "starts_reaching_h": sum(1 for g in gs if g["reaches_h"]),
                    "starts_reaching_h_excluding_degenerate":
                        sum(1 for g in gs if g["reaches_h"] and not g["degenerate"]),
                    "n_degenerate_endpoints": sum(1 for g in gs if g["degenerate"]),
                    "stops": [g["stop"] for g in gs],
                    "final_folded_point_counts": [g["final_folded_point_count"] for g in gs],
                    "final_multiplicities": [g["final_multiplicities"] for g in gs],
                    "axes_per_level": [[lv["axes_found_after"] for lv in g["levels"]] for g in gs],
                    "refused_on_the_shape":
                        [g["folded_mirrors_the_shape_does_not_have"] for g in gs],
                    "trials_total": trials, "trials_sum": sum(trials),
                    "projection_evals_sum": sum(g["projection_evals_total"] for g in gs),
                    "score_evals_sum": sum(g["score_evals_total"] for g in gs),
                    "candidate_axes_sum": sum(g["candidate_axes_total"] for g in gs),
                    "anchor_evals_sum": sum(g["anchor_evals"] for g in gs)}
        per = {v[0]: summ(v[0]) for v in VARIANTS}
        fref = (ref.get("fold_shape") or {}).get(name)
        shapes.append({
            "shape": name, "h": h, "n_vertices": len(verts),
            "compare": {"target_h": h, "count_rel_tol": FS.COUNT_TOL,
                        "preferred_points_on_axis": want,
                        "first_anchor_points_on_axis":
                            [r["first_anchor_points_on_axis"] for r in rows],
                        "re_anchored_starts": [r["start_index"] for r in rows if r["re_anchored"]],
                        "variants": per,
                        "fold_only_repeats_the_fold_shape_run":
                            (per["fold_only"]["final_counts"] == fref["final_counts"][:starts])
                            if fref else None,
                        "references": {k: ((v or {}).get(name)) for k, v in ref.items()}},
            "starts": rows})
    return {"asked": ASKED,
            "bound": {"starts_per_shape": starts, "seed": seed, "level_trials": level_trials,
                      "max_levels": FS.MAX_LEVELS, "snap_try": FS.SNAP_TRY,
                      "sep_deg": FS.SEP_DEG, "merge_tol": FS.MERGE_TOL,
                      "anchor_tol": FS.ANCHOR_TOL, "count_rel_tol": FS.COUNT_TOL,
                      "anchor_max_evals": anchor_max_evals, "max_attempts": max_attempts,
                      "preferred_points_on_axis": S.PREFERRED,
                      "preferred_from_seed": S.PREFERRED_FROM_SEED,
                      "holdout": seed != S.PREFERRED_FROM_SEED,
                      "score": "smallest residual among candidate axes (not a sum)",
                      "rule_axes": "a + k*180/h, not found, and inside the wedge when folded",
                      "weighing": "match within equal multiplicity; weight by multiplicity",
                      "shapes_not_run": ["24-cell"],
                      "shapes_not_run_reason": "no preferred anchor type has been read for it",
                      "references_not_like_for_like": "constrained_generate and anchor_selection "
                          "descend a sum over all rule axes for up to 2500 trials"},
            "exact_weighing_check": exact_weighing_check(),
            "variant_table": [{"variant": v[0], "anchor": v[1], "axes": v[2], "fold": v[3],
                               "multiplicities_weighed": v[4]} for v in VARIANTS],
            "shapes": shapes}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["starts"], SMALL["level_trials"], SMALL["anchor_max_evals"],
               SMALL["max_attempts"]) if small else run(seed))
    symbols = [str(min(c, 9)) for s in out["shapes"] for v in VARIANTS
               for c in s["compare"]["variants"][v[0]]["final_counts"]]
    comp = emlib.compressibility(
        symbols, "final mirror-axis count per start (digit, capped at 9): per shape, variants "
                 "in table order")
    angles = {s["shape"]: {v[0]: [r["variants"][v[0]]["final_mirror_axes_deg"]
                                  for r in s["starts"]] for v in VARIANTS}
              for s in out["shapes"]}
    b = out["bound"]
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "shapes": [[n, h] for n, _, h in SHAPES],
                         "starts": b["starts_per_shape"], "level_trials": b["level_trials"],
                         "variants": [v[0] for v in VARIANTS],
                         "preferred_points_on_axis": S.PREFERRED,
                         "holdout": b["holdout"],
                         "references": {k: v.format(seed=seed) for k, v in REFERENCES.items()},
                         "asked": ASKED, "status_before_run": "PROPOSED"},
                        out, angles, None, comp, None, gz=True)
    print(p)
    print("exact:", out["exact_weighing_check"])
    for s in out["shapes"]:
        c = s["compare"]
        print(s["shape"], "h =", s["h"], "| fold_only repeats fold_shape:",
              c["fold_only_repeats_the_fold_shape_run"], "| re-anchored", c["re_anchored_starts"])
        for v in VARIANTS:
            g = c["variants"][v[0]]
            print(f"  {v[0]:18s} counts {g['final_counts']} reach {g['starts_reaching_h']} "
                  f"(non-degenerate {g['starts_reaching_h_excluding_degenerate']}; degenerate "
                  f"{g['n_degenerate_endpoints']})  trials {g['trials_sum']}  cand axes "
                  f"{g['candidate_axes_sum']}  anchor evals {g['anchor_evals_sum']}  refused "
                  f"{sum(g['refused_on_the_shape'])}")
        for k, r in c["references"].items():
            if r:
                print(f"  ref {k}: counts {r['final_counts']} reach {r['starts_reaching_h']}")


if __name__ == "__main__":
    main(sys.argv[1:])
