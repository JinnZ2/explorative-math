"""u. pairing_reread -- the held anchor, with its pairing read afresh.

CC0. Co-authored. Stdlib only.

experiments/constrained_generate.py holds the anchor mirror under the
pairing read at the anchor (which point is the mirror image of which). A
mirror line through a set fixes the points lying on it and swaps the rest
in pairs, so a pairing has a type: how many points sit on the axis.

Hypothesis, PROPOSED by the operator: the 5-cell is weaker because of the
fixed pairing. With odd h each mirror passes through one point and pairs
the other four; the pairing read at the anchor may not carry. Proposed
change: re-read the pairing.

  kept      the logged constrained_generate run (pairing kept from the anchor)
  re-read   the same descents, every pull-back reading the pairing afresh
            at its trial position; the anchor residual must still end below
            the anchor tolerance, under whatever pairing holds there

Same shapes, starts, rules (m = h, h - 1, h + 1) and trial limits. Logged
per start: the anchor pairing's points on the axis, the final pairing's,
how many pull-backs used a pairing different from the anchor's, and the
counts. Observer slots are written null.

  python3 experiments/pairing_reread.py [seed]
  python3 experiments/pairing_reread.py --small [seed]   # phone bound; own run file
"""
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anchor_generate as A  # noqa: E402
import constrained_generate as K  # noqa: E402
import emlib  # noqa: E402
import mirror_pursuit as U  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "pairing_reread"
LIMITS = (K.MAX_TRIALS, K.BOUND_CHECK_TRIALS)      # 400 and 2500 trials
SMALL = {"starts": 2, "limits": (80,), "anchor_max_evals": 1000}
REFERENCE = "runs/constrained_generate_seed{seed}.json"   # the run of the same seed
RULES = ("h_minus_1", "h", "h_plus_1")
HYPOTHESIS = ("PROPOSED by the operator: the 5-cell is weaker because the pairing read at the "
              "anchor is kept; re-reading the pairing should let more starts reach h.")
KEEP = ("m", "final_score", "final_residual_per_axis", "final_anchor_residual", "anchor_held",
        "max_anchor_residual_on_path", "final_angles_deg", "passes", "trials",
        "trials_refused_by_constraint", "moves_accepted", "projection_evals", "stop", "counts",
        "degeneracy", "anchor_pairing", "anchor_points_on_axis", "final_pairing",
        "final_points_on_axis", "pairing_of_last_accepted_move",
        "pull_backs_with_a_different_pairing", "pairing")


def kept_reference(path=None, seed=0):
    """Counts of the logged kept-pairing run, per shape, rule and trial limit."""
    try:
        rec = emlib.load_run(Path(path) if path else emlib.ROOT / REFERENCE.format(seed=seed))
    except FileNotFoundError:
        return {}
    out = {}
    for s in rec["outputs"]["shapes"]:
        c = s["compare"]
        per = {}
        for key, limit in (("constrained", rec["outputs"]["bound"]["max_trials"]),
                           ("constrained_bound_check",
                            rec["outputs"]["bound"]["bound_check_max_trials"])):
            if c.get(key):
                per[str(limit)] = {label: {"final_counts": c[key][label]["final_counts"],
                                           "final_counts_at_1e-4":
                                               c[key][label]["final_counts_at_1e-4"],
                                           "starts_reaching_h": c[key][label]["starts_reaching_h"]}
                                   for label in RULES}
        out[s["shape"]] = {"starts": [r["start_angles_deg"] for r in s["starts"]], "kept": per}
    return out


def run(seed=0, starts=K.STARTS, limits=LIMITS, anchor_max_evals=U.MAX_EVALS,
        reference_path=None):
    rng = random.Random(seed)
    ref = kept_reference(reference_path, seed)
    shapes = []
    for name, make, h in K.SHAPES:
        verts = make()
        planes = P.all_planes(len(verts[0]))
        rows = []
        for i in range(starts):
            start = [rng.uniform(0.0, 360.0) for _ in planes]
            anc = U.descend(verts, start, max_evals=anchor_max_evals)
            a = anc["final_reading"]["axis_deg"]
            rules = {}
            for label, m in (("h_minus_1", h - 1), ("h", h), ("h_plus_1", h + 1)):
                rules[label] = {}
                for limit in limits:
                    g = K.constrained_descent(verts, planes, anc["final_angles_deg"], a, m,
                                              limit, reread_pairing=True)
                    rules[label][str(limit)] = {k: g[k] for k in KEEP}
            rows.append({"start_index": i, "start_angles_deg": start,
                         "anchor": {"final_angles_deg": anc["final_angles_deg"],
                                    "final_residual": anc["final_residual"], "axis_deg": a,
                                    "points_on_axis":
                                        rules["h"][str(limits[0])]["anchor_points_on_axis"],
                                    "counts": A.counts(anc["final_projected"]),
                                    "degeneracy": Q.degeneracy(anc["final_projected"])},
                         "rules": rules})
        r = ref.get(name)
        compare = {"target_h": h, "count_rel_tol": K.COUNT_TOLS[0],
                   "starts_equal_constrained_generate_starts":
                       (r["starts"][:starts] == [x["start_angles_deg"] for x in rows])
                       if r else None,
                   "anchor_points_on_axis": [x["anchor"]["points_on_axis"] for x in rows],
                   "reread": {}, "kept": (r["kept"] if r else None)}
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
                    "starts_holding_anchor": sum(1 for g in gs if g["anchor_held"]),
                    "final_points_on_axis": [g["final_points_on_axis"] for g in gs],
                    "starts_whose_pairing_type_changed":
                        sum(1 for g in gs
                            if g["final_points_on_axis"] != g["anchor_points_on_axis"]),
                    "pull_backs_with_a_different_pairing":
                        [g["pull_backs_with_a_different_pairing"] for g in gs],
                    "final_scores": [g["final_score"] for g in gs],
                    "stops": [g["stop"] for g in gs]}
            compare["reread"][str(limit)] = per
        shapes.append({"shape": name, "h": h, "n_vertices": len(verts),
                       "compare": compare, "starts": rows})
    return {"hypothesis": HYPOTHESIS,
            "bound": {"starts_per_shape": starts, "seed": seed, "trial_limits": list(limits),
                      "anchor_tol": K.ANCHOR_TOL, "anchor_max_evals": anchor_max_evals,
                      "count_rel_tols": list(K.COUNT_TOLS),
                      "pairing": "each point paired with the point nearest its reflection, "
                                 "read at the trial position before every pull-back",
                      "paths": "not kept for these descents; the kept-pairing run holds its own",
                      "everything_else": "as constrained_generate"},
            "shapes": shapes}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["starts"], SMALL["limits"], SMALL["anchor_max_evals"])
           if small else run(seed))
    limits = out["bound"]["trial_limits"]
    symbols = [str(min(c, 9)) for s in out["shapes"] for label in RULES
               for c in s["compare"]["reread"][str(limits[-1])][label]["final_counts"]]
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
                         "rules": ["h - 1", "h", "h + 1"],
                         "reference": REFERENCE.format(seed=seed),
                         "hypothesis": HYPOTHESIS, "status_before_run": "PROPOSED"},
                        out, angles, None, comp, None, gz=True)
    print(p)
    for s in out["shapes"]:
        c = s["compare"]
        print(s["shape"], "h =", s["h"], "anchor points on axis", c["anchor_points_on_axis"])
        for limit in limits:
            for label in RULES:
                g = c["reread"][str(limit)][label]
                k = ((c["kept"] or {}).get(str(limit)) or {}).get(label)
                print(f"  {limit:5d} {label:9s} re-read counts {g['final_counts']} (1e-4: "
                      f"{g['final_counts_at_1e-4']}) reach h {g['starts_reaching_h']}  type "
                      f"changed {g['starts_whose_pairing_type_changed']}  final on-axis "
                      f"{g['final_points_on_axis']}"
                      + (f"   | kept counts {k['final_counts']} reach h "
                         f"{k['starts_reaching_h']}" if k else ""))


if __name__ == "__main__":
    main(sys.argv[1:])
