"""r. constrained_generate -- anchor, then generate, with the anchor HELD.

CC0. Co-authored. Stdlib only.

experiments/anchor_generate.py generated axes from an anchor but left the
rotation free to give the anchor up. Here the anchor mirror is held: every
accepted position has anchor residual below ANCHOR_TOL.

  anchor      as before: descent on probes.mirror_residual from the same
              starts; a is the mirror's axis, and the pairing of each point
              with its mirror image is read off there and kept.
  constraint  c(angles) = the vector of (reflection of point i in axis a)
              minus (its paired point). c = 0 is an exact mirror on axis a.
  projection step
              after every trial move the six angles are pulled back onto
              c = 0 by damped Gauss-Newton steps (finite-difference
              Jacobian, 6 x 6 normal equations). A trial whose anchor
              residual is still >= ANCHOR_TOL afterwards is refused.
  generate    axes a + k * 180/m, k = 1..m-1 (the anchor, k = 0, is held,
              not scored). Score = sum of probes.mirror_residual_at on
              them. Coordinate moves, +-step, first improvement kept, step
              halved after a pass with no accepted move.
                m = h, and the controls m = h - 1, m = h + 1.

Same shapes and starts as mirror_pursuit, layered_pursuit and
anchor_generate (seed 0, 8 per shape). Every final projection is read with
probes.degeneracy. Observer slots are written null.

The full run takes about a minute here (the 2500-trial bound check); the
--small option is the phone bound.

  python3 experiments/constrained_generate.py [seed]
  python3 experiments/constrained_generate.py --small [seed]   # phone bound: fewer starts
                                                               # and trials, every output
                                                               # kept, own run file
"""
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import anchor_generate as A  # noqa: E402
import emlib  # noqa: E402
import mirror_pursuit as U  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "constrained_generate"
SHAPES = A.SHAPES
STARTS = 8
ANCHOR_TOL = 1e-6           # accepted positions keep the anchor residual below this
STEP0 = 30.0                # first step, degrees
STEP_MIN = 1e-6             # stop when the step is below this, degrees
SCORE_STOP = 1e-11
MAX_TRIALS = 400            # trial moves per descent (each is projected first)
BOUND_CHECK_TRIALS = 2500   # the same descent again at a higher limit, logged beside it
GN_ITERS = 8                # Gauss-Newton steps per projection, at most
GN_STOP = 1e-10             # stop projecting when max |c| is below this
FD_DEG = 1e-5               # finite-difference step for the Jacobian, degrees
DAMPING = 1e-12             # added to the normal equations' diagonal (relative to its trace)
COUNT_TOLS = A.COUNT_TOLS
SMALL = {"starts": 2, "max_trials": 80, "anchor_max_evals": 1000, "bound_check_trials": None}
REFERENCE = "runs/anchor_generate_seed0.json"


def project(verts, planes, angles_deg):
    return C.rotate_project(verts, [(pl, math.radians(a)) for pl, a in zip(planes, angles_deg)])


def pairing(pts, axis_deg):
    """For each point, the index of the point nearest its reflection in the axis."""
    a = math.radians(axis_deg)
    c2, s2 = math.cos(2 * a), math.sin(2 * a)
    cx, cy = P.centroid(pts)
    out = []
    for x, y in pts:
        x, y = x - cx, y - cy
        rx, ry = c2 * x + s2 * y, s2 * x - c2 * y
        out.append(min(range(len(pts)),
                       key=lambda j: (rx - (pts[j][0] - cx)) ** 2 + (ry - (pts[j][1] - cy)) ** 2))
    return out


def constraint(pts, axis_deg, pair):
    a = math.radians(axis_deg)
    c2, s2 = math.cos(2 * a), math.sin(2 * a)
    cx, cy = P.centroid(pts)
    out = []
    for (x, y), j in zip(pts, pair):
        x, y = x - cx, y - cy
        out.append(c2 * x + s2 * y - (pts[j][0] - cx))
        out.append(s2 * x - c2 * y - (pts[j][1] - cy))
    return out


def solve(mat, rhs):
    """Gaussian elimination with partial pivoting; None if singular."""
    n = len(rhs)
    m = [row[:] + [b] for row, b in zip(mat, rhs)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(m[r][i]))
        if abs(m[piv][i]) < 1e-300:
            return None
        m[i], m[piv] = m[piv], m[i]
        for r in range(i + 1, n):
            f = m[r][i] / m[i][i]
            for k in range(i, n + 1):
                m[r][k] -= f * m[i][k]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (m[i][n] - sum(m[i][k] * x[k] for k in range(i + 1, n))) / m[i][i]
    return x


def points_on_axis(pair):
    """How many points a pairing maps to themselves: the points on the mirror line."""
    return sum(1 for i, j in enumerate(pair) if i == j)


def pull_back(verts, planes, angles, axis_deg, pair, reread=False):
    """Gauss-Newton onto c = 0. Returns (angles, projections evaluated), and with
    reread=True also the pairing used: read afresh at the trial position (each
    point paired with the point nearest its reflection there) instead of the
    one handed in."""
    ang = list(angles)
    evals = 0
    if reread:
        pair = pairing(project(verts, planes, ang), axis_deg)
        evals += 1
    for _ in range(GN_ITERS):
        c0 = constraint(project(verts, planes, ang), axis_deg, pair)
        evals += 1
        if max(abs(v) for v in c0) < GN_STOP:
            break
        cols = []
        for i in range(len(ang)):
            t = list(ang)
            t[i] += FD_DEG
            ci = constraint(project(verts, planes, t), axis_deg, pair)
            evals += 1
            cols.append([(b - a) / FD_DEG for a, b in zip(c0, ci)])
        n = len(ang)
        jtj = [[sum(u * v for u, v in zip(cols[i], cols[j])) for j in range(n)] for i in range(n)]
        lam = DAMPING * (sum(jtj[i][i] for i in range(n)) or 1.0)
        for i in range(n):
            jtj[i][i] += lam
        step = solve(jtj, [-sum(u * v for u, v in zip(cols[i], c0)) for i in range(n)])
        if step is None:
            break
        ang = [a + d for a, d in zip(ang, step)]
    if reread:
        return ang, evals, pair
    return ang, evals


def constrained_descent(verts, planes, anchor_angles, axis_deg, m, max_trials=MAX_TRIALS,
                        reread_pairing=False):
    """reread_pairing=False: the anchor's pairing is kept for every pull-back.
    reread_pairing=True: every pull-back reads the pairing afresh at its trial
    position, so the anchor mirror may be held under a different pairing."""
    axes = A.generated_axes(axis_deg, m)
    scored = axes[1:]
    pair = pairing(project(verts, planes, anchor_angles), axis_deg)
    anchor_pair = list(pair)
    current_pair = list(pair)
    rereads_differing = 0
    proj_evals = 0

    def anchor_res(a):
        return Q.mirror_residual_at(project(verts, planes, a), [axis_deg])[0]

    def score(a):
        return sum(Q.mirror_residual_at(project(verts, planes, a), scored))

    ang, e = pull_back(verts, planes, anchor_angles, axis_deg, pair)
    proj_evals += e
    start_anchor = {"before_pull_back": anchor_res(anchor_angles),
                    "after_pull_back": anchor_res(ang)}
    feasible_start = start_anchor["after_pull_back"] < ANCHOR_TOL
    if not feasible_start:
        ang = list(anchor_angles)
    best = score(ang)
    step, n_pass, trials, refused, accepted = STEP0, 0, 0, 0, 0
    path = [{"pass": 0, "step_deg": step, "angles_deg": list(ang), "score": best,
             "anchor_residual": anchor_res(ang)}]
    stop = "step below step_min"
    while step >= STEP_MIN:
        if best < SCORE_STOP:
            stop = "score below score_stop"
            break
        if trials >= max_trials:
            stop = "max_trials reached"
            break
        moved = False
        for i in range(len(ang)):
            for d in (step, -step):
                t = list(ang)
                t[i] += d
                if reread_pairing:
                    t, e, used = pull_back(verts, planes, t, axis_deg, pair, reread=True)
                    if used != anchor_pair:
                        rereads_differing += 1
                else:
                    t, e = pull_back(verts, planes, t, axis_deg, pair)
                    used = pair
                proj_evals += e
                trials += 1
                if anchor_res(t) >= ANCHOR_TOL:
                    refused += 1
                    continue
                r = score(t)
                if r < best:
                    ang, best, moved = t, r, True
                    current_pair = list(used)
                    accepted += 1
                    break
        n_pass += 1
        path.append({"pass": n_pass, "step_deg": step, "angles_deg": list(ang), "score": best,
                     "anchor_residual": anchor_res(ang)})
        if not moved:
            step /= 2.0
    pts = project(verts, planes, ang)
    res = Q.mirror_residual_at(pts, axes)
    return {"m": m, "generated_axes_deg": axes, "anchor_pairing": anchor_pair,
            "anchor_points_on_axis": points_on_axis(anchor_pair),
            "pairing": "re-read at every pull-back" if reread_pairing else "kept from the anchor",
            "final_pairing": pairing(pts, axis_deg),
            "final_points_on_axis": points_on_axis(pairing(pts, axis_deg)),
            "pairing_of_last_accepted_move": current_pair,
            "pull_backs_with_a_different_pairing": rereads_differing if reread_pairing else None,
            "feasible_at_start": feasible_start, "anchor_residual_at_start": start_anchor,
            "start_score": path[0]["score"], "final_score": best,
            "final_residual_per_axis": res, "final_anchor_residual": res[0],
            "anchor_held": res[0] < ANCHOR_TOL,
            "max_anchor_residual_on_path": max(p["anchor_residual"] for p in path),
            "final_angles_deg": ang, "passes": n_pass, "trials": trials,
            "trials_refused_by_constraint": refused, "moves_accepted": accepted,
            "projection_evals": proj_evals, "stop": stop,
            "counts": A.counts(pts), "degeneracy": Q.degeneracy(pts),
            "final_projected": pts, "path": path}


def reference(path=None):
    try:
        rec = emlib.load_run(Path(path) if path else emlib.ROOT / REFERENCE)
    except FileNotFoundError:
        return {}
    return {s["shape"]: {"starts": [r["start_angles_deg"] for r in s["starts"]],
                         "compare": s["compare"]} for s in rec["outputs"]["shapes"]}


def run(seed=0, starts=STARTS, max_trials=MAX_TRIALS, anchor_max_evals=U.MAX_EVALS,
        reference_path=None, bound_check_trials=BOUND_CHECK_TRIALS):
    rng = random.Random(seed)
    ref = reference(reference_path)
    shapes = []
    for name, make, h in SHAPES:
        verts = make()
        planes = P.all_planes(len(verts[0]))
        rows = []
        for i in range(starts):
            start = [rng.uniform(0.0, 360.0) for _ in planes]
            anc = U.descend(verts, start, max_evals=anchor_max_evals)
            a = anc["final_reading"]["axis_deg"]
            rules = {label: constrained_descent(verts, planes, anc["final_angles_deg"], a, m,
                                                max_trials)
                     for label, m in (("h_minus_1", h - 1), ("h", h), ("h_plus_1", h + 1))}
            if bound_check_trials:
                for label, m in (("h_minus_1", h - 1), ("h", h), ("h_plus_1", h + 1)):
                    g = constrained_descent(verts, planes, anc["final_angles_deg"], a, m,
                                            bound_check_trials)
                    rules[label]["bound_check"] = {
                        k: g[k] for k in ("final_score", "final_residual_per_axis",
                                          "final_anchor_residual", "anchor_held",
                                          "max_anchor_residual_on_path", "final_angles_deg",
                                          "passes", "trials", "trials_refused_by_constraint",
                                          "moves_accepted", "projection_evals", "stop",
                                          "counts", "degeneracy", "final_pairing",
                                          "final_points_on_axis")}
                    rules[label]["bound_check"]["max_trials"] = bound_check_trials
            rows.append({"start_index": i, "start_angles_deg": start,
                         "anchor": {"final_angles_deg": anc["final_angles_deg"],
                                    "final_residual": anc["final_residual"], "axis_deg": a,
                                    "evals": anc["evals"], "stop": anc["stop"],
                                    "counts": A.counts(anc["final_projected"]),
                                    "degeneracy": Q.degeneracy(anc["final_projected"])},
                         "rules": rules})
        r = ref.get(name)

        def summary(label, check=False):
            gs = [x["rules"][label]["bound_check"] if check else x["rules"][label]
                  for x in rows]
            cnt = [g["counts"][0]["n_mirror_axes"] for g in gs]
            deg = [g["degeneracy"]["tags"]["1e-06"]["degenerate"] for g in gs]
            return {"m": rows[0]["rules"][label]["m"], "final_counts": cnt, "degenerate": deg,
                    "anchor_points_on_axis": [x["rules"][label]["anchor_points_on_axis"]
                                              for x in rows],
                    "final_counts_at_1e-4": [g["counts"][1]["n_mirror_axes"] for g in gs],
                    "starts_reaching_h": sum(1 for c in cnt if c == h),
                    "starts_reaching_h_excluding_degenerate":
                        sum(1 for c, d in zip(cnt, deg) if c == h and not d),
                    "starts_holding_anchor": sum(1 for g in gs if g["anchor_held"]),
                    "final_anchor_residuals": [g["final_anchor_residual"] for g in gs],
                    "final_scores": [g["final_score"] for g in gs],
                    "final_residual_per_axis": [g["final_residual_per_axis"] for g in gs],
                    "trials": [g["trials"] for g in gs],
                    "trials_refused_by_constraint": [g["trials_refused_by_constraint"]
                                                     for g in gs],
                    "moves_accepted": [g["moves_accepted"] for g in gs],
                    "projection_evals": [g["projection_evals"] for g in gs],
                    "anchor_evals": [x["anchor"]["evals"] for x in rows],
                    "stops": [g["stop"] for g in gs]}
        shapes.append({
            "shape": name, "h": h, "n_vertices": len(verts),
            "planes": [list(p) for p in planes],
            "compare": {
                "target_h": h, "count_rel_tol": COUNT_TOLS[0],
                "starts_equal_anchor_generate_starts":
                    (r["starts"][:starts] == [x["start_angles_deg"] for x in rows]) if r else None,
                "constrained": {label: summary(label)
                                for label in ("h_minus_1", "h", "h_plus_1")},
                "constrained_bound_check": ({label: summary(label, True)
                                             for label in ("h_minus_1", "h", "h_plus_1")}
                                            if bound_check_trials else None),
                "unconstrained_anchor_generate":
                    ({label: {k: r["compare"]["anchor_generate"][label][k]
                              for k in ("final_counts", "starts_reaching_h",
                                        "starts_keeping_anchor", "final_scores")}
                      for label in ("h_minus_1", "h", "h_plus_1")} if r else None),
                "layered_unanchored": (r["compare"]["layered_unanchored"] if r else None),
                "alt": (r["compare"]["alt"] if r else None)},
            "starts": rows})
    return {"bound": {"starts_per_shape": starts, "seed": seed,
                      "anchor_tol": ANCHOR_TOL, "step0_deg": STEP0, "step_min_deg": STEP_MIN,
                      "score_stop": SCORE_STOP, "max_trials": max_trials,
                      "bound_check_max_trials": bound_check_trials,
                      "gauss_newton_iters": GN_ITERS, "gauss_newton_stop": GN_STOP,
                      "finite_difference_deg": FD_DEG, "damping": DAMPING,
                      "anchor_max_evals": anchor_max_evals,
                      "count_rel_tols": list(COUNT_TOLS),
                      "pairing": "read at the anchor and kept; a mirror on axis a under a "
                                 "different pairing is outside the constraint",
                      "axes": "fixed in the projection plane at a + k*180/m; k = 0 held, "
                              "k = 1..m-1 scored",
                      "moves": "one angle at a time, +-step, pulled back, first improvement kept",
                      "trials_not_like_for_like": "a trial here costs a pull-back (up to "
                          "GN_ITERS * 7 projections) plus one score; an anchor_generate "
                          "evaluation is one score"},
            "shapes": shapes}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["starts"], SMALL["max_trials"], SMALL["anchor_max_evals"],
               bound_check_trials=SMALL["bound_check_trials"]) if small else run(seed))
    b = out["bound"]
    symbols = [str(min(c, 9)) for s in out["shapes"]
               for label in ("h_minus_1", "h", "h_plus_1")
               for c in s["compare"]["constrained"][label]["final_counts"]]
    comp = emlib.compressibility(
        symbols, "final mirror-axis count per start (digit, capped at 9): per shape, "
                 "rule h-1 then h then h+1")
    angles = {s["shape"]: {label: [r["rules"][label]["final_angles_deg"] for r in s["starts"]]
                           for label in ("h_minus_1", "h", "h_plus_1")}
              for s in out["shapes"]}
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "shapes": [[n, h] for n, _, h in SHAPES],
                         "starts": b["starts_per_shape"], "rules": ["h - 1", "h", "h + 1"],
                         "anchor_tol": ANCHOR_TOL, "step0_deg": STEP0, "step_min_deg": STEP_MIN,
                         "score_stop": SCORE_STOP, "max_trials": b["max_trials"],
                         "bound_check_max_trials": b["bound_check_max_trials"],
                         "gauss_newton_iters": GN_ITERS, "finite_difference_deg": FD_DEG,
                         "anchor_max_evals": b["anchor_max_evals"],
                         "count_rel_tols": list(COUNT_TOLS), "reference": REFERENCE,
                         "constraint": "anchor residual below anchor_tol at every accepted "
                                       "position, by a Gauss-Newton pull-back"},
                        out, angles, None, comp, None, gz=True)
    print(p)
    for s in out["shapes"]:
        c = s["compare"]
        print(s["shape"], "h =", s["h"], "same starts:", c["starts_equal_anchor_generate_starts"])
        for label in ("h_minus_1", "h", "h_plus_1"):
            g = c["constrained"][label]
            print(f"  held  {label:9s} m={g['m']}  counts {g['final_counts']}  reach h "
                  f"{g['starts_reaching_h']} (non-degenerate "
                  f"{g['starts_reaching_h_excluding_degenerate']})  anchor held "
                  f"{g['starts_holding_anchor']}  accepted {g['moves_accepted']}  refused "
                  f"{g['trials_refused_by_constraint']}  scores "
                  f"{[float(f'{x:.1e}') for x in g['final_scores']]}")
            k = (c["constrained_bound_check"] or {}).get(label)
            if k:
                print(f"        at {b['bound_check_max_trials']} trials: counts "
                      f"{k['final_counts']} (at 1e-4: {k['final_counts_at_1e-4']})  reach h "
                      f"{k['starts_reaching_h']} (non-degenerate "
                      f"{k['starts_reaching_h_excluding_degenerate']})  anchor held "
                      f"{k['starts_holding_anchor']}  on limit "
                      f"{sum(1 for x in k['stops'] if x == 'max_trials reached')}  scores "
                      f"{[float(f'{x:.1e}') for x in k['final_scores']]}")
            u = c["unconstrained_anchor_generate"]
            if u:
                print(f"  free  {label:9s}      counts {u[label]['final_counts']}  reach h "
                      f"{u[label]['starts_reaching_h']}  anchor kept "
                      f"{u[label]['starts_keeping_anchor']}")


if __name__ == "__main__":
    main(sys.argv[1:])
