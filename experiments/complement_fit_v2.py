"""ad. complement_fit_v2 -- equal point budget, points that respond to the field,
and rotations searched on the sphere.

CC0. Co-authored. Stdlib only. Continues experiments/complement_fit.py.
Status: PROPOSED. Concrete choices are below and in the run's bound.

  1 equal budget   a probe with more points fills more. Here every probe
                   contributes the same number of points, N: at the
                   placement the whole probe chose, its N points are picked
                   one at a time, each the point that leaves the least gap.
                   Probes with fewer than N points are left out of that
                   budget. Singles, and pairs (N points from the two placed
                   probes together).
  2 environment v2 the points RESPOND to the field before fitting: a
                   hexagonal sheet is relaxed under the gradient by the
                   field/ rule (field/hexsheet.relax), and every pattern
                   point is moved by the sheet's displacement where it sits.
                   The baseline responds; a probe is placed (size, turn) and
                   then responds. Gaps are counted in plain units. Rankings
                   are set beside the no-field ones.
  3 sphere         the probe is turned in R^4 (six angles) to sit as far as
                   it can from the baseline's directions, then the deepest
                   void of the two together is read. The joined set is set
                   beside a short catalogue of known vertex sets by a
                   signature that no rotation changes (its sorted pairwise
                   angles). A set that matches nothing is logged as matching
                   nothing.

Observer slots null.

  python3 experiments/complement_fit_v2.py [seed]
  python3 experiments/complement_fit_v2.py --small [seed]   # phone bound; own run file
"""
import itertools
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import collapse_gap as G  # noqa: E402
import complement_fit as F  # noqa: E402
import emlib  # noqa: E402
from field import hexsheet as H  # noqa: E402
from toolkit import complement as CP  # noqa: E402
from toolkit import primitives as P  # noqa: E402

NAME = "complement_fit_v2"
N_THETA, STRIDE, TURN_STEP = 1680, 20, 7
BUDGETS = (3, 5, 8)
BASELINES = F.BASELINES
PAIR_PROBES = F.COMBO_PROBES
SHEET_RADIUS, SHEET_ALPHA, SHEET_SLOPE, SHEET_SWEEPS = 6, 0.05, 1.0, 1500
V2_TURN_STEP = 14               # 3 degrees: every placement recomputes the moved points
SPHERE_PAIRS = (("5-cell", "5-cell"), ("5-cell", "16-cell"), ("5-cell", "tesseract"),
                ("16-cell", "16-cell"), ("16-cell", "tesseract"), ("16-cell", "24-cell"),
                ("tesseract", "tesseract"), ("tesseract", "24-cell"), ("24-cell", "24-cell"),
                ("5-cell", "24-cell"))
ROT_STARTS, ROT_STEP0, ROT_STEP_MIN = 12, 30.0, 1e-4
ROT_JIGGLES = 2500              # for small pairs; large pairs get fewer (see search_effort)
ROT_SMALL_PAIR = 100            # a pair is small if (points of one) x (points of the other) <= this              # random moves of all six angles after the coordinate ascent
VOID_STARTS, VOID_STEPS = 60, 300
MATCH_TOL_DEG = 0.5
SMALL = {"n_theta": 840, "budgets": (3, 5), "pair_probes": ("3-gon", "5-gon", "8-gon"),
         "sheet_radius": 4, "sheet_sweeps": 600, "v2_turn_step": 28,
         "sphere_pairs": (("5-cell", "5-cell"), ("16-cell", "16-cell"), ("16-cell", "tesseract")),
         "rot_starts": 3, "void_starts": 25, "void_steps": 200}
METHOD = ("PROPOSED: fit at an equal point budget; let the points respond to the field before "
          "fitting; search rotations on the sphere and check whether a pair completes a known "
          "polytope.")


# ---- 1 equal point budget -------------------------------------------------------

def sample_points(radii, n_theta, stride):
    return [[(r * math.cos(2 * math.pi * j / n_theta), r * math.sin(2 * math.pi * j / n_theta))
             for j in range(0, n_theta, stride)] for r in radii]


def placed(pattern, scale, turn_deg):
    a = math.radians(turn_deg)
    c, s = math.cos(a), math.sin(a)
    # the probe's field was read at sample angle + turn: the probe is turned by -turn
    return [[scale * (c * x + s * y), scale * (-s * x + c * y)] for x, y in pattern]


def greedy_points(current, candidates, n, samples, w):
    """Choose n of the candidate points one at a time, each the one that leaves
    the smallest mean gap. Returns (gap, chosen indices, new samples)."""
    cur = [row[:] for row in current]
    chosen = []
    dist = [[[math.hypot(sx - px, sy - py) for sx, sy in row] for row in samples]
            for px, py in candidates]
    for _ in range(min(n, len(candidates))):
        best = None
        for k in range(len(candidates)):
            if k in chosen:
                continue
            tot = cnt = 0.0
            for crow, drow, wrow in zip(cur, dist[k], w):
                for cv, dv, wv in zip(crow, drow, wrow):
                    tot += min(cv, dv) * wv
                    cnt += 1
            g = tot / cnt
            if best is None or g < best[0] - 1e-15:
                best = (g, k)
        chosen.append(best[1])
        cur = [[min(cv, dv) for cv, dv in zip(crow, drow)]
               for crow, drow in zip(cur, dist[best[1]])]
    return F.gap(cur, w), chosen, cur


def budget_part(seed, n_theta, budgets, baselines, pair_probes):
    shapes = {n: CP.turned_to_zero(CP.centred_unit(p)) for n, _, p in G.plane_shapes(seed)}
    radii, _ = CP.polar_grid(F.N_RINGS, n_theta)
    fields = {n: {s: CP.distance_field([[s * x, s * y] for x, y in pat], radii, n_theta)
                  for s in F.SCALES} for n, pat in shapes.items()}
    samples = sample_points(radii, n_theta, STRIDE)
    w = F.weights(radii, n_theta, STRIDE, "none")
    out = []
    for b in baselines:
        base = F.sample(fields[b][1.0], n_theta, STRIDE)
        g0 = F.gap(base, w)
        place1 = {}
        full = []
        for p in shapes:
            g, scale, turn, new = F.place(base, fields[p], n_theta, STRIDE, TURN_STEP, w)
            place1[p] = (scale, turn, new)
            full.append({"name": "self" if p == b else p, "shape": p,
                         "n_points": len(shapes[p]), "residual_gap": g})
        full = F.ranking(full)
        per_budget = []
        for n in budgets:
            singles = []
            for p, (scale, turn, _) in place1.items():
                if len(shapes[p]) < n:
                    continue
                g, idx, _ = greedy_points(base, placed(shapes[p], scale, turn), n, samples, w)
                singles.append({"name": "self" if p == b else p, "shape": p,
                                "n_points_available": len(shapes[p]), "n_points_used": n,
                                "residual_gap": g, "gap_closed": g0 - g,
                                "scale": scale, "turn_deg": turn, "points_chosen": idx})
            singles = F.ranking(singles)
            order = [r["shape"] for r in full if r["shape"] in pair_probes]
            pairs = []
            for a, c in itertools.combinations(order, 2):
                sa, ta, after_a = place1[a]
                _, sc, tc, _ = F.place(after_a, fields[c], n_theta, STRIDE, TURN_STEP, w)
                cand = placed(shapes[a], sa, ta) + placed(shapes[c], sc, tc)
                if len(cand) < n:
                    continue
                g, idx, _ = greedy_points(base, cand, n, samples, w)
                from_a = sum(1 for k in idx if k < len(shapes[a]))
                pairs.append({"name": f"{a} + {c}", "probes": [a, c],
                              "n_points_available": len(cand), "n_points_used": n,
                              "points_from_each": [from_a, n - from_a],
                              "residual_gap": g, "gap_closed": g0 - g})
            pairs = F.ranking(pairs)
            best_single = singles[0]["residual_gap"] if singles else None
            per_budget.append({
                "budget": n, "left_out_below_budget": [p for p in shapes if len(shapes[p]) < n],
                "singles": singles, "ranking": [r["name"] for r in singles],
                "pairs": pairs, "pair_ranking": [r["name"] for r in pairs],
                "best_single_gap": best_single,
                "pairs_below_the_best_single":
                    [r["name"] for r in pairs
                     if best_single is not None and r["residual_gap"] < best_single - 1e-12],
                "pairs_using_points_of_both":
                    [r["name"] for r in pairs if 0 < r["points_from_each"][0] < n]})
        out.append({"baseline": b, "baseline_gap": g0,
                    "full_probe_ranking": [r["name"] for r in full],
                    "budgets": per_budget})
    return out


# ---- 2 points respond to the field ------------------------------------------------

def response_map(radius, alpha, slope, sweeps):
    """Relax a hexagonal sheet under the gradient (the field/ rule) and return a
    function moving a plane point by the sheet's displacement there. The
    sheet spans the unit disc: one plane unit = `radius` cell spacings."""
    sheet = H.cells(radius)
    field = H.gradient_field(sheet, slope)
    pos, bound = H.relax(sheet, field, alpha, sweeps)
    cx = sum(p[0] for p in pos.values()) / len(pos)
    cy = sum(p[1] for p in pos.values()) / len(pos)
    ideal = {c: H.ideal_xy(c) for c in sheet}
    disp = [((ideal[c][0] / radius, ideal[c][1] / radius),
             ((pos[c][0] - cx - ideal[c][0]) / radius, (pos[c][1] - cy - ideal[c][1]) / radius))
            for c in sheet]

    def move(p):
        near = sorted(disp, key=lambda d: (d[0][0] - p[0]) ** 2 + (d[0][1] - p[1]) ** 2)[:3]
        ws = [1.0 / (math.hypot(d[0][0] - p[0], d[0][1] - p[1]) + 1e-9) for d in near]
        tot = sum(ws)
        return [p[0] + sum(wv * d[1][0] for wv, d in zip(ws, near)) / tot,
                p[1] + sum(wv * d[1][1] for wv, d in zip(ws, near)) / tot]
    size = [math.hypot(*d[1]) for d in disp]
    return move, {"sheet_radius_cells": radius, "n_cells": len(sheet), "alpha": alpha,
                  "slope": slope, "field_range": [min(field.values()), max(field.values())],
                  "relaxation": bound, "largest_displacement_plane_units": max(size),
                  "mean_displacement_plane_units": sum(size) / len(size),
                  "interpolation": "inverse-distance over the three nearest cells"}


def field_at(points, samples):
    return [[min(math.hypot(sx - px, sy - py) for px, py in points) for sx, sy in row]
            for row in samples]


def place_moving(current, pattern, samples, w, move, turn_step, n_theta):
    """Place a probe whose points respond after it is placed: every size and
    turn is tried, the placed points are moved, and the gap read."""
    best = None
    for scale in F.SCALES:
        for shift in range(0, n_theta, turn_step):
            turn = 360.0 * shift / n_theta
            pts = [move(p) for p in placed(pattern, scale, turn)] if move else \
                placed(pattern, scale, turn)
            f = field_at(pts, samples)
            tot = cnt = 0.0
            for crow, frow, wrow in zip(current, f, w):
                for cv, fv, wv in zip(crow, frow, wrow):
                    tot += min(cv, fv) * wv
                    cnt += 1
            g = tot / cnt
            if best is None or g < best[0] - 1e-15:
                best = (g, scale, turn, f)
    g, scale, turn, f = best
    new = [[min(cv, fv) for cv, fv in zip(crow, frow)] for crow, frow in zip(current, f)]
    return g, scale, turn, new


def response_part(seed, n_theta, baselines, pair_probes, sheet_radius, sheet_sweeps, turn_step):
    shapes = {n: CP.turned_to_zero(CP.centred_unit(p)) for n, _, p in G.plane_shapes(seed)}
    radii, _ = CP.polar_grid(F.N_RINGS, n_theta)
    samples = sample_points(radii, n_theta, STRIDE)
    w = F.weights(radii, n_theta, STRIDE, "none")
    move, bound = response_map(sheet_radius, SHEET_ALPHA, SHEET_SLOPE, sheet_sweeps)
    out = []
    for b in baselines:
        per = {}
        for env, mv in (("no_field", None), ("field_response", move)):
            base_pts = [mv(p) for p in shapes[b]] if mv else shapes[b]
            base = field_at(base_pts, samples)
            g0 = F.gap(base, w)
            singles, cache = [], {}
            for p in shapes:
                g, scale, turn, new = place_moving(base, shapes[p], samples, w, mv, turn_step,
                                                   n_theta)
                cache[p] = new
                singles.append({"name": "self" if p == b else p, "shape": p,
                                "n_points": len(shapes[p]), "residual_gap": g,
                                "gap_closed": g0 - g, "scale": scale, "turn_deg": turn})
            singles = F.ranking(singles)
            order = [r["shape"] for r in singles if r["shape"] in pair_probes]
            pairs = []
            for a, c in itertools.combinations(order, 2):
                g, scale, turn, _ = place_moving(cache[a], shapes[c], samples, w, mv, turn_step,
                                                 n_theta)
                pairs.append({"name": f"{a} + {c}", "probes": [a, c], "residual_gap": g,
                              "gap_closed": g0 - g, "second_scale": scale,
                              "second_turn_deg": turn})
            pairs = F.ranking(pairs)
            per[env] = {"baseline_gap": g0,
                        "baseline_points": base_pts, "singles": singles,
                        "ranking": [r["name"] for r in singles],
                        "pairs": pairs, "pair_ranking": [r["name"] for r in pairs]}

        def shift(key):
            ra, rc = per["no_field"][key], per["field_response"][key]
            pos = {n: i for i, n in enumerate(rc)}
            moves = [{"name": n, "rank_before": i + 1, "rank_after": pos[n] + 1,
                      "moved_by": i - pos[n]} for i, n in enumerate(ra) if n in pos]
            return {"same_order": ra == rc,
                    "n_changed_rank": sum(1 for m in moves if m["moved_by"]),
                    "largest_move": max((abs(m["moved_by"]) for m in moves), default=0),
                    "best_before": ra[0] if ra else None, "best_after": rc[0] if rc else None,
                    "moves": moves}
        out.append({"baseline": b, "environments": per,
                    "shift_singles": shift("ranking"), "shift_pairs": shift("pair_ranking"),
                    "baseline_point_moves": [
                        math.hypot(q[0] - p[0], q[1] - p[1])
                        for p, q in zip(shapes[b], per["field_response"]["baseline_points"])]})
    return {"response_bound": bound, "probe_turn_step_deg": 360.0 * turn_step / n_theta,
            "gaps_counted_in": "plain units (the field acts through the positions only)",
            "baselines": out}


# ---- 3 sphere: rotations searched, and what the pair completes ----------------------

def rotate4(dirs, angles_deg):
    planes = P.all_planes(4)
    return P.rotate_many(dirs, [(pl, math.radians(a)) for pl, a in zip(planes, angles_deg)])


def separation(a_dirs, b_dirs):
    """Smallest angle (degrees) between a direction of one set and one of the other."""
    best = -1.0
    for a in a_dirs:
        for b in b_dirs:
            d = P.dot(a, b)
            if d > best:
                best = d
    return math.degrees(math.acos(max(-1.0, min(1.0, best))))


def search_rotation(base, probe, rng, starts, step0=ROT_STEP0, step_min=ROT_STEP_MIN,
                    jiggles=ROT_JIGGLES):
    """Turn the probe (six angles) to make its smallest angle to the baseline
    as large as possible: from several random starts, coordinate ascent, then
    random moves of all six angles at once with a shrinking size (the
    smallest angle has kinks where one-angle moves stall)."""
    best = None
    evals = 0
    for s in range(starts):
        ang = [0.0] * 6 if s == 0 else [rng.uniform(0.0, 360.0) for _ in range(6)]
        val = separation(base, rotate4(probe, ang))
        evals += 1
        step = step0
        while step >= step_min:
            moved = False
            for i in range(6):
                for d in (step, -step):
                    t = list(ang)
                    t[i] += d
                    v = separation(base, rotate4(probe, t))
                    evals += 1
                    if v > val + 1e-12:
                        ang, val, moved = t, v, True
                        break
            if not moved:
                step /= 2.0
        size = 5.0
        for _ in range(jiggles):
            t = [a + rng.gauss(0.0, size) for a in ang]
            v = separation(base, rotate4(probe, t))
            evals += 1
            if v > val + 1e-12:
                ang, val = t, v
            else:
                size = max(size * 0.99, 1e-4)
        if best is None or val > best[0]:
            best = (val, ang, s)
    return {"separation_deg": best[0], "angles_deg": best[1], "start_index": best[2],
            "evals": evals}


def search_effort(n_base, n_probe, starts):
    """Starts and random moves for a pair: the full effort for small pairs,
    a third of the starts and a quarter of the moves for large ones, so the
    run stays short. Logged with every pair."""
    if n_base * n_probe <= ROT_SMALL_PAIR:
        return starts, ROT_JIGGLES
    return max(1, starts // 3), ROT_JIGGLES // 4


def signature(dirs):
    """Sorted pairwise angles (degrees) of a direction set: unchanged by any
    rotation or reflection of the whole set."""
    out = []
    for i in range(len(dirs)):
        for j in range(i + 1, len(dirs)):
            out.append(math.degrees(math.acos(max(-1.0, min(1.0, P.dot(dirs[i], dirs[j]))))))
    return sorted(out)


def same_signature(a, b, tol=MATCH_TOL_DEG):
    return len(a) == len(b) and all(abs(x - y) <= tol for x, y in zip(a, b))


def catalogue():
    sim = CP.unit_dirs(P.simplex(4))
    cross, cube, c24 = (CP.unit_dirs(P.cross_polytope(4)), CP.unit_dirs(P.hypercube(4)),
                        CP.unit_dirs(P.cell24()))
    e = [[1.0 if i == k else 0.0 for i in range(5)] for k in range(5)]
    rect = [[a + b - 0.4 for a, b in zip(e[i], e[j])] for i, j in itertools.combinations(range(5), 2)]
    basis = P.gram_schmidt([[x - 0.2 for x in v] for v in e[:4]])
    rect4 = CP.unit_dirs([[P.dot(v, bb) for bb in basis] for v in rect])
    return {"5-cell": sim, "16-cell": cross, "tesseract": cube, "24-cell": c24,
            "rectified 5-cell": rect4,
            "5-cell + its negative (as built)": sim + [[-x for x in v] for v in sim],
            "16-cell + tesseract (as built)": cross + cube}


def sphere_part(seed, pairs, rot_starts, void_starts, void_steps):
    cat = catalogue()
    sigs = {k: signature(v) for k, v in cat.items()}
    rng = random.Random(f"{seed}/sphere-rotations")
    names = list(cat)
    cross = [[x, y] for i, x in enumerate(names) for y in names[i + 1:]
             if same_signature(sigs[x], sigs[y])]
    rows = []
    for a, b in pairs:
        base, probe = cat[a], cat[b]
        as_built = separation(base, probe)
        n_starts, n_jig = search_effort(len(base), len(probe), rot_starts)
        s = search_rotation(base, probe, rng, n_starts, jiggles=n_jig)
        turned = rotate4(probe, s["angles_deg"])
        union = base + turned
        v_base = CP.sphere_voids(base, seed, void_starts, void_steps)
        v_union = CP.sphere_voids(union, seed, void_starts, void_steps)
        v_built = (CP.sphere_voids(base + probe, seed, void_starts, void_steps)
                   if as_built > 1.0 else None)
        sig = signature(union)
        matches = [k for k, sg in sigs.items() if same_signature(sig, sg)]
        same_size = [(max(abs(x - y) for x, y in zip(sig, sg)), k)
                     for k, sg in sigs.items() if len(sg) == len(sig)]
        nearest = min(same_size) if same_size else None
        angles = CP.distinct_values(sig, tol=MATCH_TOL_DEG)
        rows.append({"baseline": a, "probe": b, "n_points": len(union),
                     "separation_as_built_deg": as_built,
                     "separation_after_search_deg": s["separation_deg"],
                     "rotation_angles_deg": s["angles_deg"], "rotation_evals": s["evals"],
                     "rotation_starts": n_starts, "rotation_random_moves": n_jig,
                     "baseline_deepest_void_deg": v_base["deepest_depth_deg"],
                     "residual_gap_as_built_deg":
                         v_built["deepest_depth_deg"] if v_built else None,
                     "residual_gap_after_search_deg": v_union["deepest_depth_deg"],
                     "distinct_pairwise_angles_deg": angles,
                     "matches_in_catalogue": matches,
                     "completes_a_known_polytope": bool(matches),
                     "nearest_catalogue_entry_of_this_size":
                         ({"name": nearest[1], "largest_angle_difference_deg": nearest[0]}
                          if nearest else None),
                     "catalogue_entries_with_this_many_points":
                         [k for k, v in cat.items() if len(v) == len(union)],
                     "joined_directions": union})
    return {"catalogue": {k: {"n_points": len(v),
                              "distinct_pairwise_angles_deg":
                                  CP.distinct_values(sigs[k], tol=MATCH_TOL_DEG)}
                          for k, v in cat.items()},
            "catalogue_entries_that_match_each_other": cross,
            "search": "the probe turned to make its smallest angle to the baseline largest; "
                      "that is a stand-in for the residual gap, which is then read",
            "match": f"sorted pairwise angles equal within {MATCH_TOL_DEG} degrees",
            "pairs": rows}


def run(seed=0, n_theta=N_THETA, budgets=BUDGETS, baselines=BASELINES, pair_probes=PAIR_PROBES,
        sheet_radius=SHEET_RADIUS, sheet_sweeps=SHEET_SWEEPS, v2_turn_step=V2_TURN_STEP,
        sphere_pairs=SPHERE_PAIRS, rot_starts=ROT_STARTS, void_starts=VOID_STARTS,
        void_steps=VOID_STEPS):
    return {"method": METHOD,
            "bound": {"seed": seed, "n_theta": n_theta, "samples_per_ring": n_theta // STRIDE,
                      "budgets": list(budgets), "baselines": list(baselines),
                      "pair_probes": list(pair_probes),
                      "budget_choice": "at the whole probe's placement; points picked one at a "
                                       "time, least gap first; not a search over all subsets "
                                       "or all placements",
                      "rot_starts": rot_starts, "rot_step_min_deg": ROT_STEP_MIN,
                      "rot_jiggles": ROT_JIGGLES,
                      "rotation_search_is_a_lower_bound": "a larger separation than the one "
                          "found may exist; a pair that matches nothing may match at a "
                          "placement the search did not reach",
                      "void_starts": void_starts, "void_steps": void_steps,
                      "catalogue_is_short": "a joined set that matches nothing here may still "
                                            "be a known polytope that is not in the catalogue"},
            "equal_budget": budget_part(seed, n_theta, budgets, baselines, pair_probes),
            "field_response": response_part(seed, n_theta, baselines, pair_probes, sheet_radius,
                                            sheet_sweeps, v2_turn_step),
            "sphere": sphere_part(seed, sphere_pairs, rot_starts, void_starts, void_steps)}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["n_theta"], SMALL["budgets"], BASELINES, SMALL["pair_probes"],
               SMALL["sheet_radius"], SMALL["sheet_sweeps"], SMALL["v2_turn_step"],
               SMALL["sphere_pairs"], SMALL["rot_starts"], SMALL["void_starts"],
               SMALL["void_steps"]) if small else run(seed))
    symbols = ["m" if r["completes_a_known_polytope"] else "n" for r in out["sphere"]["pairs"]]
    comp = emlib.compressibility(
        symbols, "one letter per sphere pair: m the joined set matches a catalogue entry, n "
                 "it matches none")
    b = out["bound"]
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "n_theta": b["n_theta"], "budgets": b["budgets"],
                         "baselines": b["baselines"], "pair_probes": b["pair_probes"],
                         "sheet": {"radius": out["field_response"]["response_bound"]
                                   ["sheet_radius_cells"], "alpha": SHEET_ALPHA,
                                   "slope": SHEET_SLOPE},
                         "sphere_pairs": [[r["baseline"], r["probe"]]
                                          for r in out["sphere"]["pairs"]],
                         "method": METHOD, "status_before_run": "PROPOSED"},
                        out, None, None, comp, None, gz=True)
    print(p)
    for x in out["equal_budget"]:
        print("BUDGET baseline", x["baseline"], "| full-probe ranking", x["full_probe_ranking"])
        for bu in x["budgets"]:
            print(f"  N={bu['budget']}: singles", [(r["name"], round(r["residual_gap"], 4))
                                                   for r in bu["singles"]])
            print("        pairs", [(r["name"], round(r["residual_gap"], 4), r["points_from_each"])
                                    for r in bu["pairs"][:4]],
                  "| pairs below best single:", bu["pairs_below_the_best_single"])
    fr = out["field_response"]
    print("RESPONSE sheet: largest displacement",
          round(fr["response_bound"]["largest_displacement_plane_units"], 4))
    for x in fr["baselines"]:
        e = x["environments"]
        print("  baseline", x["baseline"], "gap", round(e["no_field"]["baseline_gap"], 4), "->",
              round(e["field_response"]["baseline_gap"], 4))
        print("    no field :", e["no_field"]["ranking"])
        print("    response :", e["field_response"]["ranking"])
        print("    singles changed rank", x["shift_singles"]["n_changed_rank"], "largest move",
              x["shift_singles"]["largest_move"], "| pairs changed rank",
              x["shift_pairs"]["n_changed_rank"], "largest move", x["shift_pairs"]["largest_move"])
    print("SPHERE  catalogue entries that match each other:",
          out["sphere"]["catalogue_entries_that_match_each_other"])
    print("   catalogue:", {k: (v["n_points"], [(round(a["value"], 1), a["count"])
                                                 for a in v["distinct_pairwise_angles_deg"]])
                            for k, v in out["sphere"]["catalogue"].items()})
    for r in out["sphere"]["pairs"]:
        print(f"  {r['baseline']:9s} + {r['probe']:9s} n={r['n_points']:2d} separation "
              f"{r['separation_as_built_deg']:.1f} -> {r['separation_after_search_deg']:.2f}  "
              f"gap {r['baseline_deepest_void_deg']:.1f} -> {r['residual_gap_after_search_deg']:.1f}"
              f"  angles {[(round(a['value'], 1), a['count']) for a in r['distinct_pairwise_angles_deg']][:6]}"
              f"  matches {r['matches_in_catalogue']}  nearest "
              f"{r['nearest_catalogue_entry_of_this_size']}")


if __name__ == "__main__":
    main(sys.argv[1:])
