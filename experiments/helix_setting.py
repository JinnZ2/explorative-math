"""ar. helix_setting -- the closed chains sorted by Coxeter's setting, eq 5.1
tested where it applies, the m = 4 search in 5 / 6 / 7 dimensions, and two
readings of the four unmatched chains.

CC0. Co-authored. Stdlib only. Status: PROPOSED.

Carried from REFERENCES/related_studies.md R5 (body read by extraction by a
co-author, parts OCR-garbled; not fetched here):

  eq 5.1   tan(n theta) = (1 + sec(2 pi / w)) tan(theta)
           n = the simplex's vertex count, 2 pi / w = its dihedral angle,
           theta = half a plane's turn; the simplex has as many vertices as
           the space its sphere sits in has dimensions (n = sphere dim + 1).

  2  SETTING_IN   chains ordered as in the setting (tetrahedra in R^4,
                  5-simplices in R^6). w is taken from each chain's OWN
                  dihedral angle: for n unit vectors with pairwise cosine c
                  the spherical simplex's dihedral angle d has
                  cos d = c / (1 + (n - 2) c)   (DERIVED: inverse of the Gram
                  matrix (1 - c) I + c J; checked against facet normals in
                  the tests). Residual per plane, in a form that stays
                  bounded:  cos d sin(n th) cos th - (1 + cos d) sin th cos(n th).
  3  SETTING_OUTSIDE   the 4-simplex chains logged in R^6, the four
                  unmatched ones among them: tagged, no residual. The span
                  of each chain's points is logged beside its tag.
  4  slack        the closed-chain search for 4-simplices again, with the
                  rotation in R^5 (two planes and a line turned by 180
                  degrees), R^6 (three planes) and R^7 (three planes and
                  such a line). Chains found in R^5 are in the setting and
                  get the item-2 residual.
  5  projection   Coxeter groups of rank 7 to 30 with h = N, A_n left out,
                  exponents computed from the diagrams: are the chain's
                  turns, or their complements, all exponents?
  6  vertex set   Gram matrix of each unmatched chain's N vertices: number
                  of distinct inner products, antipodal pairs yes / no.

Nothing is interpreted. Observer slots null.

  python3 experiments/helix_setting.py [seed]
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import coxeter_exponents as CE  # noqa: E402
import emlib  # noqa: E402
import spherical_helix as SH  # noqa: E402

NAME = "helix_setting"
TOL = 1e-9
N_MAX = 40
RANKS = (7, 30)


def dihedral(c, n):
    """Dihedral angle (degrees) of the spherical simplex on n unit vectors with
    pairwise cosine c, and w = 360 / angle (None when the angle is 0)."""
    den = 1.0 + (n - 2) * c
    if abs(den) < 1e-12:
        return {"cos": None, "deg": None, "w": None}
    cd = max(-1.0, min(1.0, c / den))
    deg = math.degrees(math.acos(cd))
    return {"cos": cd, "deg": deg, "w": 360.0 / deg if deg > 1e-12 else None}


def residual(theta_deg, n, cos_d):
    t = math.radians(theta_deg)
    return cos_d * math.sin(n * t) * math.cos(t) - (1.0 + cos_d) * math.sin(t) * math.cos(n * t)


def eq51(n_close, turns, n, edge_over_radius):
    c = 1.0 - edge_over_radius ** 2 / 2.0
    d = dihedral(c, n)
    planes = []
    for p in turns:
        th = 180.0 * p / n_close
        r = residual(th, n, d["cos"]) if d["cos"] is not None else None
        planes.append({"turn": p, "turn_deg": 360.0 * p / n_close, "theta_deg": th,
                       "theta_is_90_deg": 2 * p == n_close,
                       "residual_bounded_form": r,
                       "holds": None if r is None else abs(r) < TOL})
    w = d["w"]
    return {"pairwise_cosine": c, "dihedral_cos": d["cos"], "dihedral_deg": d["deg"], "w": w,
            "w_is_an_integer_within_1e-9": None if w is None else abs(w - round(w)) < 1e-9,
            "one_plus_sec": None if d["cos"] is None or abs(d["cos"]) < 1e-12
            else 1.0 + 1.0 / d["cos"],
            "planes": planes,
            "largest_abs_residual": None if d["cos"] is None
            else max(abs(x["residual_bounded_form"]) for x in planes),
            "holds_on_every_plane": None if d["cos"] is None
            else all(x["holds"] for x in planes),
            "planes_holding_not_at_theta_90": sum(
                1 for x in planes if x["holds"] and not x["theta_is_90_deg"])}


def span_of(n_close, turns, shares):
    return SH.rank(SH.orbit(n_close, turns, shares))


def solve(a, b):
    """Gaussian elimination with partial pivoting; None when singular."""
    k = len(a)
    m = [row[:] + [b[i]] for i, row in enumerate(a)]
    for col in range(k):
        piv = max(range(col, k), key=lambda r: abs(m[r][col]))
        if abs(m[piv][col]) < 1e-11:
            return None
        m[col], m[piv] = m[piv], m[col]
        for r in range(k):
            if r != col:
                f = m[r][col] / m[col][col]
                m[r] = [x - f * y for x, y in zip(m[r], m[col])]
    return [m[i][k] / m[i][i] for i in range(k)]


def search(n_planes_free, with_line, n_max=N_MAX, m=4):
    """Closed chains of regular m-simplices on the orbit of a rotation with
    n_planes_free planes turned by 360 p / N (p < N / 2) and, when with_line,
    one line turned by 180 degrees. Squared distance at step gap j is
    2 * sum_i share_i * (1 - cos(2 pi j p_i / N)); gaps 1 .. m must agree."""
    import itertools
    found, tried, singular, zero_share, negative = [], 0, 0, 0, 0
    supports = {}
    for n in range(m + 1, n_max + 1):
        if with_line and n % 2:
            continue
        top = (n - 1) // 2
        for combo in itertools.combinations(range(1, top + 1), n_planes_free):
            turns = list(combo) + ([n // 2] if with_line else [])
            g = n
            for p in turns:
                g = math.gcd(g, p)
            if g != 1:
                continue
            tried += 1
            f = [[1.0 - math.cos(2 * math.pi * j * p / n) for j in range(m + 1)] for p in turns]
            k = len(turns)
            last = f[-1]
            # unknowns: shares 0 .. k-2 ; share k-1 = 1 - sum
            rows = [[(f[i][j] - last[j]) - (f[i][1] - last[1]) for i in range(k - 1)]
                    for j in range(2, m + 1)]
            rhs = [last[1] - last[j] for j in range(2, m + 1)]
            if k - 1 == len(rows):
                sol = solve(rows, rhs)
                if sol is None:
                    singular += 1
                    continue
            else:                       # one more equation than unknowns: solve, then check
                sol = solve(rows[:k - 1], rhs[:k - 1])
                if sol is None:
                    singular += 1
                    continue
            shares = sol + [1.0 - sum(sol)]
            d = [sum(s * fi[j] for s, fi in zip(shares, f)) for j in range(m + 1)]
            if d[1] <= 1e-12 or any(abs(d[j] - d[1]) > 1e-9 * d[1] for j in range(2, m + 1)):
                continue
            if min(shares) < -1e-9:
                negative += 1
                continue
            if min(shares) <= 1e-9:
                zero_share += 1
                # the chain lives on the turns whose share is not zero
                key = (n, tuple(p for p, t in zip(turns, shares) if t > 1e-9))
                supports[key] = supports.get(key, 0) + 1
                continue
            found.append({"n": n, "turns": turns, "turn_deg": [360.0 * p / n for p in turns],
                          "radius_squared_shares": shares,
                          "edge_over_radius": math.sqrt(2 * d[1]),
                          "next_gap_also_equal": abs(
                              sum(s * (1.0 - math.cos(2 * math.pi * (m + 1) * p / n))
                                  for s, p in zip(shares, turns)) - d[1]) < 1e-9 * d[1]})
    return {"dimension": 2 * n_planes_free + (1 if with_line else 0),
            "rotation": f"{n_planes_free} planes" + (" and a line turned by 180 degrees"
                                                    if with_line else ""),
            "searched": f"N from {m + 1} to {n_max}" + (", N even" if with_line else "")
                        + "; plane turns 360 p / N, 1 <= p < N / 2, distinct, all turns "
                          "together coprime with N",
            "unknown_shares": len(found[0]["turns"]) - 1 if found else None,
            "equations": m - 1, "turn_sets_tried": tried,
            "singular_systems_skipped": singular,
            "solutions_with_a_negative_share": negative,
            "solutions_with_a_zero_share": zero_share,
            "zero_share_solutions_by_the_turns_they_live_on": [
                {"n": k[0], "turns_with_a_share": list(k[1]), "times_met": v}
                for k, v in sorted(supports.items())],
            "n_found": len(found), "n_values": sorted({c["n"] for c in found}),
            "listed": found}


def projection(n_close, turns, groups):
    half = n_close / 2.0
    comp_180 = [half - p for p in turns]
    comp_360 = [n_close - p for p in turns]
    rows = []
    for g in groups:
        if g["h"] != n_close:
            continue
        e = set(g["exponents"])
        rows.append({"group": g["group"], "rank": g["rank"], "h": g["h"],
                     "exponents": g["exponents"],
                     "turns_all_exponents": all(p in e for p in turns),
                     "complements_to_half_turn_all_exponents":
                         all(c == int(c) and int(c) in e for c in comp_180),
                     "complements_to_full_turn_all_exponents": all(c in e for c in comp_360),
                     "turns_not_exponents": [p for p in turns if p not in e]})
    hit = [r["group"] for r in rows if r["turns_all_exponents"]
           or r["complements_to_half_turn_all_exponents"]
           or r["complements_to_full_turn_all_exponents"]]
    return {"n": n_close, "turns": list(turns),
            "complements_to_half_turn": comp_180, "complements_to_full_turn": comp_360,
            "groups_with_h_equal_n": rows, "groups_matching": hit,
            "verdict": hit if hit else "none"}


def vertex_set(n_close, turns, shares):
    pts = SH.orbit(n_close, turns, shares)
    vals = sorted({round(sum(a * b for a, b in zip(pts[0], pts[j])), 9)
                   for j in range(1, n_close)})
    # the rotation carries vertex 0 to every vertex, so row 0 holds every value
    full = sorted({round(sum(a * b for a, b in zip(pts[i], pts[j])), 9)
                   for i in range(n_close) for j in range(i + 1, n_close)})
    anti = sum(1 for i in range(n_close) for j in range(i + 1, n_close)
               if max(abs(a + b) for a, b in zip(pts[i], pts[j])) < 1e-9)
    return {"n_vertices": n_close, "n_distinct": len({tuple(round(x, 9) for x in p) for p in pts}),
            "gram_entries_off_diagonal": n_close * (n_close - 1) // 2,
            "n_distinct_inner_products": len(full),
            "row_0_holds_every_value": vals == full,
            "distinct_inner_products": full,
            "contains_minus_1": -1.0 in full,
            "antipodal_pairs": anti, "has_antipodal_pairs": anti > 0,
            "every_vertex_has_an_antipode": anti * 2 == n_close,
            "rounding": 1e-9}


def coxeter_groups():
    out = []
    for kind in ("B", "D"):
        for n in range(RANKS[0], RANKS[1] + 1):
            h, e = CE.exponents(kind, n)
            out.append({"group": f"{kind}{n}", "rank": n, "h": h, "exponents": e,
                        "n_exponents_equals_rank": len(e) == n})
    for n in (7, 8):
        h, e = CE.exponents("E", n)
        out.append({"group": f"E{n}", "rank": n, "h": h, "exponents": e,
                    "n_exponents_equals_rank": len(e) == n})
    return out


def run(seed=0, n_max=N_MAX):
    logged = emlib.load_run(emlib.ROOT / "runs/coxeter_exponents_seed0.json")["outputs"]["chains"]
    r4 = {(c["n_tetrahedra_to_close"], c["p"], c["q"]): c for c in emlib.load_run(
        emlib.ROOT / "runs/corkscrew_seed0.json")["outputs"]["tetrahedra_chain"]["r4"]["closed_chains"]}
    r6 = emlib.load_run(emlib.ROOT / "runs/corkscrew6_seed0.json")["outputs"]["chains"][
        "chains_using_all_three_planes"]
    rows = []
    for c in logged:
        m, n_close, turns = c["simplex_m"], c["n"], c["turns"]
        if c["dimension"] == 4:
            src = r4[(n_close, turns[0], turns[1])]
            shares = [src["r1_squared_share"], 1.0 - src["r1_squared_share"]]
        else:
            src = next(x for x in r6[str(m)]["listed"]
                       if (x["n"], x["p"], x["q"], x["r"]) == (n_close, *turns))
            shares = src["radius_squared_shares"]
        span = span_of(n_close, turns, shares)
        row = {"dimension_logged": c["dimension"], "simplex_m": m, "simplex_vertices_n": m + 1,
               "n_to_close": n_close, "turns": turns, "coxeter_verdict": c["verdict"],
               "survivor_unmatched": c["verdict"] == "no match",
               "edge_over_radius": src["edge_over_radius"],
               "points_span_dimension": span, "sphere_dimension_from_span": span - 1,
               "n_equals_span": m + 1 == span}
        if (c["dimension"], m) in ((4, 3), (6, 5)):
            row["tag"] = "SETTING_IN"
            row["eq_5_1"] = eq51(n_close, turns, m + 1, src["edge_over_radius"])
        elif (c["dimension"], m) == (6, 4):
            row["tag"] = "SETTING_OUTSIDE"
            row["eq_5_1"] = None
        else:
            row["tag"] = None          # 6-simplices in R^6: not assigned by the order
            row["eq_5_1"] = None
        row["tag_agrees_with_span"] = (None if row["tag"] is None else
                                       (row["tag"] == "SETTING_IN") == row["n_equals_span"])
        rows.append(row)
    slack = {str(d): search(*a, n_max=n_max) for d, a in
             ((5, (2, True)), (6, (3, False)), (7, (3, True)))}
    for c in slack["5"]["listed"]:
        c["tag"] = "SETTING_IN"
        c["eq_5_1"] = eq51(c["n"], c["turns"], 5, c["edge_over_radius"])
    logged_m4 = {(r["n_to_close"], tuple(r["turns"])) for r in rows
                 if r["dimension_logged"] == 6 and r["simplex_m"] == 4}
    found56 = {(c["n"], tuple(c["turns"])) for d in ("5", "6") for c in slack[d]["listed"]}
    n_sets = {d: set(slack[d]["n_values"]) for d in slack}
    sub = []
    for c in slack["7"]["listed"]:
        t = c["turns"]
        sub.append({"n": c["n"], "turns": t,
                    "holds_a_5d_chain_s_turns": [list(x["turns"]) for x in slack["5"]["listed"]
                                                 if x["n"] == c["n"] and set(x["turns"]) <= set(t)],
                    "holds_a_6d_chain_s_turns": [list(x["turns"]) for x in slack["6"]["listed"]
                                                 if x["n"] == c["n"] and set(x["turns"]) <= set(t)]})
    persist = {"n_values": {d: slack[d]["n_values"] for d in slack},
               "n_in_5d_and_6d": sorted(n_sets["5"] & n_sets["6"]),
               "n_in_5d_and_7d": sorted(n_sets["5"] & n_sets["7"]),
               "n_in_6d_and_7d": sorted(n_sets["6"] & n_sets["7"]),
               "n_in_all_three": sorted(n_sets["5"] & n_sets["6"] & n_sets["7"]),
               "turn_sets_5d": [[c["n"], c["turns"]] for c in slack["5"]["listed"]],
               "turn_sets_6d": [[c["n"], c["turns"]] for c in slack["6"]["listed"]],
               "seven_d_chains_holding_a_5d_turn_set":
                   [s for s in sub if s["holds_a_5d_chain_s_turns"]],
               "seven_d_chains_holding_a_6d_turn_set":
                   [s for s in sub if s["holds_a_6d_chain_s_turns"]],
               "lower_chains_met_again_with_a_zero_share": {
                   d: {"five_d_chains": [z for z in slack[d][
                           "zero_share_solutions_by_the_turns_they_live_on"]
                           if [z["n"], z["turns_with_a_share"]] in
                           [[c["n"], c["turns"]] for c in slack["5"]["listed"]]],
                       "six_d_chains": [z for z in slack[d][
                           "zero_share_solutions_by_the_turns_they_live_on"]
                           if [z["n"], z["turns_with_a_share"]] in
                           [[c["n"], c["turns"]] for c in slack["6"]["listed"]]],
                       "other": [z for z in slack[d][
                           "zero_share_solutions_by_the_turns_they_live_on"]
                           if [z["n"], z["turns_with_a_share"]] not in
                           [[c["n"], c["turns"]] for dd in ("5", "6")
                            for c in slack[dd]["listed"]]]}
                   for d in ("6", "7")},
               "five_and_six_d_together_equal_the_logged_6d_list": found56 == logged_m4,
               "logged_6d_chains_found_in_5d": sorted(
                   [n, list(t)] for n, t in logged_m4
                   if (n, t) in {(c["n"], tuple(c["turns"])) for c in slack["5"]["listed"]}),
               "note": "a plane turned by 180 degrees holds its points on a line, so the "
                       "chains logged in R^6 with a turn of N / 2 are the R^5 chains"}
    groups = coxeter_groups()
    surv = [r for r in rows if r["survivor_unmatched"]]
    proj, vsets = [], []
    for r in surv:
        src = next(x for x in r6["4"]["listed"]
                   if (x["n"], x["p"], x["q"], x["r"]) == (r["n_to_close"], *r["turns"]))
        proj.append(projection(r["n_to_close"], r["turns"], groups))
        vsets.append({"n": r["n_to_close"], "turns": r["turns"],
                      **vertex_set(r["n_to_close"], r["turns"], src["radius_squared_shares"])})
    d13 = next(g for g in groups if g["group"] == "D13")
    s_in = [r for r in rows if r["tag"] == "SETTING_IN"]
    return {"equation": "tan(n theta) = (1 + sec(2 pi / w)) tan(theta)   (R5 eq 5.1, carried)",
            "bound": {"seed": seed, "tol": TOL, "theta": "half of each plane's turn",
                      "w": "360 / the chain's own dihedral angle in degrees; not rounded",
                      "dihedral": "cos d = c / (1 + (n - 2) c), c = 1 - (edge / radius)^2 / 2",
                      "bounded_form": "cos d sin(n th) cos th - (1 + cos d) sin th cos(n th)",
                      "search_n_max": n_max, "simplex_m_searched": 4,
                      "search_not_done": "rotations with a fixed axis (turn 0): the chain "
                                         "then lies in a hyperplane off the centre",
                      "groups": f"B and D of rank {RANKS[0]} to {RANKS[1]}, E7, E8; A left "
                                "out as ordered; exponents computed from the diagrams",
                      "complements": "N / 2 - p (180 degrees minus the turn) and N - p",
                      "gram_rounding": 1e-9},
            "chains": rows, "slack": slack, "persistence": persist,
            "groups": groups, "projection": proj, "vertex_sets": vsets,
            "prediction": {"stated": "N = 24 (1, 7, 12) within the exponents of D13",
                           "d13_h": d13["h"], "d13_exponents": d13["exponents"],
                           "holds": all(p in d13["exponents"] for p in (1, 7, 12))
                           and d13["h"] == 24},
            "summary": {"setting_in": len(s_in),
                        "setting_in_holding_on_every_plane":
                            sum(1 for r in s_in if r["eq_5_1"]["holds_on_every_plane"]),
                        "setting_outside": sum(1 for r in rows if r["tag"] == "SETTING_OUTSIDE"),
                        "tag_void": sum(1 for r in rows if r["tag"] is None),
                        "tags_not_agreeing_with_span":
                            [[r["dimension_logged"], r["simplex_m"], r["n_to_close"], r["turns"]]
                             for r in rows if r["tag_agrees_with_span"] is False],
                        "five_d_chains": slack["5"]["n_found"],
                        "five_d_chains_holding_on_every_plane":
                            sum(1 for c in slack["5"]["listed"]
                                if c["eq_5_1"]["holds_on_every_plane"]),
                        "six_d_chains": slack["6"]["n_found"],
                        "seven_d_chains": slack["7"]["n_found"],
                        "survivors": len(surv),
                        "survivors_with_a_group": sum(1 for p in proj if p["verdict"] != "none")}}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    comp = emlib.compressibility(
        ["v" if r["tag"] is None else ("i" if r["tag"] == "SETTING_IN" else "o")
         for r in out["chains"]],
        "one letter per logged chain: i SETTING_IN, o SETTING_OUTSIDE, v tag void")
    angles = {f"{r['n_to_close']}:{r['turns']}": r["eq_5_1"]["dihedral_deg"]
              for r in out["chains"] if r["eq_5_1"]}
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME, seed, {"n_max": N_MAX, "ranks": list(RANKS), "tol": TOL,
                                     "status_before_run": "PROPOSED"},
                        out, angles, None, comp, None)
    print(p)
    for r in out["chains"]:
        e = r["eq_5_1"]
        print(f"{r['dimension_logged']}D m={r['simplex_m']} N={r['n_to_close']:2d} {str(r['turns']):12s} "
              f"span {r['points_span_dimension']} tag {r['tag']} agrees {r['tag_agrees_with_span']}"
              + (f" | dihedral {e['dihedral_deg']:.6f} w {e['w']:.6f} resid "
                 f"{e['largest_abs_residual']:.1e} holds {e['holds_on_every_plane']}" if e else ""))
    for d in ("5", "6", "7"):
        s = out["slack"][d]
        print(f"{d}D: found {s['n_found']} N {s['n_values']} tried {s['turn_sets_tried']} "
              f"singular {s['singular_systems_skipped']} neg {s['solutions_with_a_negative_share']} "
              f"zero {s['solutions_with_a_zero_share']}")
        for c in s["listed"][:60]:
            e = c.get("eq_5_1")
            print("   ", c["n"], c["turns"], [round(x, 4) for x in c["radius_squared_shares"]],
                  round(c["edge_over_radius"], 6), c["next_gap_also_equal"],
                  (f"dihedral {e['dihedral_deg']:.6f} w {e['w']:.6f} resid "
                   f"{e['largest_abs_residual']:.1e} holds {e['holds_on_every_plane']} "
                   f"non90 {e['planes_holding_not_at_theta_90']}") if e else "")
    pe = dict(out["persistence"])
    print("PERSIST", {k: pe[k] for k in pe if k.startswith("n_") or k.startswith("five")
                      or k.startswith("logged")})
    print("7D holding 5D:", pe["seven_d_chains_holding_a_5d_turn_set"])
    print("7D holding 6D:", pe["seven_d_chains_holding_a_6d_turn_set"])
    for q in out["projection"]:
        print("PROJ", q["n"], q["turns"], "->", q["verdict"],
              [(g["group"], g["turns_not_exponents"]) for g in q["groups_with_h_equal_n"]])
    print("PREDICTION", out["prediction"]["holds"])
    for v in out["vertex_sets"]:
        print("GRAM", v["n"], v["turns"], "distinct", v["n_distinct_inner_products"],
              "antipodal pairs", v["antipodal_pairs"], v["distinct_inner_products"])
    print("SUMMARY", out["summary"])


if __name__ == "__main__":
    main(sys.argv[1:])
