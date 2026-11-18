"""at. conjugate_map -- every logged closed chain's turn set multiplied by the
units mod N; the sets this reaches that are not logged chains, solved; the
shares of two survivors and the distances between their circles in exact
form.

CC0. Co-authored. Stdlib only. Status: PROPOSED.

Inputs (no chain is searched for here):
  runs/helix_setting_seed0.json.gz   the logged chains with their shares, and
                                     the m = 4 chains found in R^5 / R^6 / R^7
  runs/helix_isoclinic_seed0.json.gz the distances between sub-orbit circles
  corkscrew / corkscrew6 / coxeter_exponents runs, through helix_isoclinic

SHARES are the per-plane squared radius fractions as helix_setting defines
them: for turns p_i and shares t_i (sum 1) the squared distance at step gap j
is 2 * sum_i t_i (1 - cos(2 pi j p_i / N)); gaps 1 .. m agree.

  1  orbit map   turn set times u for every unit u mod N, each turn reduced
                 to min(t, N - t), sorted; the distinct sets are the orbit;
                 each member marked FOUND (a chain of the same simplex size
                 is logged on it) or solved.
  2  verdicts    a member that is not a logged chain: the share equations are
                 solved for it. REAL (every share > 0) / NEGATIVE_SHARE /
                 NO_SOLUTION. A share of exactly 0 is neither and is kept
                 apart as ZERO_SHARE. Exact in Q(zeta_N) for the logged
                 chains and the R^5 / R^6 ones; floats, with helix_setting's
                 own solver and tolerances, for the R^7 list.
  3  exact       18 (1,5,8) in the basis 1, c, c^2 with c = cos(2 pi / 9);
                 28 (1,8,13) and 28 (5,9,12) tried in 1, c, c^2 with
                 c = cos(2 pi / 7), then in that basis times {1, sqrt 7}.
                 The map zeta -> zeta^5 applied to the shares of 28 (1,8,13).
  4  circles     cos^2 of the distance between two sub-orbit circles of the
                 smallest isoclinic power, exact:
                 (sum t_i cos(2 pi k p_i / N))^2 + (sum e_i t_i sin(2 pi k p_i / N))^2
                 with e_i the power's sign on plane i and k the step between
                 the circles' starting vertices.

DERIVED, and checked here on every exactly solved member: the equations for
the set u P are the equations for P with zeta -> zeta^u applied to every
coefficient, so where the solution is unique the shares of u P are the
images of the shares of P. A conjugate set therefore has a solution whenever
the chain's own system is not singular; its shares need not be positive.

Nothing is interpreted. Observer slots null.

  python3 experiments/conjugate_map.py [seed]
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
import helix_isoclinic as HI  # noqa: E402
import helix_setting as HS  # noqa: E402
from exact.cyclo import Field, show, solve_exact  # noqa: E402

NAME = "conjugate_map"
TOL = 1e-9
NAMED_CHECKS = ((30, (1, 11), 7, (7, 13)), (28, (1, 8, 13), 5, (5, 9, 12)))
MISSING = ((28, (3, 4, 11)), (18, (4, 5, 7)), (18, (1, 2, 7)))
PREDICTION = ("PROPOSED, under test: the missing conjugate sets have a solution with a "
              "negative share; NO_SOLUTION on any of them would break it")


def reduce_set(n, turns, u):
    """(sorted reduced set, the reduced image of each turn in the given order)."""
    img = [min((u * p) % n, n - (u * p) % n) for p in turns]
    return tuple(sorted(img)), img


def orbit_of(n, turns):
    """{set: smallest unit that reaches it}, in order of first reach."""
    out = {}
    for u in range(1, n):
        if math.gcd(u, n) == 1:
            s, _ = reduce_set(n, turns, u)
            out.setdefault(s, u)
    return out


def _system(f, k, m, sub):
    last = f[-1]
    rows = [[(f[i][j] - last[j]) - (f[i][1] - last[1]) for i in range(k - 1)]
            for j in range(2, m + 1)]
    rhs = [last[1] - last[j] for j in range(2, m + 1)]
    return rows[:k - 1], rhs[:k - 1]


def float_shares(n, turns, m):
    """helix_setting.search's equations and tolerances, for one turn set."""
    f = [[1.0 - math.cos(2 * math.pi * j * p / n) for j in range(m + 1)] for p in turns]
    rows, rhs = _system(f, len(turns), m, None)
    sol = HS.solve(rows, rhs) if rows else []
    if sol is None:
        return {"verdict": "NO_SOLUTION", "why": "singular system", "shares": None}
    shares = sol + [1.0 - sum(sol)]
    d = [sum(s * fi[j] for s, fi in zip(shares, f)) for j in range(m + 1)]
    if any(abs(d[j] - d[1]) > 1e-9 * max(abs(d[1]), 1e-12) for j in range(2, m + 1)):
        return {"verdict": "NO_SOLUTION", "why": "the gaps do not agree", "shares": shares}
    lo = min(shares)
    v = "NEGATIVE_SHARE" if lo < -1e-9 else "ZERO_SHARE" if lo <= 1e-9 else "REAL"
    return {"verdict": v, "why": None, "shares": shares, "gap_squared_over_2": d[1]}


def exact_shares(n, turns, m, field=None):
    """The same equations in Q(zeta_N). Returns the field, the shares (or
    None), the half squared distance per gap 0 .. m + 1, and the verdict."""
    k_f = field or Field(n)
    f = [[k_f.one - k_f.cos(j * p) for j in range(m + 2)] for p in turns]
    rows, rhs = _system(f, len(turns), m, None)
    sol = solve_exact(rows, rhs, k_f.zero, k_f.one) if rows else []
    if sol is None:
        return {"field": k_f, "shares": None, "d": None, "verdict": "NO_SOLUTION",
                "why": "singular system"}
    shares = sol + [k_f.one - sum(sol, k_f.zero)]
    d = [sum((s * fi[j] for s, fi in zip(shares, f)), k_f.zero) for j in range(m + 2)]
    if any(d[j] != d[1] for j in range(2, m + 1)):
        return {"field": k_f, "shares": shares, "d": d, "verdict": "NO_SOLUTION",
                "why": "the gaps do not agree"}
    signs = [s.sign() for s in shares]
    v = "NEGATIVE_SHARE" if min(signs) < 0 else "ZERO_SHARE" if 0 in signs else "REAL"
    return {"field": k_f, "shares": shares, "d": d, "verdict": v, "why": None}


def sqrt7(k28):
    """sqrt 7 in Q(zeta_28), positive: -i times the quadratic Gauss sum of 7."""
    z = lambda e: k28.zeta(4 * e)                           # noqa: E731  zeta_7^e
    g = z(1) + z(2) + z(4) - z(3) - z(5) - z(6)             # sqrt(-7)
    r = -k28.zeta(7) * g
    return r if r.sign() > 0 else -r


def bases(n):
    """The bases tried for N: [(label, names, generators)], smallest first."""
    k_f = Field(n)
    if n == 18:
        c = k_f.cos(2)
        return k_f, [("Q(cos 2pi/9): 1, c, c^2", ["1", "c", "c^2"], [k_f.one, c, c * c])]
    if n == 28:
        c, r = k_f.cos(4), sqrt7(k_f)
        return k_f, [("Q(cos 2pi/7): 1, c, c^2", ["1", "c", "c^2"], [k_f.one, c, c * c]),
                     ("Q(cos 2pi/7, sqrt7): 1, c, c^2, r, r c, r c^2  (r = sqrt 7)",
                      ["1", "c", "c^2", "r", "r*c", "r*c^2"],
                      [k_f.one, c, c * c, r, r * c, r * c * c])]
    return k_f, []


def closed_form(k_f, x, tried):
    """x in the first basis that holds it: label, text, and whether each
    smaller basis failed."""
    out = {"value": float(x.approx()), "in": None, "form": None, "not_in": []}
    for label, names, gens in tried:
        q = k_f.coordinates(x, gens)
        if q is None:
            out["not_in"].append(label)
            continue
        out["in"], out["form"] = label, show(q, names)
        out["coefficients"] = [str(a) for a in q]
        break
    return out


def exact_report(n, turns, m):
    k_f, tried = bases(n)
    e = exact_shares(n, turns, m, k_f)
    rows = [dict(turn=p, **closed_form(k_f, s, tried)) for p, s in zip(turns, e["shares"])]
    return {"n": n, "turns": list(turns), "simplex_m": m, "verdict": e["verdict"],
            "gaps_1_to_m_exactly_equal": all(e["d"][j] == e["d"][1] for j in range(2, m + 1)),
            "gap_m_plus_1_equal_too": e["d"][m + 1] == e["d"][1],
            "half_squared_edge": closed_form(k_f, e["d"][1], tried),
            "shares": rows,
            "every_share_in_the_smallest_basis": all(not r["not_in"] for r in rows),
            "shares_sum_to_1_exactly": sum(e["shares"], k_f.zero) == k_f.one}


def galois_check(n, turns, m, u, field=None):
    """Shares of P under zeta -> zeta^u against the shares solved for u P."""
    k_f = field or Field(n)
    a = exact_shares(n, turns, m, k_f)
    target, img = reduce_set(n, turns, u)
    b = exact_shares(n, target, m, k_f)
    if a["shares"] is None or b["shares"] is None:
        return {"n": n, "turns": list(turns), "u": u, "image": list(target), "equal": None}
    where = {p: s for p, s in zip(target, b["shares"])}
    pairs = [(p, q, s.galois(u) == where[q]) for p, q, s in zip(turns, img, a["shares"])]
    return {"n": n, "turns": list(turns), "u": u, "image": list(target),
            "turn_goes_to": [[p, q] for p, q, _ in pairs],
            "equal": all(ok for _, _, ok in pairs)}


def circle_distances(n, turns, m, logged_deg):
    k_f, tried = bases(n)
    e = exact_shares(n, turns, m, k_f)
    power = HI.isoclinic_powers(n, turns)[0]
    s, signs = power["s"], power["signs"]
    classes = math.gcd(s, n)
    rows, seen = [], {}
    for k in range(1, classes):
        a = sum((t * k_f.cos(k * p) for t, p in zip(e["shares"], turns)), k_f.zero)
        b_i = sum((t * k_f.sin_i(k * p) * sg for t, p, sg in zip(e["shares"], turns, signs)),
                  k_f.zero)
        c2 = a * a - b_i * b_i                               # A^2 + B^2, B = b_i / i
        deg = math.degrees(math.acos(math.sqrt(max(0.0, min(1.0, float(c2.approx()))))))
        rows.append({"step_between_circles": k, "cos_squared": closed_form(k_f, c2, tried),
                     "distance_deg": deg})
        seen.setdefault(c2, []).append(k)
    distinct = sorted({round(r["distance_deg"], 6) for r in rows})
    return {"n": n, "turns": list(turns), "isoclinic_power_s": s, "signs": signs,
            "n_circles": classes, "polygon": f"regular {power['sub_orbit_size']}-gon",
            "by_step": rows,
            "distinct_cos_squared_values": len(seen),
            "distinct_distances_deg": distinct,
            "logged_distances_deg": logged_deg,
            "agrees_with_helix_isoclinic_at_1e-6": distinct == sorted(logged_deg),
            "every_cos_squared_in_the_smallest_basis":
                all(not r["cos_squared"]["not_in"] for r in rows)}


def gather():
    """Every logged chain as {m, n, turns, shares, source}, and per simplex
    size the set of (n, turns) that are logged chains."""
    hs = emlib.load_run(emlib.ROOT / "runs/helix_setting_seed0.json.gz")["outputs"]
    chains = []
    for c in HI.logged_chains():
        chains.append({"source": f"logged {c['dimension_logged']}D", "simplex_m": c["simplex_m"],
                       "n": c["n"], "turns": list(c["turns"]), "shares": c["shares"],
                       "survivor_unmatched": c["survivor_unmatched"]})
    seen = {(c["simplex_m"], c["n"], tuple(c["turns"])) for c in chains}
    for d in ("5", "6", "7"):
        for c in hs["slack"][d]["listed"]:
            key = (4, c["n"], tuple(c["turns"]))
            if key in seen:
                continue
            seen.add(key)
            chains.append({"source": f"helix_setting search R^{d}", "simplex_m": 4, "n": c["n"],
                           "turns": list(c["turns"]), "shares": c["radius_squared_shares"],
                           "survivor_unmatched": False})
    found = {}
    for c in chains:
        found.setdefault(c["simplex_m"], set()).add((c["n"], tuple(c["turns"])))
    return chains, found, hs


def orbit_row(c, found, exact):
    n, turns, m = c["n"], tuple(c["turns"]), c["simplex_m"]
    members = []
    k_f = Field(n) if exact else None
    for s, u in orbit_of(n, turns).items():
        is_found = (n, s) in found[m]
        row = {"turns": list(s), "first_unit": u, "found_chain": is_found}
        if is_found:
            row["verdict"] = "FOUND"
        elif exact:
            e = exact_shares(n, s, m, k_f)
            row.update(verdict=e["verdict"], why=e["why"], solver="exact, Q(zeta_N)",
                       shares=[float(x.approx()) for x in e["shares"]] if e["shares"] else None,
                       shares_are_the_galois_images=galois_check(n, turns, m, u, k_f)["equal"])
        else:
            e = float_shares(n, s, m)
            row.update(verdict=e["verdict"], why=e["why"],
                       solver="float, helix_setting's solver and tolerances",
                       shares=e["shares"])
        members.append(row)
    return {"source": c["source"], "simplex_m": m, "n": n, "turns": list(turns),
            "survivor_unmatched": c["survivor_unmatched"],
            "units": sum(1 for u in range(1, n) if math.gcd(u, n) == 1),
            "orbit_size": len(members),
            "members_found": sum(1 for r in members if r["found_chain"]),
            "orbit": members}


def run(seed=0):
    chains, found, hs = gather()
    iso = emlib.load_run(emlib.ROOT / "runs/helix_isoclinic_seed0.json.gz")["outputs"]
    logged_dist = {(r["n"], tuple(r["turns"])): r["smallest"]["hopf"]["distances_between_circles_deg"]
                   for r in iso["isoclinic_logged"] if r["smallest"]
                   and r["smallest"]["hopf"].get("distances_between_circles_deg")}
    rows = [orbit_row(c, found, exact=not c["source"].endswith("R^7")) for c in chains]
    checks = []
    for n, turns, u, want in NAMED_CHECKS:
        got, _ = reduce_set(n, turns, u)
        checks.append({"n": n, "turns": list(turns), "u": u, "gives": list(got),
                       "expected": list(want), "holds": got == tuple(want)})
    missing = []
    for n, turns in MISSING:
        e = exact_shares(n, turns, 4)
        f = float_shares(n, turns, 4)
        rep = exact_report(n, turns, 4)
        parent = next(((r["n"], r["turns"]) for r in rows
                       if r["n"] == n and r["survivor_unmatched"]
                       and any(x["turns"] == list(turns) for x in r["orbit"])), None)
        missing.append({"n": n, "turns": list(turns), "conjugate_of": parent,
                        "verdict": e["verdict"], "why": e["why"],
                        "float_solver_verdict": f["verdict"],
                        "shares": [float(x.approx()) for x in e["shares"]] if e["shares"] else None,
                        "signs": [x.sign() for x in e["shares"]] if e["shares"] else None,
                        "negative_on_turns": [p for p, x in zip(turns, e["shares"])
                                              if x.sign() < 0] if e["shares"] else None,
                        "exact": rep})
    exact = {"18_1_5_8": exact_report(18, (1, 5, 8), 4),
             "28_1_8_13": exact_report(28, (1, 8, 13), 4),
             "28_5_9_12": exact_report(28, (5, 9, 12), 4)}
    k28 = Field(28)
    r7 = sqrt7(k28)
    c7 = k28.cos(4)
    galois = {
        "28_1_8_13_under_5": galois_check(28, (1, 8, 13), 4, 5, k28),
        "on_the_generators": {
            "c = cos(2pi/7) goes to cos(2pi 5/7)": c7.galois(5) == k28.cos(20),
            "cos(2pi 5/7) = cos(2pi 2/7)": k28.cos(20) == k28.cos(8),
            "sqrt7 goes to -sqrt7": r7.galois(5) == -r7,
            "sqrt7 squared is 7": r7 * r7 == 7},
        "unit_13_on_28_1_8_13": {
            "image": list(reduce_set(28, (1, 8, 13), 13)[0]),
            "turn_goes_to": [[p, q] for p, q in zip((1, 8, 13), reduce_set(28, (1, 8, 13), 13)[1])],
            "note": "13 fixes the set and swaps turns 1 and 13; it fixes cos(2pi/7) and "
                    "sends sqrt7 to -sqrt7",
            "fixes_c": c7.galois(13) == c7, "sqrt7_goes_to_minus": r7.galois(13) == -r7},
        "18_1_5_8_under_5": galois_check(18, (1, 5, 8), 4, 5),
        "18_1_5_8_under_7": galois_check(18, (1, 5, 8), 4, 7),
        "28_1_8_13_under_3": galois_check(28, (1, 8, 13), 4, 3, k28),
        "30_1_11_under_7": galois_check(30, (1, 11), 3, 7)}
    circles = {"18_1_5_8": circle_distances(18, (1, 5, 8), 4, logged_dist[(18, (1, 5, 8))]),
               "28_1_8_13": circle_distances(28, (1, 8, 13), 4, logged_dist[(28, (1, 8, 13))]),
               "28_5_9_12": circle_distances(28, (5, 9, 12), 4, logged_dist[(28, (5, 9, 12))])}
    tally = {}
    for r in rows:
        for x in r["orbit"]:
            if not x["found_chain"]:
                key = ("R^7 list" if r["source"].endswith("R^7") else "logged and R^5 / R^6")
                tally.setdefault(key, {}).setdefault(x["verdict"], set()).add(
                    (r["simplex_m"], r["n"], tuple(x["turns"])))
    exact_members = [x for r in rows for x in r["orbit"]
                     if x.get("shares_are_the_galois_images") is not None]
    not_r7 = [r for r in rows if not r["source"].endswith("R^7")]
    r7_rows = [r for r in rows if r["source"].endswith("R^7")]
    summary = {
        "chains_read": len(rows), "of_them_from_the_r7_list": len(r7_rows),
        "named_checks_hold": all(c["holds"] for c in checks),
        "orbits_with_every_member_found": sum(1 for r in rows
                                              if r["members_found"] == r["orbit_size"]),
        "orbits_with_a_member_not_found": sum(1 for r in rows
                                              if r["members_found"] < r["orbit_size"]),
        "not_found_members_by_verdict": {k: {v: len(s) for v, s in sorted(d.items())}
                                         for k, d in sorted(tally.items())},
        "logged_and_r5_r6_orbits_with_a_member_not_found": [
            [r["simplex_m"], r["n"], r["turns"],
             [[x["turns"], x["verdict"]] for x in r["orbit"] if not x["found_chain"]]]
            for r in not_r7 if r["members_found"] < r["orbit_size"]],
        "missing_three_verdicts": [[m["n"], m["turns"], m["verdict"]] for m in missing],
        "prediction": PREDICTION,
        "prediction_holds_on_the_three": all(m["verdict"] == "NEGATIVE_SHARE" for m in missing),
        "any_no_solution_anywhere": any(x["verdict"] == "NO_SOLUTION"
                                        for r in rows for x in r["orbit"]),
        "any_real_member_that_is_not_a_logged_chain": [
            [r["simplex_m"], r["n"], x["turns"]] for r in rows for x in r["orbit"]
            if x["verdict"] == "REAL"],
        "galois_image_equals_solved_shares": {
            "members_checked": len(exact_members),
            "equal": sum(1 for x in exact_members if x["shares_are_the_galois_images"])},
        "shares_18_all_in_Q_cos_2pi_9": exact["18_1_5_8"]["every_share_in_the_smallest_basis"],
        "shares_28_1_8_13_all_in_Q_cos_2pi_7":
            exact["28_1_8_13"]["every_share_in_the_smallest_basis"],
        "shares_28_1_8_13_in_Q_cos_2pi_7_by_turn": [
            [s["turn"], not s["not_in"]] for s in exact["28_1_8_13"]["shares"]],
        "shares_28_5_9_12_in_Q_cos_2pi_7_by_turn": [
            [s["turn"], not s["not_in"]] for s in exact["28_5_9_12"]["shares"]],
        "circle_distances_agree_with_helix_isoclinic":
            all(c["agrees_with_helix_isoclinic_at_1e-6"] for c in circles.values()),
        "circle_cos_squared_in_the_smallest_basis": {
            k: c["every_cos_squared_in_the_smallest_basis"] for k, c in circles.items()}}
    return {"bound": {"seed": seed,
                      "chains": "the logged chains (corkscrew, corkscrew6) and the m = 4 chains "
                                "helix_setting found in R^5, R^6 and R^7, N <= "
                                f"{HS.N_MAX}; none searched for here",
                      "found_means": "a logged chain of the same simplex size on that turn set",
                      "exact": "Q(zeta_N), Fractions; zero, equality and the sign of a share "
                               "all decided exactly (rational intervals for the sign); "
                               "floats only where a value is printed",
                      "float": "the R^7 list only: helix_setting's solver; negative below "
                               "-1e-9, zero within 1e-9, gaps equal within 1e-9 relative",
                      "bases_tried": {"18": ["Q(cos 2pi/9)"],
                                      "28": ["Q(cos 2pi/7)", "Q(cos 2pi/7, sqrt7)"]},
                      "circle_distances": "the smallest isoclinic power only, as in "
                                          "helix_isoclinic; compared at 1e-6 degrees"},
            "named_checks": checks, "orbits": rows, "missing": missing, "exact": exact,
            "galois": galois, "circles": circles, "summary": summary}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    letter = {"FOUND": "f", "NEGATIVE_SHARE": "n", "REAL": "r", "ZERO_SHARE": "z",
              "NO_SOLUTION": "x"}
    comp = emlib.compressibility(
        [letter[x["verdict"]] for r in out["orbits"] for x in r["orbit"]],
        "one letter per orbit member, chains in input order: f found chain, n negative "
        "share, r real and not logged, z zero share, x no solution")
    angles = {k: c["distinct_distances_deg"] for k, c in out["circles"].items()}
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME, seed, {"sources": ["runs/helix_setting_seed0.json.gz",
                                                 "runs/helix_isoclinic_seed0.json.gz",
                                                 "runs/corkscrew_seed0.json",
                                                 "runs/corkscrew6_seed0.json",
                                                 "runs/coxeter_exponents_seed0.json"],
                                     "tol": TOL, "status_before_run": "PROPOSED"},
                        out, angles, None, comp, None)
    print(p)
    for c in out["named_checks"]:
        print(f"CHECK {c['u']} x {tuple(c['turns'])} mod {c['n']} = {tuple(c['gives'])} "
              f"expected {tuple(c['expected'])} holds {c['holds']}")
    for r in out["orbits"]:
        if r["source"].endswith("R^7"):
            continue
        print(f"ORBIT m={r['simplex_m']} N={r['n']:2d} {str(tuple(r['turns'])):12s} "
              f"{'*' if r['survivor_unmatched'] else ' '} units {r['units']:2d} "
              f"orbit {r['orbit_size']} found {r['members_found']} | "
              + "; ".join(f"{tuple(x['turns'])} {x['verdict']}" for x in r["orbit"]))
    for m in out["missing"]:
        print(f"MISSING N={m['n']} {tuple(m['turns'])} of {m['conjugate_of']}: {m['verdict']} "
              f"(float solver {m['float_solver_verdict']}) shares "
              f"{[round(s, 6) for s in m['shares']]} negative on {m['negative_on_turns']}")
        for s in m["exact"]["shares"]:
            print(f"     turn {s['turn']:2d}  {s['value']: .9f} = {s['form']}   [{s['in']}]")
    for k, e in out["exact"].items():
        print(f"EXACT {k}: gaps equal {e['gaps_1_to_m_exactly_equal']} sum 1 "
              f"{e['shares_sum_to_1_exactly']} all in smallest basis "
              f"{e['every_share_in_the_smallest_basis']}")
        for s in e["shares"]:
            print(f"     turn {s['turn']:2d}  {s['value']: .9f} = {s['form']}   [{s['in']}]"
                  + (f"  not in {s['not_in']}" if s["not_in"] else ""))
        h = e["half_squared_edge"]
        print(f"     half squared edge {h['value']:.9f} = {h['form']}   [{h['in']}]")
    g = out["galois"]
    for k in ("28_1_8_13_under_5", "28_1_8_13_under_3", "18_1_5_8_under_5", "18_1_5_8_under_7",
              "30_1_11_under_7"):
        print(f"GALOIS {k}: image {g[k]['image']} turns {g[k]['turn_goes_to']} "
              f"shares equal {g[k]['equal']}")
    print("GALOIS generators", g["on_the_generators"])
    print("GALOIS unit 13", g["unit_13_on_28_1_8_13"])
    for k, c in out["circles"].items():
        print(f"CIRCLES {k}: s={c['isoclinic_power_s']} signs {c['signs']} {c['n_circles']} x "
              f"{c['polygon']}; distances {c['distinct_distances_deg']} logged "
              f"{c['logged_distances_deg']} agree {c['agrees_with_helix_isoclinic_at_1e-6']}")
        shown = set()
        for r in c["by_step"]:
            f = r["cos_squared"]["form"]
            if f not in shown:
                shown.add(f)
                print(f"     step {r['step_between_circles']}: cos^2 = {f} = "
                      f"{r['cos_squared']['value']:.12f}  [{r['cos_squared']['in']}]  "
                      f"{r['distance_deg']:.6f} deg")
    print("SUMMARY", out["summary"])


if __name__ == "__main__":
    main(sys.argv[1:])
