"""as. helix_isoclinic -- setting tags from the span, one chain in exact
numbers, the (7, 13) chain beside the 600-cell, and which powers of each
chain's step turn every plane by one angle.

CC0. Co-authored. Stdlib only. Status: PROPOSED. Follows helix_setting.

  1  tags     SETTING_IN when the simplex's vertex count n equals the
              dimension its chain's points span (sphere dim = n - 1);
              SETTING_OUTSIDE otherwise, with which side (sphere dim above
              or below n - 1). The tag written by helix_setting on
              2026-11-16 is kept beside the new one as a record; that run
              file is not touched. Eq 5.1 residual for every SETTING_IN row.
  2  exact    the chain N = 24, turns (1, 7, 12), in a + b sqrt2 + c sqrt3
              + d sqrt6 with rational a..d: shares, pairwise cosine,
              cos(dihedral) against (sqrt2 - 1) / 2, 1 + sec against
              3 + 2 sqrt2, and eq 5.1 as an exact equality:
              tan(5 th) / tan(th) = (sin 3t + sin 2t) / (sin 3t - sin 2t)
              with t = 2 th the plane's turn (product-to-sum). The float
              residual is read at 1e-12.
  3  (7, 13)  the chain N = 30, turns (7, 13): its inner products by step
              gap beside those of the 600-cell's Coxeter-element orbits
              taken 7 steps at a time (exact, from e8_h4_fold); which of the
              600-cell's chords its edge is; the Petrie readings.
  4  isoclinic  for a chain with step turns p_i / N: every s in 1 .. N - 1
              with s p_i = +c or -c (mod N) on all planes, c not 0. A line
              turned by 180 degrees counts as a turn of N / 2. For the
              smallest s the sub-orbits v_k, v_(k+s), ... are read: size,
              span, equal chords.
  5  Hopf     the circles of those sub-orbits, pair by pair: the two
              principal angles between their planes; equal = the circles
              keep a constant distance.
  6  all      item 4 on the 17 logged chains and the R^7 chains.

Nothing is interpreted. Observer slots null.

  python3 experiments/helix_isoclinic.py [seed]
"""
import itertools
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import e8_h4_fold as E  # noqa: E402
import emlib  # noqa: E402
import helix_setting as H  # noqa: E402
import spherical_helix as SH  # noqa: E402
from toolkit import primitives as P  # noqa: E402

NAME = "helix_isoclinic"
TOL = 1e-12
F = Fraction
S2, S3 = math.sqrt(2.0), math.sqrt(3.0)


class Q23:
    """a + b sqrt2 + c sqrt3 + d sqrt6, rational a..d. Exact."""
    __slots__ = ("v",)

    def __init__(self, a=0, b=0, c=0, d=0):
        self.v = (F(a), F(b), F(c), F(d))

    def __add__(self, o):
        o = o if isinstance(o, Q23) else Q23(o)
        return Q23(*[x + y for x, y in zip(self.v, o.v)])
    __radd__ = __add__

    def __neg__(self):
        return Q23(*[-x for x in self.v])

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Q23) else Q23(o)))

    def __rsub__(self, o):
        return Q23(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, Q23) else Q23(o)
        a, b, c, d = self.v
        e, f, g, h = o.v
        return Q23(a * e + 2 * b * f + 3 * c * g + 6 * d * h,
                   a * f + b * e + 3 * (c * h + d * g),
                   a * g + c * e + 2 * (b * h + d * f),
                   a * h + d * e + b * g + c * f)
    __rmul__ = __mul__

    def inverse(self):
        a, b, c, d = self.v
        conj2 = Q23(a, -b, c, -d)                 # sqrt2 -> -sqrt2
        n = self * conj2                          # in Q(sqrt3)
        p, _, q, _ = n.v
        den = p * p - 3 * q * q
        return conj2 * Q23(p / den, 0, -q / den, 0)

    def __truediv__(self, o):
        return self * (o if isinstance(o, Q23) else Q23(o)).inverse()

    def __eq__(self, o):
        return self.v == (o if isinstance(o, Q23) else Q23(o)).v

    def __hash__(self):
        return hash(self.v)

    def __float__(self):
        a, b, c, d = self.v
        return float(a) + float(b) * S2 + float(c) * S3 + float(d) * S2 * S3

    def __str__(self):
        out = [f"{x}{s}" for x, s in zip(self.v, ("", "*sqrt2", "*sqrt3", "*sqrt6")) if x]
        return " + ".join(out) if out else "0"


def cos15(k):
    """cos(15 k degrees), exact, by c_(k+1) = 2 c_1 c_k - c_(k-1)."""
    c1 = Q23(0, F(1, 4), 0, F(1, 4))
    a, b = Q23(1), c1
    for _ in range(abs(k) % 24):
        a, b = b, 2 * c1 * b - a
    return a


def sin15(k):
    return cos15(6 - k)


def exact_24():
    n, turns, m = 24, (1, 7, 12), 4
    f = [[Q23(1) - cos15(j * p) for j in range(m + 2)] for p in turns]
    fa, fb, fc = f
    a11 = (fa[2] - fc[2]) - (fa[1] - fc[1])
    a12 = (fb[2] - fc[2]) - (fb[1] - fc[1])
    a21 = (fa[3] - fc[3]) - (fa[1] - fc[1])
    a22 = (fb[3] - fc[3]) - (fb[1] - fc[1])
    b1, b2 = fc[1] - fc[2], fc[1] - fc[3]
    det = a11 * a22 - a12 * a21
    t1 = (b1 * a22 - a12 * b2) / det
    t2 = (a11 * b2 - b1 * a21) / det
    t3 = Q23(1) - t1 - t2
    d = [t1 * fa[j] + t2 * fb[j] + t3 * fc[j] for j in range(m + 2)]
    c = Q23(1) - d[1]
    cos_d = c / (Q23(1) + 3 * c)
    one_plus_sec = Q23(1) + cos_d.inverse()
    planes = []
    for p in (1, 7):
        lhs = sin15(3 * p) + sin15(2 * p)
        rhs = one_plus_sec * (sin15(3 * p) - sin15(2 * p))
        th = 180.0 * p / n
        planes.append({"turn": p, "turn_deg": 15.0 * p, "theta_deg": th,
                       "sin_3t_plus_sin_2t": str(lhs),
                       "one_plus_sec_times_sin_3t_minus_sin_2t": str(rhs),
                       "exactly_equal": lhs == rhs,
                       "float_residual_bounded_form": H.residual(th, 5, float(cos_d)),
                       "float_holds_at_1e-12": abs(H.residual(th, 5, float(cos_d))) < TOL})
    return {"chain": "N = 24, turns (1, 7, 12); 4-simplices; spans 5 dimensions",
            "field": "a + b sqrt2 + c sqrt3 + d sqrt6, rational a..d",
            "radius_squared_shares_exact": [str(t1), str(t2), str(t3)],
            "radius_squared_shares": [float(t1), float(t2), float(t3)],
            "gaps_1_to_4_exactly_equal": d[1] == d[2] == d[3] == d[4],
            "gap_5_equal_too": d[5] == d[1],
            "pairwise_cosine_exact": str(c), "pairwise_cosine": float(c),
            "edge_over_radius_squared_exact": str(2 * d[1]),
            "cos_dihedral_exact": str(cos_d), "cos_dihedral": float(cos_d),
            "cos_dihedral_equals_sqrt2_minus_1_over_2": cos_d == Q23(F(-1, 2), F(1, 2)),
            "one_plus_sec_exact": str(one_plus_sec),
            "one_plus_sec_equals_3_plus_2_sqrt2": one_plus_sec == Q23(3, 2),
            "dihedral_deg": math.degrees(math.acos(float(cos_d))),
            "w": 360.0 / math.degrees(math.acos(float(cos_d))),
            "eq_5_1_planes": planes,
            "line_plane": {"theta_deg": 90.0, "note": "cos(theta) = 0 and cos(5 theta) = 0: "
                           "both tangents unbounded; the bounded form is 0",
                           "float_residual_bounded_form": H.residual(90.0, 5, float(cos_d))},
            "eq_5_1_exact_on_both_planes": all(x["exactly_equal"] for x in planes)}


PHI_F = (1 + 5 ** 0.5) / 2


def zf(v):
    return float(v.a) + float(v.b) * PHI_F


def chain_7_13(share_7):
    n = 30
    cg = [share_7 * math.cos(2 * math.pi * 7 * g / n)
          + (1 - share_7) * math.cos(2 * math.pi * 13 * g / n) for g in range(n)]
    gram = E.e8_gram()
    rs = E.roots(gram)
    fold = {r: E.folded(r) for r in rs}
    back = {v: k for k, v in fold.items()}
    left, orbits = set(rs), []
    while left:
        r0 = min(left)
        orb, x = [], fold[r0]
        for _ in range(30):
            orb.append(x)
            x = E.coxeter_h4(x)
        orbits.append(orb)
        left -= {back[v] for v in orb}
    rows = []
    for i, orb in enumerate(orbits):
        n2 = E.h4_form(orb[0], orb[0])
        c = [zf(E.h4_form(orb[0], orb[g]) / n2) for g in range(30)]
        by_k = {}
        for k in (1, 7, 11, 13):          # the steps coprime with 30, up to direction
            by_k[k] = max(abs(a - c[(k * g) % 30]) for g, a in enumerate(cg))
        best = min(by_k, key=by_k.get)
        rows.append({"orbit": i, "squared_length": str(n2),
                     "on_the_unit_600_cell": str(n2) == "2",
                     "step_deg": math.degrees(math.acos(max(-1.0, min(1.0, c[1])))),
                     "largest_difference_taken_k_at_a_time": {str(k): v for k, v in by_k.items()},
                     "matches_taken_k_at_a_time": [k for k, v in by_k.items() if v < 1e-9],
                     "matches_taken_7_at_a_time": by_k[7] < 1e-9,
                     "closest_k": best})
    chords = sorted({round(math.degrees(math.acos(max(-1.0, min(1.0, zf(
        E.h4_form(orbits[k][0], y) / E.h4_form(orbits[k][0], orbits[k][0])))))), 6)
        for k in range(len(orbits)) if str(E.h4_form(orbits[k][0], orbits[k][0])) == "2"
        for o in orbits if str(E.h4_form(o[0], o[0])) == "2" for y in o} - {0.0})
    edge_deg = math.degrees(math.acos(cg[1]))
    win = 2
    while win < 8 and abs(cg[win] - cg[1]) < 1e-9:
        win += 1
    d = H.dihedral(cg[1], 4)
    match = [r["orbit"] for r in rows if r["matches_taken_k_at_a_time"]]
    return {"chain": "N = 30, turns (7, 13); tetrahedra in R^4", "share_on_turn_7": share_7,
            "cosine_by_step_gap": cg[1:16],
            "orbits_of_the_600_cells": rows,
            "orbits_matching": match,
            "vertex_set_is_30_vertices_of_a_600_cell":
                any(r["matches_taken_k_at_a_time"] and r["on_the_unit_600_cell"] for r in rows),
            "matches": [{"orbit": r["orbit"], "squared_length": r["squared_length"],
                         "orbit_step_deg": r["step_deg"], "k": r["matches_taken_k_at_a_time"]}
                        for r in rows if r["matches_taken_k_at_a_time"]],
            "matching_orbits_step_deg": [rows[i]["step_deg"] for i in match],
            "how": "equal inner products at every step gap = congruent point sets; the chain's "
                   "k-th vertex is the orbit's (7 k)-th",
            "n_vertices_of_the_600_cell": 120, "n_vertices_of_the_chain": 30,
            "chords_of_the_600_cell_deg": chords,
            "edge_angle_deg": edge_deg,
            "edge_is_chord_number": (chords.index(round(edge_deg, 6)) + 1
                                     if round(edge_deg, 6) in chords else None),
            "edge_over_radius": math.sqrt(2 - 2 * cg[1]),
            "edge_over_radius_equals_phi": abs(math.sqrt(2 - 2 * cg[1]) - PHI_F) < 1e-9,
            "edge_of_the_600_cell_deg": chords[0],
            "consecutive_vertices_forming_a_regular_simplex": win,
            "four_consecutive_are_a_regular_tetrahedron": win >= 4,
            "five_consecutive_are_a_regular_simplex": win >= 5,
            "dihedral_deg": d["deg"], "w": d["w"],
            "name_carried_not_verified": "{3,3,5/2}: the name is carried; what is checked "
                                         "here is the list of readings above",
            "petrie_polygon_of_3_3_5_over_2": None}


def isoclinic_powers(n, turns):
    out = []
    for s in range(1, n):
        r = [(s * p) % n for p in turns]
        c = min(r[0], n - r[0])
        if c == 0:
            continue
        if all(x == c or x == n - c for x in r):
            out.append({"s": s, "c": c, "angle_deg": 360.0 * c / n,
                        "signs": [1 if x == c else -1 for x in r] if 2 * c != n
                        else [0] * len(r),
                        "sub_orbit_size": n // math.gcd(c, n),
                        "n_sub_orbits": math.gcd(c, n),
                        "is_minus_identity": 2 * c == n})
    return out


def plane_basis(pts):
    return P.gram_schmidt([list(p) for p in pts], tol=1e-9)


def read_sub_orbits(n, turns, shares, power):
    pts = SH.orbit(n, turns, shares)
    s, size, count = power["s"], power["sub_orbit_size"], power["n_sub_orbits"]
    circles, ok_span, ok_chord, worst = [], True, True, 0.0
    for start in range(count):
        # the sub-orbits are the classes of k mod gcd(s, N) -- take one start in each
        sub = [pts[(start + j * s) % n] for j in range(size)]
        basis = plane_basis(sub)
        chords = [math.dist(sub[j], sub[(j + 1) % size]) for j in range(size)]
        spread = max(chords) - min(chords)
        worst = max(worst, spread)
        ok_span = ok_span and len(basis) == (2 if size > 2 else 1)
        ok_chord = ok_chord and spread < 1e-9
        circles.append(basis)
    starts = sorted({k % math.gcd(s, n) for k in range(n)})
    out = {"s": s, "c": power["c"], "angle_deg": power["angle_deg"],
           "sub_orbit_size": size, "n_sub_orbits": count,
           "sub_orbits_cover_every_vertex_once": len(starts) == count and count * size == n,
           "every_sub_orbit_spans_2_dimensions": ok_span if size > 2 else None,
           "every_sub_orbit_has_equal_chords": ok_chord if size > 2 else None,
           "largest_chord_spread": worst,
           "regular_polygons_on_great_circles": (ok_span and ok_chord) if size > 2 else None,
           "polygon": f"regular {size}-gon" if size > 2 else "two antipodal points"}
    if size > 2 and ok_span:
        pairs, par, worst_gap = [], True, 0.0
        for i, j in itertools.combinations(range(count), 2):
            a, b = circles[i], circles[j]
            mm = [[P.dot(a[r], b[c]) for c in range(2)] for r in range(2)]
            g11 = mm[0][0] ** 2 + mm[1][0] ** 2
            g22 = mm[0][1] ** 2 + mm[1][1] ** 2
            g12 = mm[0][0] * mm[0][1] + mm[1][0] * mm[1][1]
            tr, dt = g11 + g22, g11 * g22 - g12 * g12
            disc = math.sqrt(max(0.0, tr * tr / 4 - dt))
            s1 = math.sqrt(max(0.0, min(1.0, tr / 2 + disc)))
            s2 = math.sqrt(max(0.0, min(1.0, tr / 2 - disc)))
            # equal principal angles <=> M^T M is a multiple of the identity
            gap = max(abs(g11 - g22), abs(g12))
            if gap < 1e-9:
                s1 = s2 = math.sqrt(max(0.0, min(1.0, tr / 2)))
            a1, a2 = math.degrees(math.acos(s1)), math.degrees(math.acos(s2))
            worst_gap = max(worst_gap, gap)
            par = par and gap < 1e-9
            pairs.append({"circles": [i, j], "principal_angles_deg": [a1, a2],
                          "equal": gap < 1e-9})
        out["hopf"] = {"n_circles": count, "pairs_checked": len(pairs),
                       "pairwise_clifford_parallel": par if pairs else None,
                       "largest_departure_of_MtM_from_a_multiple_of_I": worst_gap,
                       "pairs_with_equal_angles": sum(1 for p in pairs if p["equal"]),
                       "distances_between_circles_deg": sorted(
                           {round(p["principal_angles_deg"][0], 6) for p in pairs
                            if p["equal"]}),
                       "principal_angles_of_unequal_pairs_deg": sorted(
                           {(round(p["principal_angles_deg"][0], 6),
                             round(p["principal_angles_deg"][1], 6)) for p in pairs
                            if not p["equal"]}),
                       "pairs": pairs}
    else:
        out["hopf"] = {"n_circles": None, "pairs_checked": 0,
                       "pairwise_clifford_parallel": None,
                       "note": "the sub-orbits are pairs of opposite points, not circles"
                       if size == 2 else "sub-orbits do not span planes"}
    return out


def logged_chains():
    logged = emlib.load_run(emlib.ROOT / "runs/coxeter_exponents_seed0.json")["outputs"]["chains"]
    r4 = {(c["n_tetrahedra_to_close"], c["p"], c["q"]): c for c in emlib.load_run(
        emlib.ROOT / "runs/corkscrew_seed0.json")["outputs"]["tetrahedra_chain"]["r4"]["closed_chains"]}
    r6 = emlib.load_run(emlib.ROOT / "runs/corkscrew6_seed0.json")["outputs"]["chains"][
        "chains_using_all_three_planes"]
    out = []
    for c in logged:
        m, n, turns = c["simplex_m"], c["n"], c["turns"]
        if c["dimension"] == 4:
            src = r4[(n, turns[0], turns[1])]
            shares = [src["r1_squared_share"], 1.0 - src["r1_squared_share"]]
        else:
            src = next(x for x in r6[str(m)]["listed"]
                       if (x["n"], x["p"], x["q"], x["r"]) == (n, *turns))
            shares = src["radius_squared_shares"]
        out.append({"dimension_logged": c["dimension"], "simplex_m": m, "n": n, "turns": turns,
                    "shares": shares, "edge_over_radius": src["edge_over_radius"],
                    "survivor_unmatched": c["verdict"] == "no match"})
    return out


def isoclinic_row(n, turns, shares, with_geometry=True):
    powers = isoclinic_powers(n, turns)
    row = {"n": n, "turns": list(turns), "has_an_isoclinic_power": bool(powers),
           "smallest_s": powers[0]["s"] if powers else None,
           "all_isoclinic_powers": [{k: p[k] for k in ("s", "c", "sub_orbit_size",
                                                        "n_sub_orbits", "is_minus_identity")}
                                    for p in powers],
           "only_minus_identity": bool(powers) and all(p["is_minus_identity"] for p in powers)}
    row["smallest"] = read_sub_orbits(n, turns, shares, powers[0]) \
        if powers and with_geometry else None
    return row


def run(seed=0):
    old = {(r["dimension_logged"], r["simplex_m"], r["n_to_close"], tuple(r["turns"])): r
           for r in emlib.load_run(emlib.ROOT / "runs/helix_setting_seed0.json.gz")["outputs"]["chains"]}
    chains = logged_chains()
    tags, iso = [], []
    for c in chains:
        n_v = c["simplex_m"] + 1
        span = SH.rank(SH.orbit(c["n"], c["turns"], c["shares"]))
        side = ("sphere dim = n - 1" if span == n_v else
                "sphere dim > n - 1" if span > n_v else "sphere dim < n - 1")
        o = old[(c["dimension_logged"], c["simplex_m"], c["n"], tuple(c["turns"]))]
        tag = "SETTING_IN" if span == n_v else "SETTING_OUTSIDE"
        e = H.eq51(c["n"], c["turns"], n_v, c["edge_over_radius"]) if tag == "SETTING_IN" else None
        if e:
            for x in e["planes"]:
                x["holds_at_1e-12"] = abs(x["residual_bounded_form"]) < TOL
            e["holds_on_every_plane_at_1e-12"] = all(x["holds_at_1e-12"] for x in e["planes"])
        tags.append({"dimension_logged": c["dimension_logged"], "simplex_m": c["simplex_m"],
                     "simplex_vertices_n": n_v, "n_to_close": c["n"], "turns": c["turns"],
                     "survivor_unmatched": c["survivor_unmatched"],
                     "points_span_dimension": span, "sphere_dimension": span - 1,
                     "relation": side, "tag": tag,
                     "tag_2026_11_16_as_ordered": o["tag"],
                     "tag_changed": o["tag"] != tag,
                     "eq_5_1": e,
                     "eq_5_1_2026_11_16_kept_as_record":
                         None if tag == "SETTING_IN" or not o["eq_5_1"] else
                         {"w": o["eq_5_1"]["w"], "dihedral_deg": o["eq_5_1"]["dihedral_deg"],
                          "largest_abs_residual": o["eq_5_1"]["largest_abs_residual"],
                          "read_with_n": n_v}})
        iso.append({"dimension_logged": c["dimension_logged"], "simplex_m": c["simplex_m"],
                    "survivor_unmatched": c["survivor_unmatched"],
                    **isoclinic_row(c["n"], c["turns"], c["shares"])})
    r7 = emlib.load_run(emlib.ROOT / "runs/helix_setting_seed0.json.gz")["outputs"]["slack"]["7"]
    iso7 = [isoclinic_row(c["n"], c["turns"], c["radius_squared_shares"]) for c in r7["listed"]]
    share7 = next(c["shares"][0] for c in chains if c["n"] == 30 and c["turns"] == [7, 13])
    with7 = [r for r in iso7 if r["has_an_isoclinic_power"]]
    return {"bound": {"seed": seed, "tol": TOL, "span_tolerance": 1e-9,
                      "tags": "SETTING_IN when simplex vertex count = dimension spanned",
                      "isoclinic": "s from 1 to N - 1; s p = +c or -c (mod N) on every plane, "
                                   "c not 0; a line turned by 180 degrees is a turn of N / 2",
                      "geometry_read_for": "the smallest s only; the others are listed",
                      "hopf": "principal angles between the planes of every pair of "
                              "sub-orbit circles; equal within 1e-9 in cosine",
                      "r7_chains": f"the {r7['n_found']} listed by helix_setting, N <= 40",
                      "exact_field": "Q(sqrt2, sqrt3)"},
            "tags": tags, "exact_n24": exact_24(), "chain_7_13": chain_7_13(share7),
            "isoclinic_logged": iso, "isoclinic_r7": iso7,
            "summary": {
                "tags_changed": [[t["dimension_logged"], t["simplex_m"], t["n_to_close"],
                                  t["turns"], t["tag_2026_11_16_as_ordered"], t["tag"],
                                  t["relation"]] for t in tags if t["tag_changed"]],
                "setting_in": sum(1 for t in tags if t["tag"] == "SETTING_IN"),
                "setting_in_holding_at_1e-12": sum(
                    1 for t in tags if t["eq_5_1"]
                    and t["eq_5_1"]["holds_on_every_plane_at_1e-12"]),
                "setting_outside": sum(1 for t in tags if t["tag"] == "SETTING_OUTSIDE"),
                "logged_with_an_isoclinic_power": sum(1 for r in iso
                                                      if r["has_an_isoclinic_power"]),
                "logged_without": [[r["dimension_logged"], r["simplex_m"], r["n"], r["turns"]]
                                   for r in iso if not r["has_an_isoclinic_power"]],
                "survivors_with": [[r["n"], r["turns"], r["smallest_s"]] for r in iso
                                   if r["survivor_unmatched"] and r["has_an_isoclinic_power"]],
                "survivors_without": [[r["n"], r["turns"]] for r in iso
                                      if r["survivor_unmatched"]
                                      and not r["has_an_isoclinic_power"]],
                "r7_chains": len(iso7), "r7_with_an_isoclinic_power": len(with7),
                "r7_with": [[r["n"], r["turns"], r["smallest_s"],
                             r["smallest"]["polygon"]] for r in with7],
                "r7_without": len(iso7) - len(with7)}}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    comp = emlib.compressibility(
        ["y" if r["has_an_isoclinic_power"] else "n" for r in out["isoclinic_logged"]],
        "one letter per logged chain: y a power of the step is isoclinic, n none is")
    angles = {f"{r['n']}:{r['turns']}": (r["smallest"]["angle_deg"] if r["smallest"] else None)
              for r in out["isoclinic_logged"]}
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME, seed, {"sources": ["runs/helix_setting_seed0.json.gz",
                                                 "runs/corkscrew_seed0.json",
                                                 "runs/corkscrew6_seed0.json",
                                                 "runs/coxeter_exponents_seed0.json"],
                                     "tol": TOL, "status_before_run": "PROPOSED"},
                        out, angles, None, comp, None)
    print(p)
    for t in out["tags"]:
        e = t["eq_5_1"]
        print(f"{t['dimension_logged']}D m={t['simplex_m']} N={t['n_to_close']:2d} {str(t['turns']):12s}"
              f" span {t['points_span_dimension']} {t['relation']:18s} {t['tag']:15s} was "
              f"{t['tag_2026_11_16_as_ordered']}"
              + (f" | w {e['w']:.6f} resid {e['largest_abs_residual']:.1e} "
                 f"1e-12 {e['holds_on_every_plane_at_1e-12']}" if e else ""))
    x = out["exact_n24"]
    print("EXACT", {k: x[k] for k in x if k not in ("eq_5_1_planes",)})
    for q in x["eq_5_1_planes"]:
        print("   ", q)
    c = out["chain_7_13"]
    print("7,13", {k: c[k] for k in c if k not in ("orbits_of_the_600_cells",
                                                    "cosine_by_step_gap")})
    for r in c["orbits_of_the_600_cells"]:
        print("   ", r["orbit"], r["squared_length"], r["step_deg"], r["matches_taken_k_at_a_time"])
    for r in out["isoclinic_logged"]:
        s = r["smallest"]
        print(f"ISO {r['dimension_logged']}D m={r['simplex_m']} N={r['n']:2d} {str(r['turns']):12s} "
              f"{'*' if r['survivor_unmatched'] else ' '} powers "
              f"{[(q['s'], q['c'], q['sub_orbit_size']) for q in r['all_isoclinic_powers']]}"
              + (f" | {s['n_sub_orbits']} x {s['polygon']} great-circle "
                 f"{s['regular_polygons_on_great_circles']} parallel "
                 f"{s['hopf']['pairwise_clifford_parallel']} dist "
                 f"{s['hopf'].get('distances_between_circles_deg')} unequal "
                 f"{s['hopf'].get('principal_angles_of_unequal_pairs_deg')} signs "
                 f"{isoclinic_powers(r['n'], r['turns'])[0]['signs']}" if s else ""))
    sm = out["summary"]
    print("SUMMARY", {k: sm[k] for k in sm if k != "r7_with"})
    for r in out["isoclinic_r7"]:
        if r["has_an_isoclinic_power"]:
            print("R7", r["n"], r["turns"], r["all_isoclinic_powers"])


if __name__ == "__main__":
    main(sys.argv[1:])
