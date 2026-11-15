"""aj. corkscrew6 -- triple rotation in R^6: return, miss, partial closure, chains.

CC0. Co-authored. Stdlib only. Status: PROPOSED. Extends
experiments/corkscrew.py to three planes; that file and its runs are not
changed.

  1 triple rotation   a shape in R^6 is turned by alpha, beta, gamma in the
                      planes (0,1), (2,3), (4,5), step after step. Shapes:
                      6-simplex, 6-cube, 6-orthoplex.
  2 closure           return step and smallest miss. For rational turns the
                      predicted return is the least common multiple of the
                      three planes' periods. Exact where the angles allow
                      (exact/corkscrew.py, turns by multiples of 36 degrees).
  3 sweep             rate triples a : b : c with all three rational, one
                      irrational, two irrational. Point return and set
                      return are logged separately.
  4 miss              one plane's angle off by eps: is the miss confined to
                      that plane? Two planes off: does each plane carry what
                      it carries alone?
  5 partial closure   two rates rational, one not: at the step where the two
                      rational planes are back, how far is each plane?
  6 chains            closed chains of regular simplices on a triple-rotation
                      orbit, v_k = (r1 e^(ik a), r2 e^(ik b), r3 e^(ik c)),
                      a, b, c = 360 p/N, 360 q/N, 360 r/N, N up to N_MAX.
                      m + 1 consecutive points form a regular m-simplex when
                      d_1 = ... = d_m.

Positions come from the angle k * alpha directly, not from k products.
Observer slots null.

  python3 experiments/corkscrew6.py --small [seed]   # phone bound; own run file
  python3 experiments/corkscrew6.py [seed]           # full (chains to N = 120: ~1 min)
"""
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from exact import corkscrew as XK  # noqa: E402
from toolkit import primitives as P  # noqa: E402

NAME = "corkscrew6"
N_STEPS, N_MAX_CHAIN = 6000, 120
ALPHA = 36.0
RETURN_TOL = 1e-9
SET_DECIMALS = 7              # a set is back when its points, rounded to this, are the same set
PHI, R2, PI3 = (1 + 5 ** 0.5) / 2, 2 ** 0.5, math.pi / 3
F = Fraction
TRIPLES = (("1 : 1 : 1", (F(1), F(1), F(1)), None, "all rational"),
           ("1 : 2 : 3", (F(1), F(2), F(3)), None, "all rational"),
           ("1 : 2 : 4", (F(1), F(2), F(4)), None, "all rational"),
           ("2 : 3 : 5", (F(2), F(3), F(5)), None, "all rational"),
           ("1 : 3/2 : 5/3", (F(1), F(3, 2), F(5, 3)), None, "all rational"),
           ("1 : 13/8 : 7/4", (F(1), F(13, 8), F(7, 4)), None, "all rational"),
           ("1 : 2 : phi", (F(1), F(2), None), (1.0, 2.0, PHI), "one irrational"),
           ("1 : 3/2 : sqrt(2)", (F(1), F(3, 2), None), (1.0, 1.5, R2), "one irrational"),
           ("1 : 1 : pi/3", (F(1), F(1), None), (1.0, 1.0, PI3), "one irrational"),
           ("1 : phi : sqrt(2)", (F(1), None, None), (1.0, PHI, R2), "two irrational"),
           ("1 : sqrt(2) : pi/3", (F(1), None, None), (1.0, R2, PI3), "two irrational"))
SHAPES = (("6-simplex", lambda: P.simplex(6)), ("6-orthoplex", lambda: P.cross_polytope(6)),
          ("6-cube", lambda: P.hypercube(6)))
EPSILONS = (1e-6, 1e-4, 1e-2)
SMALL = {"n_steps": 1200, "n_max_chain": 40, "shapes": ("6-simplex", "6-orthoplex")}
ASKED = ("PROPOSED: a triple rotation returns at the least common multiple of its three "
         "periods and only for rational rates; a wrong angle shows in its own plane; with one "
         "irrational rate the orbit closes in four dimensions and not in the other two.")


def turned(verts, angles_deg, k):
    cs = [(math.cos(math.radians(a * k)), math.sin(math.radians(a * k))) for a in angles_deg]
    out = []
    for v in verts:
        w = list(v)
        for i, (c, s) in enumerate(cs):
            x, y = v[2 * i], v[2 * i + 1]
            w[2 * i], w[2 * i + 1] = c * x - s * y, s * x + c * y
        out.append(w)
    return out


def plane_misses(verts, cur, radius):
    """Largest displacement of a point within each plane, over the radius."""
    return [max(math.hypot(a[2 * i] - b[2 * i], a[2 * i + 1] - b[2 * i + 1])
                for a, b in zip(verts, cur)) / radius for i in range(3)]


def key(points):
    return frozenset(tuple(round(x, SET_DECIMALS) + 0.0 for x in p) for p in points)


def follow(verts, angles, n_steps):
    radius = max(P.norm(v) for v in verts)
    start_key = key(verts)
    point_step = set_step = None
    plane_first = [None, None, None]
    best = (float("inf"), None)
    records = []
    for k in range(1, n_steps + 1):
        cur = turned(verts, angles, k)
        pm = plane_misses(verts, cur, radius)
        m = max(math.dist(a, b) for a, b in zip(verts, cur)) / radius
        for i in range(3):
            if plane_first[i] is None and pm[i] < RETURN_TOL:
                plane_first[i] = k
        if m < best[0]:
            best = (m, k)
            records.append([k, m])
        if set_step is None and key(cur) == start_key:
            set_step = k
        if m < RETURN_TOL:
            point_step = k
            break
    return {"point_return_step": point_step, "set_return_step": set_step,
            "first_return_step_per_plane": plane_first,
            "smallest_miss": best[0], "smallest_miss_step": best[1], "steps_run": k,
            "record_low_misses": records[-12:]}


def period(alpha_deg, rate):
    """Steps for one plane to complete whole turns (rate a Fraction)."""
    return (Fraction(alpha_deg).limit_denominator(10 ** 6) * rate / 360).denominator


def lcm(values):
    out = 1
    for v in values:
        out = out * v // math.gcd(out, v)
    return out


def sweep_part(n_steps, shape_names):
    shapes = dict((n, m()) for n, m in SHAPES)
    out = []
    for sname in shape_names:
        rows = []
        for label, fr, irr, kind in TRIPLES:
            rates = [float(x) for x in fr] if irr is None else list(irr)
            angles = [ALPHA * r for r in rates]
            f = follow(shapes[sname], angles, n_steps)
            periods = [period(ALPHA, x) if x is not None else None for x in fr]
            pred = lcm(periods) if all(p is not None for p in periods) else None
            rational_planes = [i for i, x in enumerate(fr) if x is not None]
            partial = lcm([periods[i] for i in rational_planes])
            rows.append({"rates": label, "kind": kind, "angles_deg": angles,
                         "period_per_plane": periods, "predicted_point_return_lcm": pred,
                         "point_return_matches_lcm":
                             (f["point_return_step"] == pred) if pred is not None
                             and pred <= n_steps else None,
                         "closes": f["point_return_step"] is not None,
                         "set_returns": f["set_return_step"] is not None,
                         "lcm_of_the_rational_planes": partial, **f})
        kinds = {}
        for r in rows:
            k = kinds.setdefault(r["kind"], {"n": 0, "point_return": 0, "set_return": 0})
            k["n"] += 1
            k["point_return"] += r["closes"]
            k["set_return"] += r["set_returns"]
        out.append({"shape": sname, "n_vertices": len(shapes[sname]), "alpha_deg": ALPHA,
                    "n_steps": n_steps, "rows": rows, "by_kind": kinds})
    return out


def miss_part():
    """Nominal rates 1 : 2 : 3 at 36 degrees: returns at step 10."""
    verts = P.hypercube(6)
    radius = max(P.norm(v) for v in verts)
    nominal = [36.0, 72.0, 108.0]
    K = 10
    rho = [max(math.hypot(v[2 * i], v[2 * i + 1]) for v in verts) / radius for i in range(3)]

    def at(eps):
        cur = turned(verts, [a + e for a, e in zip(nominal, eps)], K)
        return plane_misses(verts, cur, radius)
    single, double = [], []
    for e in EPSILONS:
        for i in range(3):
            eps = [0.0, 0.0, 0.0]
            eps[i] = e
            pm = at(eps)
            lin = K * math.radians(e) * rho[i]
            single.append({"eps_deg": e, "plane_perturbed": i, "miss_per_plane": pm,
                           "miss_in_the_other_planes": max(pm[j] for j in range(3) if j != i),
                           "confined_to_its_plane":
                               max(pm[j] for j in range(3) if j != i) < 1e-12,
                           "linear_estimate": lin, "miss_over_linear_estimate": pm[i] / lin})
        for i, j in ((0, 1), (0, 2), (1, 2)):
            eps = [0.0, 0.0, 0.0]
            eps[i], eps[j] = e, 2 * e
            pm = at(eps)
            alone_i, alone_j = [0.0] * 3, [0.0] * 3
            alone_i[i], alone_j[j] = e, 2 * e
            mi, mj = at(alone_i)[i], at(alone_j)[j]
            k = [x for x in range(3) if x not in (i, j)][0]
            double.append({"eps_deg": [e, 2 * e], "planes_perturbed": [i, j],
                           "miss_per_plane": pm,
                           "each_plane_alone": [mi, mj],
                           "difference_from_alone": [pm[i] - mi, pm[j] - mj],
                           "same_as_alone": abs(pm[i] - mi) < 1e-12 and abs(pm[j] - mj) < 1e-12,
                           "untouched_plane_miss": pm[k]})
    return {"shape": "6-cube", "nominal_deg": nominal, "return_step": K,
            "nominal_miss_per_plane": at([0.0, 0.0, 0.0]),
            "one_plane": single, "two_planes": double,
            "all_single_confined": all(s["confined_to_its_plane"] for s in single),
            "all_double_same_as_alone": all(d["same_as_alone"] for d in double),
            "miss_over_linear_estimate_range":
                [min(s["miss_over_linear_estimate"] for s in single),
                 max(s["miss_over_linear_estimate"] for s in single)]}


def partial_part(n_steps):
    """Rates 1 : 2 : r3 with r3 irrational: read each plane at the steps where
    the two rational planes are back."""
    verts = P.hypercube(6)
    radius = max(P.norm(v) for v in verts)
    out = []
    for label, r3 in (("phi", PHI), ("sqrt(2)", R2), ("pi/3", PI3)):
        angles = [36.0, 72.0, 36.0 * r3]
        K = 10
        rows = []
        for lap in (1, 2, 3, 5, 10, 50, 100):
            if K * lap > n_steps:
                continue
            pm = plane_misses(verts, turned(verts, angles, K * lap), radius)
            rows.append({"step": K * lap, "miss_per_plane": pm,
                         "rational_planes_back": pm[0] < RETURN_TOL and pm[1] < RETURN_TOL,
                         "third_plane_miss": pm[2]})
        f = follow(verts, angles, n_steps)
        best3 = min((plane_misses(verts, turned(verts, angles, K * lap), radius)[2], K * lap)
                    for lap in range(1, n_steps // K + 1))
        out.append({"rates": f"1 : 2 : {label}", "angles_deg": angles,
                    "return_step_of_the_two_rational_planes": K, "at_laps": rows,
                    "rational_planes_back_at_every_lap":
                        all(r["rational_planes_back"] for r in rows),
                    "third_plane_ever_back_at_a_lap": best3[0] < RETURN_TOL,
                    "third_plane_smallest_miss_at_a_lap": best3[0],
                    "third_plane_smallest_miss_step": best3[1],
                    "first_return_step_per_plane": f["first_return_step_per_plane"],
                    "full_point_return_step": f["point_return_step"],
                    "full_smallest_miss": f["smallest_miss"]})
    return {"shape": "6-cube", "closes_in": "planes (0,1) and (2,3): four of six dimensions",
            "cases": out}


def chain_part(n_max):
    """Closed chains of regular m-simplices, m = 3..6, on a triple-rotation orbit."""
    found = {m: [] for m in (3, 4, 5, 6)}
    lower = {m: 0 for m in (3, 4, 5, 6)}
    tried = 0
    for n in range(4, n_max + 1):
        f = [[1.0 - math.cos(2 * math.pi * j * p / n) for j in range(7)] for p in range(n // 2 + 1)]
        for p in range(1, n // 2 + 1):
            for q in range(p + 1, n // 2 + 1):
                g = math.gcd(math.gcd(p, q), n)
                for r in range(q + 1, n // 2 + 1):
                    if math.gcd(g, r) != 1:
                        continue
                    tried += 1
                    fa, fb, fc = f[p], f[q], f[r]
                    # d_j = t1 fa_j + t2 fb_j + (1 - t1 - t2) fc_j ; d_2 = d_1 and d_3 = d_1
                    a11 = (fa[2] - fc[2]) - (fa[1] - fc[1])
                    a12 = (fb[2] - fc[2]) - (fb[1] - fc[1])
                    a21 = (fa[3] - fc[3]) - (fa[1] - fc[1])
                    a22 = (fb[3] - fc[3]) - (fb[1] - fc[1])
                    b1, b2 = fc[1] - fc[2], fc[1] - fc[3]
                    det = a11 * a22 - a12 * a21
                    if abs(det) < 1e-12:
                        continue
                    t1 = (b1 * a22 - a12 * b2) / det
                    t2 = (a11 * b2 - b1 * a21) / det
                    t3 = 1.0 - t1 - t2
                    if min(t1, t2, t3) < -1e-9:
                        continue
                    d = [t1 * fa[j] + t2 * fb[j] + t3 * fc[j] for j in range(7)]
                    if d[1] <= 1e-12:
                        continue
                    m = 3
                    while m < 6 and abs(d[m + 1] - d[1]) <= 1e-9 * d[1]:
                        m += 1
                    six_d = min(t1, t2, t3) > 1e-9
                    for mm in range(3, m + 1):
                        if six_d:
                            if mm >= 4 or len(found[3]) < 40:
                                found[mm].append({"n": n, "p": p, "q": q, "r": r,
                                                  "angles_deg": [360.0 * x / n for x in (p, q, r)],
                                                  "radius_squared_shares": [t1, t2, t3],
                                                  "edge_over_radius": math.sqrt(2 * d[1])})
                            if mm == 3 and len(found[3]) >= 40:
                                lower.setdefault("m3_more", 0)
                                lower["m3_more"] += 1
                        else:
                            lower[mm] += 1
    return {"searched": f"N from 4 to {n_max}; turns 360 p/N, 360 q/N, 360 r/N with "
                        "1 <= p < q < r <= N/2, coprime with N",
            "triples_tried": tried,
            "method": "d_1 = d_2 = d_3 fixes the three radius shares; larger simplices are "
                      "then checked (relative tolerance 1e-9)",
            "chains_using_all_three_planes": {
                str(m): {"n_found": len(found[m]) + (lower.get("m3_more", 0) if m == 3 else 0),
                         "n_values": sorted({c["n"] for c in found[m]}),
                         "listed": found[m]} for m in (3, 4, 5, 6)},
            "chains_with_a_zero_radius_plane": {str(m): lower[m] for m in (3, 4, 5, 6)},
            "note": "m = 3 (tetrahedra): two equations for two free shares, so a solution is "
                    "expected for many turns; only the first 40 are listed, all are counted"}


def exact_part():
    out = {}
    for name, pts in (("three pentagons in orthogonal planes", XK.three_pentagons()),
                      ("(zeta^k, zeta^2k, zeta^3k)", XK.lifted_orbit())):
        rows = []
        for turns in ((2, 4, 6), (1, 2, 3), (2, 2, 2), (2, 4, 5), (1, 5, 10), (4, 6, 8)):
            c = XK.closure_multi(pts, turns)
            periods = [10 // math.gcd(t, 10) for t in turns]
            c["predicted_lcm_of_plane_periods"] = lcm(periods)
            c["point_return_matches_lcm"] = c["point_return_step"] == lcm(periods)
            rows.append(c)
        out[name] = {"n_points": len(pts), "closures": rows,
                     "all_match_the_lcm": all(r["point_return_matches_lcm"] for r in rows)}
    return out


def run(seed=0, n_steps=N_STEPS, n_max_chain=N_MAX_CHAIN, shape_names=None):
    names = shape_names or [n for n, _ in SHAPES]
    return {"asked": ASKED,
            "bound": {"seed": seed, "n_steps": n_steps, "alpha_deg": ALPHA,
                      "return_tol": RETURN_TOL, "set_decimals": SET_DECIMALS,
                      "n_max_chain": n_max_chain, "shapes": list(names),
                      "planes": [[0, 1], [2, 3], [4, 5]],
                      "irrationals_are_floats": "held as 64-bit floats; within the steps run "
                                                "they do not return",
                      "not_closing_means": "not within the steps run",
                      "positions": "from the angle k * alpha directly"},
            "sweep": sweep_part(n_steps, names), "miss": miss_part(),
            "partial_closure": partial_part(n_steps), "chains": chain_part(n_max_chain),
            "exact": exact_part()}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = (run(seed, SMALL["n_steps"], SMALL["n_max_chain"], SMALL["shapes"]) if small
           else run(seed))
    symbols = [("p" if r["closes"] else "-") + ("s" if r["set_returns"] else "-")
               for s in out["sweep"] for r in s["rows"]]
    comp = emlib.compressibility(
        symbols, "per rate triple: p point return or -, then s set return or -; shape by shape")
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "n_steps": out["bound"]["n_steps"],
                         "n_max_chain": out["bound"]["n_max_chain"],
                         "shapes": out["bound"]["shapes"], "alpha_deg": ALPHA,
                         "rate_triples": [t[0] for t in TRIPLES], "asked": ASKED,
                         "status_before_run": "PROPOSED"},
                        out, None, None, comp, None)
    print(p)
    for s in out["sweep"]:
        print("SWEEP", s["shape"], s["by_kind"])
        for r in s["rows"]:
            print(f"   {r['rates']:18s} periods {r['period_per_plane']} lcm "
                  f"{r['predicted_point_return_lcm']}  point {r['point_return_step']}  set "
                  f"{r['set_return_step']}  planes {r['first_return_step_per_plane']}  miss "
                  f"{r['smallest_miss']:.2e}")
    m = out["miss"]
    print("MISS one plane confined:", m["all_single_confined"], "| two planes same as alone:",
          m["all_double_same_as_alone"], "| miss/linear",
          [round(x, 6) for x in m["miss_over_linear_estimate_range"]])
    for c in out["partial_closure"]["cases"]:
        print("PARTIAL", c["rates"], "rational planes back at every lap:",
              c["rational_planes_back_at_every_lap"], "| third plane at laps",
              [(r["step"], float(f"{r['third_plane_miss']:.3g}")) for r in c["at_laps"]],
              "| smallest", float(f"{c['third_plane_smallest_miss_at_a_lap']:.3g}"),
              "| full return", c["full_point_return_step"])
    ch = out["chains"]
    print("CHAINS tried", ch["triples_tried"], {m: (v["n_found"], v["n_values"][:12])
                                                for m, v in ch["chains_using_all_three_planes"].items()},
          "| zero-radius plane:", ch["chains_with_a_zero_radius_plane"])
    for m in ("4", "5", "6"):
        print("   m =", m, [(c["n"], c["p"], c["q"], c["r"], round(c["edge_over_radius"], 5))
                            for c in ch["chains_using_all_three_planes"][m]["listed"]][:12])
    for n, e in out["exact"].items():
        print("EXACT", n, "all match lcm:", e["all_match_the_lcm"],
              [(c["turns_of_36_deg"], c["set_return_step"], c["point_return_step"],
                c["plane_return_steps"]) for c in e["closures"]])


if __name__ == "__main__":
    main(sys.argv[1:])
