"""ag. corkscrew -- double rotation: does the shape come back, when, and how
does a wrong angle show?

CC0. Co-authored. Stdlib only. Carried method; status PROPOSED.

  1 corkscrew    a shape in R^4 or R^5 is turned in two planes at once, by
                 angle alpha in plane (0, 1) and beta in plane (2, 3), step
                 after step, until it is back or N steps are used. (In R^5
                 the fifth axis is fixed.)
  2 closure      the miss: how far the points are from where they started.
                 point return = every point back; set return = the shape
                 back on itself as a set. Exact where the shape allows it
                 (exact/corkscrew.py: turns by multiples of 36 degrees).
  3 sweep        alpha fixed, beta = alpha * r for ratios r: small
                 rationals, rational approximations, and irrationals. Does
                 it close only at rational ratios? (PROPOSED)
  4 error test   one angle is perturbed slightly. Read at the nominal return
                 step: how large is the miss, how does it grow with the
                 perturbation and with the number of laps, and which plane
                 carries it? (the "did it right" signal)
  5 tetrahedra   the Boerdijk-Coxeter chain: regular tetrahedra face to
                 face. In R^3 the chain is a helix that never returns.
                 CARRIED, to verify: in R^4 a chain of 30 closes. Here every
                 closed chain of N <= N_MAX regular tetrahedra lying on a
                 double-rotation orbit is searched for.

The step-k position is computed from the angle k * alpha directly, not by
k multiplications, so the miss is not built up from rounding.

Observer slots null.

  python3 experiments/corkscrew.py [seed]
  python3 experiments/corkscrew.py --small [seed]   # phone bound; own run file
"""
import math
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from exact import corkscrew as XK  # noqa: E402
from toolkit import primitives as P  # noqa: E402

NAME = "corkscrew"
N_STEPS = 6000
ALPHA = 36.0
RETURN_TOL = 1e-9
PHI = (1 + 5 ** 0.5) / 2
RATIOS = (("1", Fraction(1), None), ("2", Fraction(2), None), ("3", Fraction(3), None),
          ("1/2", Fraction(1, 2), None), ("3/2", Fraction(3, 2), None),
          ("5/3", Fraction(5, 3), None), ("7/4", Fraction(7, 4), None),
          ("13/8", Fraction(13, 8), None), ("1618/1000", Fraction(1618, 1000), None),
          ("3/2 + 1/10000", Fraction(3, 2) + Fraction(1, 10000), None),
          ("phi", None, PHI), ("sqrt(2)", None, 2 ** 0.5), ("pi/3", None, math.pi / 3),
          ("e/2", None, math.e / 2))
SHAPES = (("5-cell", 4, lambda: P.simplex(4)), ("16-cell", 4, lambda: P.cross_polytope(4)),
          ("tesseract", 4, lambda: P.hypercube(4)), ("24-cell", 4, P.cell24),
          ("5-simplex", 5, lambda: P.simplex(5)), ("5-orthoplex", 5, lambda: P.cross_polytope(5)),
          ("5-cube", 5, lambda: P.hypercube(5)))
SWEEP_SHAPES = ("tesseract", "5-cube")
EPSILONS = (1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
LAPS = (1, 2, 3, 5, 10)
N_MAX_CHAIN = 60
SMALL = {"n_steps": 1500, "n_max_chain": 32, "sweep_shapes": ("tesseract",)}
ASKED = ("PROPOSED, carried: a double rotation closes only at rational ratios of its two "
         "angles; a slightly wrong angle shows as a miss that grows in a readable way; a "
         "chain of regular tetrahedra that never closes in R^3 closes in R^4 at 30.")


def turned(verts, alpha_deg, beta_deg, k):
    a, b = math.radians(alpha_deg * k), math.radians(beta_deg * k)
    ca, sa, cb, sb = math.cos(a), math.sin(a), math.cos(b), math.sin(b)
    out = []
    for v in verts:
        w = list(v)
        w[0], w[1] = ca * v[0] - sa * v[1], sa * v[0] + ca * v[1]
        w[2], w[3] = cb * v[2] - sb * v[3], sb * v[2] + cb * v[3]
        out.append(w)
    return out


def point_miss(verts, cur, radius):
    return max(math.dist(a, b) for a, b in zip(verts, cur)) / radius


def set_miss(verts, cur, radius):
    return max(min(math.dist(a, b) for b in verts) for a in cur) / radius


def follow(verts, alpha, beta, n_steps, want_set=True):
    radius = max(P.norm(v) for v in verts)
    point_step = set_step = None
    best = (float("inf"), None)
    records, head = [], []
    for k in range(1, n_steps + 1):
        cur = turned(verts, alpha, beta, k)
        m = point_miss(verts, cur, radius)
        if k <= 24:
            head.append(m)
        if m < best[0]:
            best = (m, k)
            records.append([k, m])
        if want_set and set_step is None and set_miss(verts, cur, radius) < RETURN_TOL:
            set_step = k
        if m < RETURN_TOL:
            point_step = k
            break
    return {"point_return_step": point_step, "set_return_step": set_step,
            "smallest_miss": best[0], "smallest_miss_step": best[1],
            "steps_run": k, "record_low_misses": records[-16:],
            "miss_first_24_steps": head}


def predicted_return(alpha_deg, ratio):
    """Steps for both planes to complete whole turns, for a rational alpha and
    ratio (Fractions): the least k with k*alpha/360 and k*alpha*ratio/360 whole."""
    a = Fraction(alpha_deg).limit_denominator(10 ** 6) / 360
    b = a * ratio
    k = a.denominator * b.denominator // math.gcd(a.denominator, b.denominator)
    return k


def shapes_part(n_steps):
    """Every shape under two fixed double rotations."""
    out = []
    for name, dim, make in SHAPES:
        verts = make()
        rows = []
        for alpha, beta in ((90.0, 90.0), (72.0, 144.0), (36.0, 36.0 * PHI)):
            f = follow(verts, alpha, beta, min(n_steps, 400))
            rows.append({"alpha_deg": alpha, "beta_deg": beta, **f})
        out.append({"shape": name, "dimension": dim, "n_vertices": len(verts),
                    "planes": [[0, 1], [2, 3]],
                    "fixed_axes": [4] if dim == 5 else [], "rotations": rows})
    return out


def sweep_part(n_steps, sweep_shapes):
    shapes = {n: m() for n, _, m in SHAPES}
    out = []
    for sname in sweep_shapes:
        rows = []
        for label, frac, irr in RATIOS:
            r = float(frac) if frac is not None else irr
            f = follow(shapes[sname], ALPHA, ALPHA * r, n_steps, want_set=False)
            pred = predicted_return(ALPHA, frac) if frac is not None else None
            rows.append({"ratio": label, "ratio_value": r, "rational": frac is not None,
                         "predicted_point_return_step": pred,
                         "predicted_within_the_steps_run":
                             (pred <= n_steps) if pred is not None else None,
                         "closes": f["point_return_step"] is not None,
                         "return_matches_prediction":
                             (f["point_return_step"] == pred) if pred is not None
                             and pred <= n_steps else None,
                         **{k: f[k] for k in ("point_return_step", "smallest_miss",
                                              "smallest_miss_step", "steps_run",
                                              "record_low_misses")}})
        rat = [r for r in rows if r["rational"]]
        irr_rows = [r for r in rows if not r["rational"]]
        out.append({"shape": sname, "alpha_deg": ALPHA, "n_steps": n_steps, "rows": rows,
                    "rational_ratios": len(rat),
                    "rational_closing": sum(1 for r in rat if r["closes"]),
                    "rational_not_closing_within_steps":
                        [r["ratio"] for r in rat if not r["closes"]],
                    "irrational_ratios": len(irr_rows),
                    "irrational_closing": sum(1 for r in irr_rows if r["closes"]),
                    "smallest_miss_among_irrationals":
                        min(r["smallest_miss"] for r in irr_rows)})
    return out


def error_part():
    """Nominal (36, 72) degrees returns at step 10. One angle is off by eps."""
    verts = P.hypercube(4)
    radius = max(P.norm(v) for v in verts)
    nominal = (36.0, 72.0)
    K = 10
    out = []
    for which in ("alpha", "beta"):
        rows = []
        for eps in EPSILONS:
            a = nominal[0] + (eps if which == "alpha" else 0.0)
            b = nominal[1] + (eps if which == "beta" else 0.0)
            laps = []
            for lap in LAPS:
                cur = turned(verts, a, b, K * lap)
                m = point_miss(verts, cur, radius)
                # how much of the squared miss lies in each plane
                s01 = sum((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 for p, q in zip(verts, cur))
                s23 = sum((p[2] - q[2]) ** 2 + (p[3] - q[3]) ** 2 for p, q in zip(verts, cur))
                tot = s01 + s23
                rho = max(math.hypot(v[0], v[1]) if which == "alpha" else math.hypot(v[2], v[3])
                          for v in verts) / radius
                linear = K * lap * math.radians(eps) * rho
                laps.append({"laps": lap, "step": K * lap, "miss": m,
                             "linear_estimate": linear, "miss_over_linear_estimate": m / linear,
                             "share_in_plane_01": s01 / tot if tot else None,
                             "share_in_plane_23": s23 / tot if tot else None,
                             "plane_carrying_the_miss":
                                 ("01" if s01 > s23 else "23") if tot else None})
            rows.append({"eps_deg": eps, "laps": laps})
        one = [r["laps"][0]["miss"] for r in rows]
        slope = (math.log(one[-1] / one[0])
                 / math.log(EPSILONS[-1] / EPSILONS[0]))
        out.append({"perturbed": which, "nominal_deg": list(nominal), "return_step": K,
                    "rows": rows,
                    "log_log_slope_of_miss_against_eps": slope,
                    "miss_over_linear_estimate_range":
                        [min(l["miss_over_linear_estimate"] for r in rows for l in r["laps"]),
                         max(l["miss_over_linear_estimate"] for r in rows for l in r["laps"])],
                    "plane_carrying_the_miss": sorted({l["plane_carrying_the_miss"]
                                                       for r in rows for l in r["laps"]}),
                    "expected_plane": "01" if which == "alpha" else "23"})
    cur = turned(verts, nominal[0], nominal[1], K)
    return {"shape": "tesseract", "nominal_miss_at_return": point_miss(verts, cur, radius),
            "linear_estimate": "laps * return_step * eps (radians) * (largest radius in the "
                               "perturbed plane / radius)",
            "perturbations": out}


def chain_part(n_max):
    """Regular-tetrahedra chains on a double-rotation orbit in R^4:
    v_k = (r1 cos k a, r1 sin k a, r2 cos k b, r2 sin k b), a = 360 p / N,
    b = 360 q / N. Four consecutive points are a regular tetrahedron when
    d1 = d2 = d3, with d_j^2 = t (1 - cos j a) + (1 - t)(1 - cos j b) (up to a
    factor) and t = r1^2 / (r1^2 + r2^2). d1 = d2 fixes t; d3 is then checked."""
    found = []
    for n in range(4, n_max + 1):
        for p in range(1, n // 2 + 1):
            for q in range(p + 1, n // 2 + 1):
                if math.gcd(math.gcd(p, q), n) != 1:
                    continue
                a, b = 2 * math.pi * p / n, 2 * math.pi * q / n
                fa = [1 - math.cos(j * a) for j in (1, 2, 3)]
                fb = [1 - math.cos(j * b) for j in (1, 2, 3)]
                den = (fa[0] - fb[0]) - (fa[1] - fb[1])
                if abs(den) < 1e-12:
                    continue
                t = (fb[1] - fb[0]) / den
                if not 1e-9 < t < 1 - 1e-9:
                    continue
                d = [t * x + (1 - t) * y for x, y in zip(fa, fb)]
                if abs(d[0] - d[2]) > 1e-9 * d[0]:
                    continue
                edge = math.sqrt(2 * d[0])          # for r1^2 + r2^2 = 1
                found.append({"n_tetrahedra_to_close": n, "p": p, "q": q,
                              "alpha_deg": 360.0 * p / n, "beta_deg": 360.0 * q / n,
                              "r1_squared_share": t, "edge_over_radius": edge,
                              "d1_d2_d3_squared": d,
                              "largest_relative_edge_difference":
                                  max(abs(x - d[0]) for x in d) / d[0]})
    # the same chain in R^3: v_k = (r cos k th, r sin k th, k h); d1 = d2 = d3
    th = math.acos(-2.0 / 3.0)
    r, h = 3 * math.sqrt(3) / 10, 1 / math.sqrt(10)
    pts = [(r * math.cos(k * th), r * math.sin(k * th), k * h) for k in range(4)]
    d3 = [math.dist(pts[0], pts[j]) for j in (1, 2, 3)]
    best = (float("inf"), None)
    records = []
    for k in range(1, 2001):
        m = abs((k * math.degrees(th) + 180.0) % 360.0 - 180.0)
        if m < best[0]:
            best = (m, k)
            records.append([k, m])
    by_n = {}
    for f in found:
        by_n.setdefault(f["n_tetrahedra_to_close"], []).append([f["p"], f["q"]])
    return {"r4": {"searched": f"every N from 4 to {n_max}, every pair of turns "
                               "(360 p / N, 360 q / N), 1 <= p < q <= N / 2, coprime with N",
                   "closed_chains": found, "n_values_with_a_closed_chain": sorted(by_n),
                   "turn_pairs_by_n": {str(k): v for k, v in by_n.items()},
                   "closes_at_30": 30 in by_n,
                   "edge_over_radius_at_30":
                       [f["edge_over_radius"] for f in found
                        if f["n_tetrahedra_to_close"] == 30],
                   "one_over_phi": 1 / PHI},
            "r3": {"helix": "v_k = (r cos k theta, r sin k theta, k h), theta = acos(-2/3), "
                            "r = 3 sqrt(3) / 10, h = 1 / sqrt(10), edge 1",
                   "twist_deg": math.degrees(th),
                   "edges_d1_d2_d3": d3,
                   "moves_along_its_axis_every_step": h,
                   "returns_to_its_start": False,
                   "twist_back_to_zero_within_2000_steps": best[0] < 1e-9,
                   "smallest_twist_miss_deg": best[0], "at_step": best[1],
                   "record_low_twist_misses": records[-10:]}}


def exact_part():
    f = XK.five_cell()
    duo = XK.duoprism_5_5()
    rows = [XK.closure(f, p, q) for p, q in ((2, 4), (1, 2), (2, 2), (1, 3), (4, 3), (5, 5),
                                              (1, 1), (3, 7))]
    return {"five_cell": {"vertices": "(zeta^k, zeta^(2k)), k = 0..4",
                          "is_a_regular_simplex_exactly": XK.is_regular_simplex(f),
                          "closures": rows},
            "duoprism_5_5": {"n_points": len(duo),
                             "closures": [XK.closure(duo, p, q) for p, q in ((2, 4), (1, 3))]},
            "limit": "turns by multiples of 36 degrees only; an irrational ratio cannot be "
                     "written in these forms, so the exact part reads rational turns only"}


def run(seed=0, n_steps=N_STEPS, n_max_chain=N_MAX_CHAIN, sweep_shapes=SWEEP_SHAPES):
    return {"asked": ASKED,
            "bound": {"seed": seed, "n_steps": n_steps, "alpha_deg": ALPHA,
                      "return_tol": RETURN_TOL, "n_max_chain": n_max_chain,
                      "miss": "largest distance of a point from where it started, over the "
                              "shape's radius",
                      "positions": "from the angle k * alpha directly, not by k products",
                      "irrationals_are_floats": "an irrational ratio is held as a 64-bit float, "
                                                "which is itself a rational; within N steps it "
                                                "does not return",
                      "not_closing_means": "not within the steps run"},
            "shapes": shapes_part(n_steps),
            "sweep": sweep_part(n_steps, sweep_shapes),
            "error_test": error_part(),
            "tetrahedra_chain": chain_part(n_max_chain),
            "exact": exact_part()}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = (run(seed, SMALL["n_steps"], SMALL["n_max_chain"], SMALL["sweep_shapes"])
           if small else run(seed))
    symbols = ["c" if r["closes"] else "o" for s in out["sweep"] for r in s["rows"]]
    comp = emlib.compressibility(
        symbols, "one letter per ratio in the sweep: c closes within the steps run, o does "
                 "not; shape by shape")
    angles = {str(f["n_tetrahedra_to_close"]) + f":{f['p']},{f['q']}":
              [f["alpha_deg"], f["beta_deg"]]
              for f in out["tetrahedra_chain"]["r4"]["closed_chains"]}
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "n_steps": out["bound"]["n_steps"],
                         "alpha_deg": ALPHA, "ratios": [r[0] for r in RATIOS],
                         "epsilons_deg": list(EPSILONS), "laps": list(LAPS),
                         "n_max_chain": out["bound"]["n_max_chain"], "asked": ASKED,
                         "status_before_run": "PROPOSED"},
                        out, angles, None, comp, None)
    print(p)
    print("SHAPES (72, 144):", [(s["shape"], s["rotations"][1]["set_return_step"],
                                 s["rotations"][1]["point_return_step"]) for s in out["shapes"]])
    print("SHAPES (90, 90):", [(s["shape"], s["rotations"][0]["set_return_step"],
                                s["rotations"][0]["point_return_step"]) for s in out["shapes"]])
    for s in out["sweep"]:
        print("SWEEP", s["shape"], "rational closing", s["rational_closing"], "of",
              s["rational_ratios"], "| not closing:", s["rational_not_closing_within_steps"],
              "| irrational closing", s["irrational_closing"], "of", s["irrational_ratios"])
        for r in s["rows"]:
            print(f"   {r['ratio']:14s} predicted {r['predicted_point_return_step']}  returned "
                  f"{r['point_return_step']}  smallest miss {r['smallest_miss']:.3e} at "
                  f"{r['smallest_miss_step']}")
    for e in out["error_test"]["perturbations"]:
        print("ERROR", e["perturbed"], "slope", round(e["log_log_slope_of_miss_against_eps"], 6),
              "miss/linear", [round(x, 6) for x in e["miss_over_linear_estimate_range"]],
              "plane", e["plane_carrying_the_miss"], "expected", e["expected_plane"])
    c = out["tetrahedra_chain"]
    print("CHAIN R4 closes at N =", c["r4"]["n_values_with_a_closed_chain"], "pairs",
          c["r4"]["turn_pairs_by_n"], "| at 30:", c["r4"]["closes_at_30"], "edge/radius",
          c["r4"]["edge_over_radius_at_30"], "1/phi", c["r4"]["one_over_phi"])
    print("CHAIN R3 twist", c["r3"]["twist_deg"], "edges", c["r3"]["edges_d1_d2_d3"],
          "smallest twist miss", c["r3"]["smallest_twist_miss_deg"], "at", c["r3"]["at_step"])
    print("EXACT 5-cell regular:", out["exact"]["five_cell"]["is_a_regular_simplex_exactly"],
          [(x["turns_of_36_deg"], x["set_return_step"], x["point_return_step"])
           for x in out["exact"]["five_cell"]["closures"]])


if __name__ == "__main__":
    main(sys.argv[1:])
