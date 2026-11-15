"""an. spherical_helix -- the four unmatched 5-cell chains, vertex by vertex,
beside the 30-tetrahedron ring.

CC0. Co-authored. Stdlib only. Status: PROPOSED.

The four closed chains of 5-cells logged by corkscrew6 (N = 18, 24, 28, 28)
matched no Coxeter group. Each is rebuilt from its logged turns and radius
shares, v_k = (r1 e^(ik a), r2 e^(ik b), r3 e^(ik c)), and read:

  1  do every 5 consecutive vertices form a regular 5-cell (all ten edges
     equal)? is the step v_k -> v_(k+1) one isometry for every k?
  2  the same readings on the ring of 30 tetrahedra in R^4 (corkscrew:
     N = 30, turns 1 and 11), and the two set side by side: simplex
     dimension against sphere dimension, shared facets, the angle between
     consecutive simplices at their shared facet, edge over radius, and
     each plane's turn beside the turns of the FLAT simplicial helix with
     the same n (roots of the bracket in REFERENCES R5, computed here).

Nothing is interpreted. Observer slots null.

  python3 experiments/spherical_helix.py [seed]
"""
import cmath
import itertools
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from toolkit import primitives as P  # noqa: E402

NAME = "spherical_helix"
TOL = 1e-9


def orbit(n, turns, shares):
    out = []
    for k in range(n):
        v = []
        for t, s in zip(turns, shares):
            a = 2 * math.pi * t * k / n
            v += [math.sqrt(s) * math.cos(a), math.sqrt(s) * math.sin(a)]
        out.append(v)
    return out


def step_matrix_apply(v, n, turns):
    w = []
    for i, t in enumerate(turns):
        a = 2 * math.pi * t / n
        x, y = v[2 * i], v[2 * i + 1]
        w += [math.cos(a) * x - math.sin(a) * y, math.sin(a) * x + math.cos(a) * y]
    return w


def rank(vectors, tol=1e-9):
    return len(P.gram_schmidt([list(v) for v in vectors], tol=tol))


def hinge_angle(pts, k, m, n):
    """Angle (degrees) between simplex k (vertices k..k+m) and simplex k+1 at the
    facet they share (vertices k+1..k+m): between the two apexes' heights
    above that facet."""
    facet = [pts[(k + j) % n] for j in range(1, m + 1)]
    base = facet[0]
    basis = P.gram_schmidt([P.sub(f, base) for f in facet[1:]], tol=1e-12)

    def height(apex):
        w = P.sub(apex, base)
        for b in basis:
            w = P.sub(w, P.scale(b, P.dot(w, b)))
        return w
    h1, h2 = height(pts[k % n]), height(pts[(k + m + 1) % n])
    c = P.dot(h1, h2) / (P.norm(h1) * P.norm(h2))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def flat_helix_turns(n):
    """Turns 2 theta (degrees) of the flat simplicial helix with n-vertex
    simplices: the unit roots lambda = exp(2 theta i) of
    sum_{k=1}^{n-1} k (n-k) lambda^(n-1-k) (REFERENCES R5, corrected form)."""
    coef = [k * (n - k) for k in range(1, n)]          # highest power first
    deg = n - 2

    def g(phi):
        lam = cmath.exp(1j * phi)
        p = sum(c * lam ** (deg - i) for i, c in enumerate(coef))
        return (cmath.exp(-1j * deg * phi / 2) * p).real      # real: the list is palindromic
    roots, steps = [], 20000
    prev_phi, prev = 1e-9, g(1e-9)
    for s in range(1, steps + 1):
        phi = math.pi * s / steps
        cur = g(phi)
        if prev == 0.0 or prev * cur < 0:
            lo, hi = prev_phi, phi
            for _ in range(80):
                mid = (lo + hi) / 2
                if g(lo) * g(mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            roots.append(math.degrees((lo + hi) / 2))
        prev_phi, prev = phi, cur
    if abs(g(math.pi)) < 1e-9 and not any(abs(r - 180.0) < 1e-6 for r in roots):
        roots.append(180.0)
    return {"n": n, "coefficients": coef, "turns_deg": sorted(roots),
            "n_turns_expected": (n - 2 + 1) // 2}


def read(label, n, turns, shares, m):
    pts = orbit(n, turns, shares)
    dim = len(pts[0])
    windows = []
    worst = 0.0
    for k in range(n):
        idx = [(k + j) % n for j in range(m + 1)]
        d = [math.dist(pts[a], pts[b]) for a, b in itertools.combinations(idx, 2)]
        worst = max(worst, (max(d) - min(d)) / min(d))
        windows.append(max(d) - min(d))
    edge = math.dist(pts[0], pts[1])
    d_by_gap = [math.dist(pts[0], pts[j % n]) for j in range(1, m + 3)]
    step_err = max(math.dist(step_matrix_apply(pts[k], n, turns), pts[(k + 1) % n])
                   for k in range(n))
    back = pts[0]
    for _ in range(n):
        back = step_matrix_apply(back, n, turns)
    hinges = [hinge_angle(pts, k, m, n) for k in range(n)]
    flat = flat_helix_turns(m + 1)
    turn_deg = [360.0 * t / n for t in turns]
    nearest = []
    for ft in flat["turns_deg"]:
        j = min(range(len(turn_deg)), key=lambda i: abs(turn_deg[i] - ft))
        nearest.append({"flat_turn_deg": ft, "nearest_chain_turn_deg": turn_deg[j],
                        "difference_deg": turn_deg[j] - ft})
    used = {min(range(len(turn_deg)), key=lambda i: abs(turn_deg[i] - ft))
            for ft in flat["turns_deg"]}
    left = [turn_deg[i] for i in range(len(turn_deg)) if i not in used]
    return {"chain": label, "n_to_close": n, "turns": list(turns), "turn_deg": turn_deg,
            "radius_squared_shares": list(shares), "simplex_m": m,
            "simplex_vertices": m + 1, "ambient_dimension": dim,
            "points_span_dimension": rank(pts),
            # the sphere the chain actually lies on: one less than what its points
            # span (a plane turned by 180 degrees holds its points on a line)
            "sphere_dimension": rank(pts) - 1,
            "sphere_dimension_minus_simplex_dimension": rank(pts) - 1 - m,
            "planes_turned_by_180_deg": [i for i, t in enumerate(turns) if 2 * t == n],
            "all_on_one_sphere": max(abs(P.norm(p) - 1.0) for p in pts) < TOL,
            "n_distinct_vertices": len({tuple(round(x, 9) for x in p) for p in pts}),
            "every_window_is_a_regular_simplex": worst < TOL,
            "windows_checked": n, "largest_relative_edge_spread": worst,
            "edge_over_radius": edge,
            "distance_by_step_gap": d_by_gap,
            "next_vertex_is_at_the_edge_distance": abs(d_by_gap[m] - edge) < TOL * edge,
            "step_is_one_isometry": step_err < TOL,
            "step_map": "the same rotation for every k: turns " + ", ".join(
                f"{a:.4f}" for a in turn_deg) + " degrees in three orthogonal planes"
            if len(turns) == 3 else "the same rotation for every k: turns " + ", ".join(
                f"{a:.4f}" for a in turn_deg) + " degrees in two orthogonal planes",
            "step_largest_error": step_err,
            "step_has_a_translation": False,
            "n_steps_return_to_start": math.dist(back, pts[0]) < 1e-9,
            "consecutive_simplices_share_a_facet_of": m,
            "hinge_angle_deg": {"min": min(hinges), "max": max(hinges),
                                "same_at_every_facet": max(hinges) - min(hinges) < 1e-7},
            "flat_helix_same_n": flat,
            "chain_turns_beside_flat_turns": nearest,
            "chain_turns_left_over_deg": left}


def run(seed=0):
    r6 = emlib.load_run(emlib.RUNS / "corkscrew6_seed0.json")["outputs"]["chains"]
    unmatched = {(r["n"], tuple(r["turns"]))
                 for r in emlib.load_run(emlib.RUNS / "coxeter_exponents_seed0.json")
                 ["outputs"]["chains"] if r["verdict"] == "no match"}
    survivors = []
    for c in r6["chains_using_all_three_planes"]["4"]["listed"]:
        if (c["n"], (c["p"], c["q"], c["r"])) in unmatched:
            survivors.append(read(f"6D 5-cells N={c['n']} ({c['p']}, {c['q']}, {c['r']})",
                                  c["n"], (c["p"], c["q"], c["r"]),
                                  c["radius_squared_shares"], 4))
    r4 = emlib.load_run(emlib.RUNS / "corkscrew_seed0.json")["outputs"]["tetrahedra_chain"]["r4"]
    ring_src = next(c for c in r4["closed_chains"]
                    if c["n_tetrahedra_to_close"] == 30 and c["p"] == 1)
    t = ring_src["r1_squared_share"]
    ring = read("4D tetrahedra N=30 (1, 11): the 30-ring", 30, (1, 11), (t, 1 - t), 3)
    keys = ("simplex_m", "ambient_dimension", "points_span_dimension", "sphere_dimension",
            "sphere_dimension_minus_simplex_dimension", "every_window_is_a_regular_simplex",
            "step_is_one_isometry", "step_has_a_translation", "edge_over_radius",
            "consecutive_simplices_share_a_facet_of")
    side = [{"chain": x["chain"], **{k: x[k] for k in keys},
             "hinge_angle_deg": x["hinge_angle_deg"]["min"],
             "flat_turns_deg": x["flat_helix_same_n"]["turns_deg"],
             "chain_turns_deg": x["turn_deg"],
             "differences_from_flat_turns_deg":
                 [d["difference_deg"] for d in x["chain_turns_beside_flat_turns"]],
             "turns_left_over_deg": x["chain_turns_left_over_deg"]}
            for x in [ring] + survivors]
    return {"bound": {"seed": seed, "tol": TOL,
                      "vertices": "rebuilt from the logged turns and radius shares",
                      "flat_turns": "unit roots of the bracket of REFERENCES R5 for n = simplex "
                                    "vertex count, found by scanning and bisection",
                      "hinge_angle": "between the two apexes' heights above the shared facet",
                      "only_these_chains": "the four with no Coxeter match, and the 30-ring"},
            "ring_reference": ring, "survivors": survivors, "side_by_side": side,
            "summary": {"survivors": len(survivors),
                        "every_window_regular": sum(1 for x in survivors
                                                    if x["every_window_is_a_regular_simplex"]),
                        "step_is_one_isometry": sum(1 for x in survivors
                                                    if x["step_is_one_isometry"]),
                        "sphere_minus_simplex_dimension":
                            {x["chain"]: x["sphere_dimension_minus_simplex_dimension"]
                             for x in [ring] + survivors},
                        "ring_every_window_regular": ring["every_window_is_a_regular_simplex"],
                        "ring_step_is_one_isometry": ring["step_is_one_isometry"]}}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    comp = emlib.compressibility(
        [("y" if x["every_window_is_a_regular_simplex"] else "n")
         + ("y" if x["step_is_one_isometry"] else "n") for x in out["survivors"]],
        "per survivor: regular windows y/n, then one isometry y/n")
    angles = {x["chain"]: x["turn_deg"] for x in [out["ring_reference"]] + out["survivors"]}
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME, seed, {"sources": ["runs/corkscrew6_seed0.json",
                                                 "runs/corkscrew_seed0.json",
                                                 "runs/coxeter_exponents_seed0.json"],
                                     "status_before_run": "PROPOSED"},
                        out, angles, None, comp, None)
    print(p)
    for x in out["side_by_side"]:
        print(x["chain"])
        print("   regular windows", x["every_window_is_a_regular_simplex"], "| one isometry",
              x["step_is_one_isometry"], "| translation", x["step_has_a_translation"],
              "| simplex dim", x["simplex_m"], "on S^%d" % x["sphere_dimension"], "(gap",
              x["sphere_dimension_minus_simplex_dimension"], ") | edge/R",
              round(x["edge_over_radius"], 5), "| hinge", round(x["hinge_angle_deg"], 3))
        print("   chain turns", [round(a, 3) for a in x["chain_turns_deg"]], "| flat turns",
              [round(a, 3) for a in x["flat_turns_deg"]], "| differences",
              [round(a, 3) for a in x["differences_from_flat_turns_deg"]], "| left over",
              [round(a, 3) for a in x["turns_left_over_deg"]])
    print("SUMMARY", out["summary"])


if __name__ == "__main__":
    main(sys.argv[1:])
