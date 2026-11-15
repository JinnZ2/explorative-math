"""ao. ring_tiling -- can copies of the N = 24 ring of 5-cells tile S^4?

CC0. Co-authored. Stdlib only. Status: PROPOSED.

The N = 24 chain of spherical_helix: 24 regular 4-simplices whose vertices
v_k = (r1 e^(ik a), r2 e^(ik b), r3 (-1)^k) lie on S^4 in R^5. Each simplex
is pushed out to the sphere (the cone over it from the centre); a direction
is covered by a simplex when it lies in that cone.

  ring      how much of S^4 the 24 simplices cover, and whether they
            overlap each other
  cliques   every regular 4-simplex with the ring's edge length among the
            ring's own 24 vertices
  copies    congruent copies of the ring under the rotations its own step
            is built from, taken one plane at a time: plane (2,3) turned by
            j * 15 degrees, and the fifth axis kept or flipped. (Turning
            both planes together by the step only moves the ring along
            itself.) Copies are placed one at a time, in order, and a copy
            is kept only if it overlaps nothing already placed.
  compare   the boundary of the 5-orthoplex: 32 regular 4-simplices that
            tile S^4, read the same way.

Coverage is read on random directions (seeded): the fraction covered 0, 1,
2 or more times. Observer slots null.

  python3 experiments/ring_tiling.py [seed]
  python3 experiments/ring_tiling.py --small [seed]   # phone bound; own run file
"""
import itertools
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from toolkit import primitives as P  # noqa: E402

NAME = "ring_tiling"
N, TURNS = 24, (1, 7, 12)
SAMPLES = 6000
SMALL = {"samples": 1500}
EDGE_TOL = 1e-7


def ring_vertices(shares, phase2=0, flip=False):
    """The 24 vertices in R^5; phase2 turns plane (2,3) by phase2 * 15 degrees,
    flip negates the fifth axis."""
    out = []
    for k in range(N):
        a = 2 * math.pi * TURNS[0] * k / N
        b = 2 * math.pi * (TURNS[1] * k + phase2) / N
        s = math.sqrt(shares[2]) * (-1) ** k * (-1 if flip else 1)
        out.append([math.sqrt(shares[0]) * math.cos(a), math.sqrt(shares[0]) * math.sin(a),
                    math.sqrt(shares[1]) * math.cos(b), math.sqrt(shares[1]) * math.sin(b), s])
    return out


def inverse(m):
    n = len(m)
    a = [list(row) + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(m)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(a[r][i]))
        if abs(a[piv][i]) < 1e-12:
            return None
        a[i], a[piv] = a[piv], a[i]
        d = a[i][i]
        a[i] = [x / d for x in a[i]]
        for r in range(n):
            if r != i and a[r][i]:
                f = a[r][i]
                a[r] = [x - f * y for x, y in zip(a[r], a[i])]
    return [row[n:] for row in a]


def cones(simplices):
    """Per simplex: the matrix taking a direction to its coordinates in the
    simplex's vertices (all >= 0 inside the cone), and its unit centroid."""
    out = []
    for verts in simplices:
        cols = [[verts[j][i] for j in range(5)] for i in range(5)]     # columns = vertices
        inv = inverse(cols)
        c = P.unit([sum(v[i] for v in verts) for i in range(5)])
        reach = min(P.dot(c, P.unit(v)) for v in verts)
        out.append((inv, c, reach))
    return out


def covered_by(x, cone_list):
    """Indices of the cones that contain direction x."""
    hits = []
    for i, (inv, c, reach) in enumerate(cone_list):
        if inv is None or P.dot(x, c) < reach - 1e-9:
            continue
        if all(sum(r * xv for r, xv in zip(row, x)) >= -1e-9 for row in inv):
            hits.append(i)
    return hits


def coverage(samples, cone_list):
    counts = [len(covered_by(x, cone_list)) for x in samples]
    n = len(samples)
    return {"samples": n, "uncovered": sum(1 for c in counts if c == 0) / n,
            "covered_once": sum(1 for c in counts if c == 1) / n,
            "covered_twice_or_more": sum(1 for c in counts if c >= 2) / n,
            "largest_cover_count": max(counts),
            "mean_cover_count": sum(counts) / n}, counts


def dihedral(simplex):
    """Dihedral angle (degrees) between two facets of the spherical simplex over
    these five vertices: from the inverse of the vertices' Gram matrix,
    cos = -g_ij / sqrt(g_ii g_jj). For a regular simplex every pair agrees."""
    u = [P.unit(v) for v in simplex]
    g = inverse([[P.dot(a, b) for b in u] for a in u])
    vals = sorted({round(math.degrees(math.acos(-g[i][j] / math.sqrt(g[i][i] * g[j][j]))), 9)
                   for i in range(5) for j in range(i + 1, 5)})
    return vals


def run(seed=0, n_samples=SAMPLES):
    src = emlib.load_run(emlib.RUNS / "spherical_helix_seed0.json")
    chain = next(x for x in src["outputs"]["survivors"] if x["n_to_close"] == N)
    shares = chain["radius_squared_shares"]
    rng = random.Random(f"{seed}/ring-tiling")
    samples = [P.unit([rng.gauss(0, 1) for _ in range(5)]) for _ in range(n_samples)]

    verts = ring_vertices(shares)
    edge = math.dist(verts[0], verts[1])
    ring = [[verts[(k + j) % N] for j in range(5)] for k in range(N)]
    ring_cones = cones(ring)
    ring_cov, ring_counts = coverage(samples, ring_cones)
    one = coverage(samples, ring_cones[:1])[0]

    near = [[abs(math.dist(verts[i], verts[j]) - edge) < EDGE_TOL for j in range(N)]
            for i in range(N)]
    cliques = [c for c in itertools.combinations(range(N), 5)
               if all(near[a][b] for a, b in itertools.combinations(c, 2))]
    ring_sets = {tuple(sorted((k + j) % N for j in range(5))) for k in range(N)}
    clique_cov = coverage(samples, cones([[verts[i] for i in c] for c in cliques]))[0]

    placed = [{"phase2": 0, "flip": False}]
    counts = list(ring_counts)
    tried, rejected = 0, []
    for flip in (False, True):
        for j in range(N):
            if j == 0 and not flip:
                continue
            tried += 1
            v = ring_vertices(shares, j, flip)
            cs = cones([[v[(k + i) % N] for i in range(5)] for k in range(N)])
            own = [len(covered_by(x, cs)) for x in samples]
            clash = sum(1 for a, b in zip(counts, own) if a and b)
            if clash == 0 and any(own):
                placed.append({"phase2": j, "flip": flip})
                counts = [a + b for a, b in zip(counts, own)]
            else:
                rejected.append({"phase2": j, "flip": flip, "overlapping_samples": clash,
                                 "overlap_fraction": clash / n_samples})
    n = n_samples
    tiling = {"copies_tried": tried, "copies_placed": len(placed),
              "simplices_placed": len(placed) * N, "placed": placed,
              "rejected": rejected,
              "uncovered_fraction_gaps": sum(1 for c in counts if c == 0) / n,
              "covered_once_fraction": sum(1 for c in counts if c == 1) / n,
              "overlap_fraction": sum(1 for c in counts if c >= 2) / n,
              "tiles_the_sphere": all(c == 1 for c in counts)}

    orth = P.cross_polytope(5)
    facets = [[orth[2 * i + (0 if sg[i] > 0 else 1)] for i in range(5)]
              for sg in itertools.product((1, -1), repeat=5)]
    orth_cov = coverage(samples, cones(facets))[0]
    orth_one = coverage(samples, cones(facets[:1]))[0]

    f1 = one["covered_once"]
    d_ring, d_orth = dihedral(ring[0]), dihedral(facets[0])

    def around(d):
        return {"dihedral_deg": d, "one_value": len(d) == 1,
                "simplices_round_a_shared_triangle": 360.0 / d[0],
                "whole_number": abs(360.0 / d[0] - round(360.0 / d[0])) < 1e-6}
    return {"bound": {"seed": seed, "samples": n_samples,
                      "coverage": "fraction of random directions on S^4 inside the cone over "
                                  "a simplex; a statistical reading, not exact",
                      "one_sample_is": 1.0 / n_samples,
                      "copies": "plane (2,3) turned by j * 15 degrees, j = 0..23, fifth axis "
                                "kept or flipped: 47 candidates beside the ring; placed in "
                                "that order, a copy kept only if no sample is covered twice",
                      "not_tried": "other rotations; other placement orders; an exact test"},
            "ring": {"n_simplices": N, "edge_over_radius": edge,
                     "radius_squared_shares": shares,
                     "fraction_of_sphere_per_simplex": f1,
                     "simplices_needed_at_that_fraction": (1.0 / f1) if f1 else None,
                     "coverage": ring_cov,
                     "ring_overlaps_itself": ring_cov["covered_twice_or_more"] > 0},
            "regular_simplices_among_the_ring_vertices": {
                "n": len(cliques), "all_are_ring_simplices": set(cliques) == ring_sets,
                "n_not_in_the_ring": len(set(cliques) - ring_sets),
                "coverage": clique_cov},
            "copies": tiling,
            "orthoplex_5": {"n_facets": len(facets), "edge_over_radius": math.sqrt(2.0),
                            "fraction_of_sphere_per_facet": orth_one["covered_once"],
                            "exact_fraction_per_facet": 1.0 / 32, "coverage": orth_cov,
                            "tiles_the_sphere": orth_cov["covered_once"] == 1.0},
            "dihedral": {"ring_simplex": around(d_ring), "orthoplex_facet": around(d_orth),
                         "derived": "congruent simplices meeting face to face round a "
                                    "shared triangle close up only if 360 / dihedral is a "
                                    "whole number"},
            "side_by_side": {"simplices": [len(placed) * N, 32],
                             "dihedral_deg": [d_ring[0], d_orth[0]],
                             "simplices_round_a_triangle": [360.0 / d_ring[0], 360.0 / d_orth[0]],
                             "edge_over_radius": [edge, math.sqrt(2.0)],
                             "fraction_per_simplex": [f1, orth_one["covered_once"]],
                             "covered_once": [tiling["covered_once_fraction"],
                                              orth_cov["covered_once"]],
                             "gaps": [tiling["uncovered_fraction_gaps"], orth_cov["uncovered"]],
                             "overlaps": [tiling["overlap_fraction"],
                                          orth_cov["covered_twice_or_more"]],
                             "order": ["N = 24 ring and its placed copies", "5-orthoplex boundary"]}}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = run(seed, SMALL["samples"] if small else SAMPLES)
    comp = emlib.compressibility(
        ["p"] + ["r" for _ in out["copies"]["rejected"]]
        + ["p" for _ in out["copies"]["placed"][1:]],
        "the ring (p), then one letter per candidate copy: r rejected, p placed")
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "samples": out["bound"]["samples"], "n": N,
                         "turns": list(TURNS), "source": "runs/spherical_helix_seed0.json",
                         "status_before_run": "PROPOSED"},
                        out, None, None, comp, None)
    print(p)
    r, c, o = out["ring"], out["copies"], out["orthoplex_5"]
    print("RING per simplex", r["fraction_of_sphere_per_simplex"], "needed",
          r["simplices_needed_at_that_fraction"], "| coverage", r["coverage"])
    print("CLIQUES", out["regular_simplices_among_the_ring_vertices"])
    print("COPIES tried", c["copies_tried"], "placed", c["copies_placed"], c["placed"],
          "| gaps", c["uncovered_fraction_gaps"], "once", c["covered_once_fraction"],
          "overlap", c["overlap_fraction"])
    print("   rejected overlap fractions",
          sorted({round(x["overlap_fraction"], 4) for x in c["rejected"]})[:12])
    print("DIHEDRAL", out["dihedral"]["ring_simplex"], out["dihedral"]["orthoplex_facet"])
    print("ORTHOPLEX facets", o["n_facets"], "per facet", o["fraction_of_sphere_per_facet"],
          "coverage", o["coverage"])


if __name__ == "__main__":
    main(sys.argv[1:])
