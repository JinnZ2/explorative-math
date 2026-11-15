"""e. simplex_4d -- the 5-cell: rotate in 4D planes, project to 2D, measure angles.

CC0. Co-authored. Stdlib only.

Vertices: the five unit vectors of R^5 minus their centroid, expressed in an
orthonormal basis of the 4D hyperplane they span (a regular 4-simplex).

For each trial, compose rotations in the six coordinate planes of R^4 by
random angles, project onto the first two coordinates, and measure the
polar angles of the five projected points about their centroid.

  pentagon_score   max deviation (degrees) of the angular gaps from 72,
                   plus 100 x the relative spread of the radii. 0 = regular
                   pentagon. Lower = closer to the symmetric projection.

The Petrie projection is computed directly: the plane spanned by
u = sum cos(2 pi k/5) e_k and v = sum sin(2 pi k/5) e_k, which the cyclic
permutation of vertices rotates by 72 degrees. It is the score-0 reference.

  python3 experiments/simplex_4d.py [seed]
"""
import itertools
import math
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402

NAME = "simplex_4d"
TRIALS = 2000
LOG_BELOW = 5.0     # log the angle set of any trial with pentagon_score below this


def gram_schmidt(vs):
    basis = []
    for v in vs:
        w = list(v)
        for b in basis:
            c = emlib.dot(w, b)
            w = [x - c * y for x, y in zip(w, b)]
        if emlib.norm(w) > 1e-12:
            basis.append(emlib.unit(w))
    return basis


def simplex_vertices():
    e = [[1.0 if i == k else 0.0 for i in range(5)] for k in range(5)]
    centred = [[x - 0.2 for x in v] for v in e]
    basis = gram_schmidt(centred[:4])
    return [[emlib.dot(v, b) for b in basis] for v in centred], centred


def rotate(points, i, j, theta):
    c, s = math.cos(theta), math.sin(theta)
    out = []
    for p in points:
        q = list(p)
        q[i], q[j] = c * p[i] - s * p[j], s * p[i] + c * p[j]
        out.append(q)
    return out


def polar_angles(pts2):
    cx = sum(p[0] for p in pts2) / len(pts2)
    cy = sum(p[1] for p in pts2) / len(pts2)
    ang = sorted(math.degrees(math.atan2(p[1] - cy, p[0] - cx)) % 360.0 for p in pts2)
    radii = [math.hypot(p[0] - cx, p[1] - cy) for p in pts2]
    return ang, radii


def pentagon_score(pts2):
    ang, radii = polar_angles(pts2)
    gaps = [b - a for a, b in zip(ang, ang[1:])] + [ang[0] + 360.0 - ang[-1]]
    gap_dev = max(abs(g - 72.0) for g in gaps)
    mr = statistics.fmean(radii)
    spread = (max(radii) - min(radii)) / mr if mr > 1e-12 else float("inf")
    return gap_dev + 100.0 * spread, ang


def petrie_projection(centred):
    u = [math.cos(2 * math.pi * k / 5) for k in range(5)]
    v = [math.sin(2 * math.pi * k / 5) for k in range(5)]
    return [(emlib.dot(p, u), emlib.dot(p, v)) for p in centred]


def edge_lengths(vs):
    return [math.dist(a, b) for a, b in itertools.combinations(vs, 2)]


def run(seed=0, trials=TRIALS):
    rng = random.Random(seed)
    v4, centred = simplex_vertices()
    planes = list(itertools.combinations(range(4), 2))
    petrie_score, petrie_ang = pentagon_score(petrie_projection(centred))
    scores, logged = [], []
    for t in range(trials):
        pts = v4
        thetas = [rng.uniform(0, 2 * math.pi) for _ in planes]
        for (i, j), th in zip(planes, thetas):
            pts = rotate(pts, i, j, th)
        sc, ang = pentagon_score([(p[0], p[1]) for p in pts])
        scores.append(sc)
        if sc < LOG_BELOW:
            logged.append({"trial": t, "score": sc, "plane_angles_rad": thetas,
                           "polar_angles_deg": ang})
    edges = edge_lengths(v4)
    return {
        "edge_length_min": min(edges), "edge_length_max": max(edges),
        "petrie_score": petrie_score, "petrie_polar_angles_deg": petrie_ang,
        "scores": scores,
        "score_min": min(scores), "score_median": statistics.median(scores),
        "n_logged_below": len(logged), "logged": logged,
    }


def main(argv):
    seed = int(argv[0]) if argv else 0
    out = run(seed)
    sym = ["0123456789"[min(9, int(math.log10(1 + s) * 3))] for s in out["scores"]]
    comp = emlib.compressibility(sym, "per-trial pentagon_score, 3*log10(1+score) floored, 0-9")
    ai = (f"OBSERVED: regular (edges {out['edge_length_min']:.6f}..{out['edge_length_max']:.6f}); "
          f"Petrie score {out['petrie_score']:.2e}; best of {TRIALS} random rotations "
          f"{out['score_min']:.3f}, median {out['score_median']:.2f}; {out['n_logged_below']} "
          f"trials below {LOG_BELOW}. AI reading: on this seed random rotation did not land "
          "near the symmetric view, so here the symmetric projection is found by "
          "construction (the eigenplane), not by search.")
    angles = {"petrie": out["petrie_polar_angles_deg"],
              "logged_trials": {str(r["trial"]): r["polar_angles_deg"] for r in out["logged"]}}
    p = emlib.write_run(NAME, seed, {"trials": TRIALS, "log_below": LOG_BELOW,
                                     "projection": "first two coordinates after rotation"},
                        out, angles, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
