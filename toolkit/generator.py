"""generator -- a random combination of composers becomes a new experiment.

CC0. Co-authored. Stdlib only.

From a seed, pick one of three pipelines and its parameters, run it, apply
the probes, and log the run through loghook (interest slots void). The spec
that was drawn is stored in `params`, so a run can be re-made from its file.

  shape     source shape -> optional phi_stack -> rotate in random planes
            -> project to 2D -> symmetry, mirror, polar-angle histogram, zlib
  cut       Z^2 cut at a random slope -> projected chain -> gap lengths,
            period profile, zlib
  orbit     iterate x -> x + alpha (mod 1) for a random alpha -> symbol
            sequence (x < 1/2 ?) -> period profile, zlib

  cut_weighted   NOT drawn from the seed: asked for by name. Z^2 cut at a
            random slope with a Gaussian weight by distance to the cut
            (random sigma) -> weighted points, read at several weight levels
            -> gap lengths per level. Written to its own run files, so the
            runs drawn from the three pipelines above do not move.

Nothing here decides whether a run is interesting. Every generated run is
written, including the ones that show nothing.

  python3 -m toolkit.generator SEED [SEED ...]
  python3 -m toolkit.generator --pipeline cut_weighted SEED [SEED ...]
"""
import math
import random
import sys

from . import composers as C
from . import loghook
from . import primitives as P
from . import probes as Q

NAME = "generated"
PIPELINES = ("shape", "cut", "orbit")      # drawn from the seed; order is pinned
OPTION_PIPELINES = ("cut_weighted",)        # only by name; never drawn
WEIGHT_FLOOR = 1e-6                         # recording limit for weighted points
WEIGHT_LEVELS = (0.9, math.exp(-0.5), 0.5, 0.1, 0.01)
SHAPE_SOURCES = ("simplex", "hypercube", "cross_polytope", "cell24", "icosahedron")


def _round(points, nd=6):
    return [[round(x, nd) for x in p] for p in points]


def draw_spec(seed, pipeline=None):
    rng = random.Random(seed)
    if pipeline == "cut_weighted":
        return {"pipeline": pipeline, "slope": rng.uniform(0.05, 1.5), "extent": 60,
                "sigma_over_strip_width": rng.uniform(0.25, 1.0),
                "weight_floor": WEIGHT_FLOOR, "levels": list(WEIGHT_LEVELS)}
    if pipeline is not None and pipeline not in PIPELINES:
        raise ValueError(f"unknown pipeline {pipeline!r}")
    kind = rng.choice(PIPELINES)
    if pipeline is not None and pipeline != kind:
        raise ValueError(f"seed {seed} draws {kind!r}, not {pipeline!r}; "
                         "drawn pipelines cannot be forced")
    if kind == "shape":
        src = rng.choice(SHAPE_SOURCES)
        if src in P.POLYTOPES:
            dim = rng.randint(3, 6)
        else:
            dim = {"cell24": 4, "icosahedron": 3}[src]
        n_rot = rng.randint(0, len(P.all_planes(dim)))
        planes = rng.sample(P.all_planes(dim), n_rot)
        rots = [[list(pl), rng.uniform(0, 2 * math.pi)] for pl in planes]
        stack = rng.choice([[0], [0, 1], [0, -1]])
        return {"pipeline": kind, "source": src, "dim": dim, "rotations": rots,
                "phi_stack": stack}
    if kind == "cut":
        return {"pipeline": kind, "slope": rng.uniform(0.05, 1.5), "extent": 120}
    return {"pipeline": kind, "alpha": rng.random(), "x0": rng.random(), "steps": 512}


def _shape(spec):
    src = spec["source"]
    pts = P.POLYTOPES[src](spec["dim"]) if src in P.POLYTOPES else P.FIXED_SHAPES[src]()
    stacked = [r["point"] for r in C.phi_stack(pts, spec["phi_stack"])]
    proj = C.rotate_project(stacked, [(tuple(pl), th) for pl, th in spec["rotations"]])
    proj = _round(proj, 9)
    ang = Q.polar_angles(proj)
    outputs = {
        "n_points": len(proj), "projected": proj,
        "rotational_symmetry": Q.rotational_symmetry(proj),
        "mirror_axes_deg": Q.mirror_symmetry(proj),
        "polar_angle_histogram": Q.angle_histogram(ang, width=15.0, hi=360.0),
    }
    symbols = [chr(ord("a") + int(a // 15)) for a in ang]
    return outputs, ang, symbols, "polar angles of projected points in 15-degree bins (a..x)"


def _gap_classes(xs):
    """Sorted positions -> (gaps, mean length of each gap class, letters A..)."""
    gaps = [b - a for a, b in zip(xs, xs[1:])]
    lengths = []
    for g in sorted(gaps):
        if not lengths or g - lengths[-1][-1] > 1e-6:
            lengths.append([g])
        else:
            lengths[-1].append(g)
    means = [sum(c) / len(c) for c in lengths]
    letters = "".join(chr(ord("A") + min(range(len(means)), key=lambda i: abs(g - means[i])))
                      for g in gaps)
    return gaps, means, letters


def _cut(spec):
    par, perp = C.line_basis_2d(spec["slope"])
    pts = C.cut_and_project(P.lattice_zn(2, spec["extent"]), par, perp,
                            C.strip_window_1d(perp[0]))
    xs = sorted(p[0] for p in pts)
    margin = spec["extent"] * 0.5 * par[0][0]
    xs = [x for x in xs if -margin <= x <= margin]
    gaps, means, letters = _gap_classes(xs)
    outputs = {"n_gaps": len(gaps), "distinct_gap_lengths": means,
               "word": letters, "period_profile": Q.period_profile(letters),
               "zlib_ratio": Q.zlib_ratio(letters)}
    return outputs, None, list(letters), "gap classes as letters A.. (shortest first)"


def _cut_weighted(spec):
    par, perp = C.line_basis_2d(spec["slope"])
    lo, hi = C.strip_bounds_1d(perp[0])
    sigma = spec["sigma_over_strip_width"] * (hi - lo)
    recs = C.cut_and_project(P.lattice_zn(2, spec["extent"]), par, perp,
                             weight=C.gaussian_weight(sigma))
    margin = spec["extent"] * 0.5 * par[0][0]
    inside = [r for r in recs if -margin <= r["point"][0] <= margin]
    total = sum(r["weight"] for r in inside)
    in_strip = sum(r["weight"] for r in inside if lo <= r["perp"][0] < hi)
    kept = sorted((r for r in inside if r["weight"] >= spec["weight_floor"]),
                  key=lambda r: r["point"][0])
    levels = []
    for t in spec["levels"]:
        xs = [r["point"][0] for r in kept if r["weight"] >= t]
        gaps, means, letters = _gap_classes(xs)
        levels.append({"level": t,
                       "half_width_in_q": sigma * math.sqrt(2.0 * math.log(1.0 / t)),
                       "n_points": len(xs), "distinct_gap_lengths": means,
                       "word": letters, "period_profile": Q.period_profile(letters),
                       "zlib_ratio": Q.zlib_ratio(letters)})
    outputs = {"sigma": sigma, "hard_strip": [lo, hi], "margin": margin,
               "n_lattice_points_in_margin": len(inside),
               "weight_sum_all": total,
               "weight_share_inside_hard_strip": in_strip / total,
               "n_points_above_floor": len(kept),
               "points_x_q_weight": [[r["point"][0], r["perp"][0], r["weight"]]
                                     for r in kept],
               "levels": levels}
    word = levels[1]["word"]
    return (outputs, None, list(word),
            "gap classes as letters A.. (shortest first) of the points with weight "
            ">= e^-1/2 (|q| <= sigma)")


def _orbit(spec):
    a = spec["alpha"]
    xs = C.iterate(lambda x: (x + a) % 1.0, spec["x0"], spec["steps"])
    word = "".join("0" if x < 0.5 else "1" for x in xs)
    outputs = {"orbit_head": [round(x, 9) for x in xs[:32]], "word": word,
               "period_profile": Q.period_profile(word), "zlib_ratio": Q.zlib_ratio(word)}
    return outputs, None, list(word), "orbit coded 0 if x < 1/2 else 1"


def run(seed, pipeline=None):
    spec = draw_spec(seed, pipeline)
    outputs, angles, symbols, definition = {"shape": _shape, "cut": _cut,
                                            "orbit": _orbit,
                                            "cut_weighted": _cut_weighted}[spec["pipeline"]](spec)
    return spec, outputs, angles, symbols, definition


def generate(seed, out_dir=None, pipeline=None):
    spec, outputs, angles, symbols, definition = run(seed, pipeline)
    name = NAME if pipeline in (None,) + PIPELINES else f"{NAME}_{pipeline}"
    return loghook.log_run(name, seed, {"spec": spec}, outputs, angles, symbols,
                           definition, out_dir=out_dir)


if __name__ == "__main__":
    args = sys.argv[1:]
    pipe = None
    if args[:1] == ["--pipeline"]:
        if len(args) < 2 or args[1] not in OPTION_PIPELINES:
            sys.exit(f"usage: --pipeline {{{', '.join(OPTION_PIPELINES)}}} SEED [SEED ...]")
        pipe, args = args[1], args[2:]
    for s in args or ["0"]:
        print(generate(int(s), pipeline=pipe))
