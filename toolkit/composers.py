"""composers -- ways to combine primitives into a shape or a sequence.

CC0. Co-authored. Stdlib only.

  rotate_project   rotate in chosen planes, then project
  phi_scale        scale a shape by phi^k (optionally overlay the copies)
  cut_and_project  keep lattice points inside a window in the perpendicular
                   space, project them onto the parallel space; optionally
                   weight each kept point by a function of its perpendicular
                   coordinates (gaussian_weight) beside or instead of the
                   hard window
  iterate          iterate a map from a start point, keep the orbit
  overlay          union of two shapes, each point tagged with its source
"""
import math

from . import primitives as P


def rotate_project(points, plane_angles, coords=(0, 1), basis=None):
    return P.project(P.rotate_many(points, plane_angles), basis=basis, coords=coords)


def phi_scale(points, k=1):
    f = P.PHI ** k
    return [P.scale(p, f) for p in points]


def phi_stack(points, ks=(0, 1)):
    """Overlay copies scaled by phi^k for each k."""
    out = []
    for k in ks:
        out.extend(overlay_tagged(phi_scale(points, k), f"phi^{k}"))
    return out


def cut_and_project(lattice, parallel, perpendicular, window=None, weight=None):
    """lattice: points; parallel / perpendicular: orthonormal bases (lists of vectors).

    window: a function taking the perpendicular coordinates and returning True
    if the point is kept (hard window). None keeps every point.
    weight: a function taking the perpendicular coordinates and returning a
    number (weighting-function window). None gives the hard-window result.

    Without weight: returns the kept points' parallel coordinates.
    With weight: returns one {"point", "perp", "weight"} record per kept point;
    nothing is dropped for having a small weight."""
    out = []
    for p in lattice:
        q = [P.dot(p, b) for b in perpendicular]
        if window is not None and not window(q):
            continue
        x = [P.dot(p, b) for b in parallel]
        if weight is None:
            out.append(x)
        else:
            out.append({"point": x, "perp": q, "weight": weight(q)})
    return out


def gaussian_weight(sigma, center=None):
    """Weighting-function window: exp(-d^2 / (2 sigma^2)), d = distance of the
    perpendicular coordinates from `center` (default: the origin, so d is the
    distance from the lattice point to the cut)."""
    if sigma <= 0:
        raise ValueError("sigma must be positive")

    def w(q):
        c = center if center is not None else [0.0] * len(q)
        d2 = sum((a - b) ** 2 for a, b in zip(q, c))
        return math.exp(-d2 / (2.0 * sigma * sigma))
    return w


def strip_bounds_1d(perp_vec, offset=1e-7):
    """(lo, hi) of the hard strip window used by strip_window_1d."""
    corners = [P.dot(c, perp_vec) for c in ((0, 0), (1, 0), (0, 1), (1, 1))]
    return min(corners) + offset, max(corners) + offset


def strip_window_1d(perp_vec, offset=1e-7):
    """Window for a line in Z^2: the projection of the unit square onto perp_vec."""
    lo, hi = strip_bounds_1d(perp_vec, offset)
    return lambda q: lo <= q[0] < hi


def line_basis_2d(slope):
    n = math.hypot(1.0, slope)
    return [[1.0 / n, slope / n]], [[-slope / n, 1.0 / n]]


def iterate(f, x0, n):
    xs = [x0]
    for _ in range(n):
        xs.append(f(xs[-1]))
    return xs


def overlay_tagged(points, tag):
    return [{"tag": tag, "point": list(p)} for p in points]


def overlay(a, b, tags=("A", "B")):
    return overlay_tagged(a, tags[0]) + overlay_tagged(b, tags[1])
