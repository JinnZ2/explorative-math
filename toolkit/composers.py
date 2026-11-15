"""composers -- ways to combine primitives into a shape or a sequence.

CC0. Co-authored. Stdlib only.

  rotate_project   rotate in chosen planes, then project
  phi_scale        scale a shape by phi^k (optionally overlay the copies)
  cut_and_project  keep lattice points inside a window in the perpendicular
                   space, project them onto the parallel space
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


def cut_and_project(lattice, parallel, perpendicular, window):
    """lattice: points; parallel / perpendicular: orthonormal bases (lists of vectors).

    window: a function taking the perpendicular coordinates and returning True
    if the point is kept. Returns the kept points' parallel coordinates."""
    out = []
    for p in lattice:
        q = [P.dot(p, b) for b in perpendicular]
        if window(q):
            out.append([P.dot(p, b) for b in parallel])
    return out


def strip_window_1d(perp_vec, offset=1e-7):
    """Window for a line in Z^2: the projection of the unit square onto perp_vec."""
    corners = [P.dot(c, perp_vec) for c in ((0, 0), (1, 0), (0, 1), (1, 1))]
    lo, hi = min(corners) + offset, max(corners) + offset
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
