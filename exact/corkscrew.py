"""corkscrew -- a double rotation in R^4, carried exactly.

CC0. Co-authored. Stdlib only. No floats.

A point of R^4 is a pair (z1, z2) of plane points (Cyc5). A double rotation
turns the first plane by one angle and the second by another. Turns by
multiples of 36 degrees are multiplications by powers of Z10, so the
rotation, the return step and the miss are exact.

The regular 5-cell is exact in this form: its vertices are (zeta^k,
zeta^(2k)), k = 0..4.
"""
from .complement import Z10
from .cyc5 import Cyc5, ZETA, norm2
from .zphi import ZERO

ONE = Cyc5(1, 0)


def power(z, n):
    out = ONE
    for _ in range(n % 10 if z == Z10 else n):
        out = out * z
    return out


def five_cell():
    """Vertices (zeta^k, zeta^(2k)), k = 0..4."""
    out, a, b = [], ONE, ONE
    for _ in range(5):
        out.append((a, b))
        a, b = a * ZETA, b * ZETA * ZETA
    return out


def duoprism_5_5():
    """The 25 points (zeta^j, zeta^k)."""
    p = [power(ZETA, j) for j in range(5)]
    return [(a, b) for a in p for b in p]


def dist2(p, q):
    """Squared distance between two points of R^4, a ZPhi."""
    return norm2(p[0] - q[0]) + norm2(p[1] - q[1])


def turn(points, p, q):
    """One step: the first plane turned by p * 36 degrees, the second by q * 36."""
    a, b = power(Z10, p), power(Z10, q)
    return [(z1 * a, z2 * b) for z1, z2 in points]


def closure(points, p, q, max_steps=40):
    """Step the double rotation (p, q) x 36 degrees until every point is back
    where it started. Returns the first step at which the SET is back, the
    first step at which every POINT is back, and the exact squared miss of
    the first point at every step up to that."""
    start = list(points)
    sset = set(start)
    cur = start
    set_step = point_step = None
    misses = []
    for k in range(1, max_steps + 1):
        cur = turn(cur, p, q)
        misses.append(dist2(cur[0], start[0]))
        if set_step is None and set(cur) == sset:
            set_step = k
        if cur == start:
            point_step = k
            break
    return {"turns_of_36_deg": [p, q], "set_return_step": set_step,
            "point_return_step": point_step,
            "miss_squared_of_first_point_per_step": [m.pair() for m in misses],
            "miss_is_exactly_zero_at_return": bool(misses) and misses[-1] == ZERO
            and point_step is not None,
            "arithmetic": "exact a + b*phi forms; equality is ==, no tolerance"}


def is_regular_simplex(points):
    """True exactly when every pair of points is at the same distance."""
    d = {dist2(points[i], points[j]) for i in range(len(points))
         for j in range(i + 1, len(points))}
    return len(d) == 1
