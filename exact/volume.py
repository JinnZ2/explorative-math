"""volume -- exact 4-volumes of polytopes with coordinates in a + b*phi forms.

CC0. Co-authored. Stdlib only. No floats.

sqrt(5) = 2*phi - 1, so anything with rational-and-sqrt(5) coordinates is
exact here. The regular 5-cell inscribed in the unit sphere of R^4 is such
a set.
"""
from fractions import Fraction

from .zphi import ZPhi, PHI, ZERO

SQRT5 = PHI * 2 - 1


def five_cell_unit():
    """The regular 5-cell on the unit sphere: (s, s, s, -1/4), (s, -s, -s, -1/4),
    (-s, s, -s, -1/4), (-s, -s, s, -1/4), (0, 0, 0, 1), with s = sqrt(5)/4."""
    s = SQRT5 * Fraction(1, 4)
    q = ZPhi(Fraction(-1, 4), 0)
    return [(s, s, s, q), (s, -s, -s, q), (-s, s, -s, q), (-s, -s, s, q),
            (ZERO, ZERO, ZERO, ZPhi(1, 0))]


def negated(points):
    return [tuple(-x for x in p) for p in points]


def dot(a, b):
    out = ZERO
    for x, y in zip(a, b):
        out = out + x * y
    return out


def det(rows):
    """Determinant by cofactor expansion (exact)."""
    n = len(rows)
    if n == 1:
        return rows[0][0]
    out = ZERO
    for j in range(n):
        if not rows[0][j]:
            continue
        minor = [list(r[:j]) + list(r[j + 1:]) for r in rows[1:]]
        term = rows[0][j] * det(minor)
        out = out + term if j % 2 == 0 else out - term
    return out


def volume_from_facets(points, facets):
    """4-volume of a polytope that contains the origin, from its tetrahedral
    facets (4 point indices each): the sum of |det(a, b, c, d)| / 24."""
    total = ZERO
    for f in facets:
        if len(f) != 4:
            raise ValueError("tetrahedral facets only")
        d = det([list(points[i]) for i in f])
        total = total + abs(d)
    return total * Fraction(1, 24)
