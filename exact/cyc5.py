"""Cyc5 -- points of the plane built on fifth roots of unity, exactly.

CC0. Co-authored. Stdlib only. No floats.

A point is u + v*w with u, v in ZPhi, where w = 2i*sin(72 degrees) and so
w^2 = -(2 + phi). The real part is u; the imaginary part is v times the
fixed irrational |w|, which is never evaluated. The fifth root of unity is

  zeta = cos 72 + i sin 72 = (phi - 1)/2 + w/2

so every vertex of a regular pentagon, a pentagram or a Penrose tiling
drawn from fifth roots is a Cyc5, and every squared length, dot product
and (cross product / |w|) between such points is a ZPhi.
"""
from fractions import Fraction

from .zphi import ZPhi, PHI, ONE, ZERO

_W2 = -(ZPhi(2, 0) + PHI)        # w^2
_NEG_W2 = ZPhi(2, 0) + PHI       # |w|^2 = 2 + phi


class Cyc5:
    __slots__ = ("u", "v")

    def __init__(self, u=0, v=0):
        object.__setattr__(self, "u", ZPhi.of(u))
        object.__setattr__(self, "v", ZPhi.of(v))

    def __setattr__(self, *_):
        raise AttributeError("Cyc5 is immutable")

    @staticmethod
    def of(x):
        return x if isinstance(x, Cyc5) else Cyc5(x, 0)

    def __add__(self, o):
        o = Cyc5.of(o)
        return Cyc5(self.u + o.u, self.v + o.v)

    __radd__ = __add__

    def __neg__(self):
        return Cyc5(-self.u, -self.v)

    def __sub__(self, o):
        o = Cyc5.of(o)
        return Cyc5(self.u - o.u, self.v - o.v)

    def __mul__(self, o):
        if not isinstance(o, Cyc5):
            o = ZPhi.of(o)
            return Cyc5(self.u * o, self.v * o)          # scaling by a real ZPhi
        return Cyc5(self.u * o.u + _W2 * (self.v * o.v), self.u * o.v + self.v * o.u)

    __rmul__ = __mul__

    def conjugate(self):
        """Complex conjugate: w -> -w."""
        return Cyc5(self.u, -self.v)

    def __eq__(self, o):
        if isinstance(o, Cyc5):
            return self.u == o.u and self.v == o.v
        return NotImplemented

    def __hash__(self):
        return hash((self.u, self.v))

    def __float__(self):
        raise TypeError("Cyc5 does not convert to float")

    def pair(self):
        return [self.u.pair(), self.v.pair()]

    def __repr__(self):
        return f"Cyc5({self.u!r}, {self.v!r})"


W = Cyc5(ZERO, ONE)
ZETA = Cyc5((PHI - 1) * Fraction(1, 2), ZPhi(Fraction(1, 2), 0))


def norm2(z):
    """|z|^2, a ZPhi."""
    return z.u * z.u + _NEG_W2 * (z.v * z.v)


def dot(a, b):
    """Real dot product of two points, a ZPhi."""
    return a.u * b.u + _NEG_W2 * (a.v * b.v)


def cross(a, b):
    """(a x b) / |w|, a ZPhi: 0 exactly when a and b are parallel, and its sign
    is the orientation. The true cross product is this times the fixed |w|."""
    return a.u * b.v - a.v * b.u
