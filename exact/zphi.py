"""ZPhi -- numbers of the form a + b*phi, exactly. CC0. Co-authored. Stdlib only.

a and b are fractions.Fraction; phi is the golden ratio, entering only
through phi^2 = phi + 1. Nothing here converts to a float: a float handed
in is refused, and there is no __float__. Comparison is decided by exact
integer / rational arithmetic (the sign of p + q*sqrt(5) from p, q and
p^2 against 5*q^2).

  add, sub, neg, mul    closed
  inverse, truediv      closed (the field Q(phi)); division by zero raises
  conjugate             phi -> 1 - phi
  norm                  x * conjugate(x) = a^2 + a*b - b^2, a Fraction
  compare               <, <=, >, >=, ==, sign()
"""
from fractions import Fraction
from numbers import Rational


def _frac(x):
    if isinstance(x, bool) or not isinstance(x, Rational):
        raise TypeError(f"exact arithmetic takes int or Fraction, not {type(x).__name__}")
    return Fraction(x)


class ZPhi:
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        object.__setattr__(self, "a", _frac(a))
        object.__setattr__(self, "b", _frac(b))

    def __setattr__(self, *_):
        raise AttributeError("ZPhi is immutable")

    @staticmethod
    def of(x):
        return x if isinstance(x, ZPhi) else ZPhi(x, 0)

    @staticmethod
    def _other(x):
        """ZPhi for a ZPhi / int / Fraction; None for anything else, so that the
        operator returns NotImplemented and the other operand can answer."""
        if isinstance(x, ZPhi):
            return x
        if isinstance(x, Rational) and not isinstance(x, bool):
            return ZPhi(x, 0)
        if isinstance(x, float):
            raise TypeError("exact arithmetic takes int or Fraction, not float")
        return None

    # ---- ring / field operations
    def __add__(self, o):
        o = ZPhi._other(o)
        if o is None:
            return NotImplemented
        return ZPhi(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return ZPhi(-self.a, -self.b)

    def __sub__(self, o):
        o = ZPhi._other(o)
        if o is None:
            return NotImplemented
        return ZPhi(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        return ZPhi.of(o) - self

    def __mul__(self, o):
        o = ZPhi._other(o)
        if o is None:
            return NotImplemented
        bd = self.b * o.b
        # (a + b phi)(c + d phi) = ac + (ad + bc) phi + bd phi^2,  phi^2 = phi + 1
        return ZPhi(self.a * o.a + bd, self.a * o.b + self.b * o.a + bd)

    __rmul__ = __mul__

    def conjugate(self):
        """phi -> 1 - phi."""
        return ZPhi(self.a + self.b, -self.b)

    def norm(self):
        return self.a * self.a + self.a * self.b - self.b * self.b

    def inverse(self):
        n = self.norm()
        if n == 0:
            raise ZeroDivisionError("ZPhi zero has no inverse")
        c = self.conjugate()
        return ZPhi(c.a / n, c.b / n)

    def __truediv__(self, o):
        return self * ZPhi.of(o).inverse()

    def __rtruediv__(self, o):
        return ZPhi.of(o) * self.inverse()

    def __pow__(self, n):
        if isinstance(n, bool) or not isinstance(n, int):
            raise TypeError("integer powers only")
        base, n = (self.inverse(), -n) if n < 0 else (self, n)
        out = ONE
        while n:
            if n & 1:
                out = out * base
            base = base * base
            n >>= 1
        return out

    # ---- comparison, exact
    def sign(self):
        """-1, 0 or 1: the sign of a + b*phi = p + q*sqrt(5), p = a + b/2, q = b/2."""
        p, q = self.a + self.b / 2, self.b / 2
        sp = (p > 0) - (p < 0)
        sq = (q > 0) - (q < 0)
        if sp == sq or sq == 0:
            return sp
        if sp == 0:
            return sq
        lhs, rhs = p * p, 5 * q * q          # opposite signs: the larger magnitude wins
        if lhs == rhs:
            return 0                         # cannot happen for rational p, q != 0
        return sp if lhs > rhs else sq

    def __eq__(self, o):
        if isinstance(o, ZPhi):
            return self.a == o.a and self.b == o.b
        if isinstance(o, Rational) and not isinstance(o, bool):
            return self.b == 0 and self.a == o
        return NotImplemented

    def __hash__(self):
        return hash((self.a, self.b))

    def __lt__(self, o):
        return (self - o).sign() < 0

    def __le__(self, o):
        return (self - o).sign() <= 0

    def __gt__(self, o):
        return (self - o).sign() > 0

    def __ge__(self, o):
        return (self - o).sign() >= 0

    def __abs__(self):
        return -self if self.sign() < 0 else self

    def __bool__(self):
        return self.a != 0 or self.b != 0

    # ---- no floats
    def __float__(self):
        raise TypeError("ZPhi does not convert to float; use pair() or bounds()")

    def pair(self):
        """(a, b) as [numerator, denominator] pairs, for logging."""
        return [[self.a.numerator, self.a.denominator], [self.b.numerator, self.b.denominator]]

    def bounds(self, k=40):
        """Two Fractions lo <= self <= hi, from consecutive Fibonacci quotients
        around phi (F(k+1)/F(k)). Exact rationals, not floats."""
        f0, f1 = 1, 1
        for _ in range(k):
            f0, f1 = f1, f0 + f1
        x, y = Fraction(f1, f0), Fraction(f0 + f1, f1)   # phi lies between them
        lo, hi = min(x, y), max(x, y)
        if self.b >= 0:
            return self.a + self.b * lo, self.a + self.b * hi
        return self.a + self.b * hi, self.a + self.b * lo

    def __repr__(self):
        return f"ZPhi({self.a}, {self.b})"

    def __str__(self):
        if self.b == 0:
            return str(self.a)
        if self.a == 0:
            return f"{self.b}*phi"
        return f"{self.a} {'+' if self.b > 0 else '-'} {abs(self.b)}*phi"


ZERO = ZPhi(0, 0)
ONE = ZPhi(1, 0)
PHI = ZPhi(0, 1)
