"""cyclo -- the field Q(zeta_N), exactly. CC0. Co-authored. Stdlib only. No floats
in the arithmetic.

An element is a polynomial in zeta = exp(2 pi i / N) with Fraction
coefficients, reduced modulo the N-th cyclotomic polynomial (degree phi(N)).

  K = Field(18)
  K.cos(k)        cos(2 pi k / N)                 real, exact
  K.sin_i(k)      i sin(2 pi k / N)               so sin a sin b = -sin_i(a) sin_i(b)
  x.galois(u)     zeta -> zeta^u, gcd(u, N) = 1
  x.sign()        -1 / 0 / +1 of a real element, exactly: zero by equality,
                  otherwise a rational interval around the value, narrowed
                  until it excludes 0
  x.approx(d)     a Fraction within 10^-d of a real element's value
  K.coordinates(x, gens)   Fractions q with x = sum q_i gens_i, or None when x
                           is not in the rational span of gens (exact)

No float anywhere: values are read as Fractions (approx) and turned into
floats by the caller. The intervals come from Machin's formula for pi and the
cosine series, both with rational bounds.
"""
from fractions import Fraction


def gcd(a, b):
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a


def _atan_inv(q, terms):
    """(low, high) for atan(1 / q), q > 1: an alternating series, so two
    consecutive partial sums bracket it."""
    x = Fraction(1, q)
    total, power, sums = Fraction(0), x, []
    for k in range(terms + 1):
        total += (-1) ** k * power / (2 * k + 1)
        power *= x * x
        sums.append(total)
    return min(sums[-1], sums[-2]), max(sums[-1], sums[-2])


def pi_bounds(terms):
    """(low, high) for pi: 16 atan(1/5) - 4 atan(1/239)."""
    a_lo, a_hi = _atan_inv(5, terms)
    b_lo, b_hi = _atan_inv(239, terms)
    return 16 * a_lo - 4 * b_hi, 16 * a_hi - 4 * b_lo


def cos_bounds(k, n, terms, _memo={}):
    """(low, high) for cos(2 pi k / N), rational."""
    k %= n
    k = min(k, n - k)                                       # angle in [0, pi]
    key = (k, n, terms)
    if key not in _memo:
        p_lo, p_hi = pi_bounds(terms)
        lo, hi = 2 * p_lo * k / n, 2 * p_hi * k / n
        mid, width = (lo + hi) / 2, hi - lo
        total, term = Fraction(0), Fraction(1)
        for j in range(terms + 4):                          # 4 > pi: the terms shrink from here
            total += term
            term = -term * mid * mid / ((2 * j + 1) * (2 * j + 2))
        slack = abs(term) + width / 2                       # series remainder + |cos'| <= 1
        _memo[key] = (max(Fraction(-1), total - slack), min(Fraction(1), total + slack))
    return _memo[key]


def cyclotomic(n, _memo={}):
    """Coefficients of the n-th cyclotomic polynomial, lowest degree first."""
    if n in _memo:
        return _memo[n]
    num = [-1] + [0] * (n - 1) + [1]                    # x^n - 1
    for d in range(1, n):
        if n % d == 0:
            den = cyclotomic(d)
            q = [0] * (len(num) - len(den) + 1)
            num = num[:]
            for i in range(len(q) - 1, -1, -1):
                q[i] = num[i + len(den) - 1] // den[-1]
                for j, c in enumerate(den):
                    num[i + j] -= q[i] * c
            num = q
    _memo[n] = num
    return num


def solve_exact(rows, rhs, zero, one):
    """Gauss-Jordan over any exact field: rows (k x k), rhs (k). Returns the
    solution list, or None when the matrix is singular. Zero is decided by ==."""
    k = len(rows)
    m = [list(r) + [rhs[i]] for i, r in enumerate(rows)]
    for col in range(k):
        piv = next((r for r in range(col, k) if m[r][col] != zero), None)
        if piv is None:
            return None
        m[col], m[piv] = m[piv], m[col]
        inv = one / m[col][col]
        m[col] = [x * inv for x in m[col]]
        for r in range(k):
            if r != col and m[r][col] != zero:
                f = m[r][col]
                m[r] = [x - f * y for x, y in zip(m[r], m[col])]
    return [m[i][k] for i in range(k)]


class Field:
    def __init__(self, n):
        if n < 1:
            raise ValueError("N >= 1")
        self.n = n
        self.phi = cyclotomic(n)
        self.deg = len(self.phi) - 1
        self.zero = Elt(self, (Fraction(0),) * self.deg)
        self.one = self.of(1)

    def _reduce(self, coeffs):
        c = [Fraction(x) for x in coeffs]
        d, phi = self.deg, self.phi
        for i in range(len(c) - 1, d - 1, -1):
            f = c[i]
            if f:
                for j in range(d + 1):
                    c[i - d + j] -= f * phi[j]
        c = c[:d] + [Fraction(0)] * (d - len(c))
        return Elt(self, tuple(c))

    def of(self, x):
        if isinstance(x, Elt):
            if x.k is not self and x.k.n != self.n:
                raise ValueError("element of another field")
            return x
        return self._reduce([Fraction(x)])

    def zeta(self, k=1):
        c = [Fraction(0)] * self.n
        c[k % self.n] = Fraction(1)
        return self._reduce(c)

    def cos(self, k):
        return (self.zeta(k) + self.zeta(-k)) * Fraction(1, 2)

    def sin_i(self, k):
        """i sin(2 pi k / N)."""
        return (self.zeta(k) - self.zeta(-k)) * Fraction(1, 2)

    def units(self):
        return [u for u in range(1, self.n) if gcd(u, self.n) == 1]

    def coordinates(self, x, gens):
        """Fractions q with x = sum q_i gens_i, or None. Exact."""
        x = self.of(x)
        cols = [self.of(g).c for g in gens]
        rows = [[col[r] for col in cols] + [x.c[r]] for r in range(self.deg)]
        k, piv_cols, r0 = len(cols), [], 0
        for col in range(k):
            piv = next((r for r in range(r0, len(rows)) if rows[r][col]), None)
            if piv is None:
                continue
            rows[r0], rows[piv] = rows[piv], rows[r0]
            inv = 1 / rows[r0][col]
            rows[r0] = [v * inv for v in rows[r0]]
            for r in range(len(rows)):
                if r != r0 and rows[r][col]:
                    f = rows[r][col]
                    rows[r] = [a - f * b for a, b in zip(rows[r], rows[r0])]
            piv_cols.append(col)
            r0 += 1
        if any(row[k] for row in rows[r0:]):
            return None                                   # not in the span
        if len(piv_cols) < k:
            return None                                   # gens are dependent: no unique answer
        return [rows[i][k] for i in range(k)]


class Elt:
    __slots__ = ("k", "c")

    def __init__(self, k, c):
        self.k, self.c = k, c

    def _o(self, o):
        return self.k.of(o)

    def __add__(self, o):
        o = self._o(o)
        return Elt(self.k, tuple(a + b for a, b in zip(self.c, o.c)))

    __radd__ = __add__

    def __neg__(self):
        return Elt(self.k, tuple(-a for a in self.c))

    def __sub__(self, o):
        return self + (-self._o(o))

    def __rsub__(self, o):
        return self._o(o) - self

    def __mul__(self, o):
        if isinstance(o, (int, Fraction)):
            return Elt(self.k, tuple(a * o for a in self.c))
        o = self._o(o)
        out = [Fraction(0)] * (2 * self.k.deg - 1)
        for i, a in enumerate(self.c):
            if a:
                for j, b in enumerate(o.c):
                    if b:
                        out[i + j] += a * b
        return self.k._reduce(out)

    __rmul__ = __mul__

    def inverse(self):
        """Solve self * y = 1 as a linear system over Q."""
        d, k = self.k.deg, self.k
        cols = [(self * k.zeta(j)).c for j in range(d)]
        rows = [[cols[j][r] for j in range(d)] for r in range(d)]
        sol = solve_exact(rows, list(k.one.c), Fraction(0), Fraction(1))
        if sol is None:
            raise ZeroDivisionError("zero has no inverse")
        return Elt(k, tuple(sol))

    def __truediv__(self, o):
        if isinstance(o, (int, Fraction)):
            return self * (Fraction(1) / o)
        return self * self._o(o).inverse()

    def __rtruediv__(self, o):
        return self._o(o) * self.inverse()

    def __eq__(self, o):
        if isinstance(o, (int, Fraction)):
            o = self.k.of(o)
        return isinstance(o, Elt) and o.k.n == self.k.n and self.c == o.c

    def __hash__(self):
        return hash((self.k.n, self.c))

    def galois(self, u):
        if gcd(u, self.k.n) != 1:
            raise ValueError("u is not a unit mod N")
        out = [Fraction(0)] * self.k.n
        for i, a in enumerate(self.c):
            if a:
                out[(i * u) % self.k.n] += a
        return self.k._reduce(out)

    def is_real(self):
        return self.galois(self.k.n - 1) == self if self.k.n > 2 else True

    def bounds(self, terms):
        """(low, high) for the value of a real element."""
        lo = hi = Fraction(0)
        for i, a in enumerate(self.c):
            if a:
                c_lo, c_hi = cos_bounds(i, self.k.n, terms)     # real part of zeta^i
                lo += min(a * c_lo, a * c_hi)
                hi += max(a * c_lo, a * c_hi)
        return lo, hi

    def _real(self):
        if not self.is_real():
            raise ValueError("not a real element")

    def sign(self):
        if self == self.k.zero:
            return 0
        self._real()
        terms = 8
        while True:
            lo, hi = self.bounds(terms)
            if lo > 0:
                return 1
            if hi < 0:
                return -1
            terms *= 2

    def approx(self, digits=18):
        """A Fraction within 10^-digits of the value of a real element."""
        self._real()
        terms, eps = 8, Fraction(1, 10 ** digits)
        while True:
            lo, hi = self.bounds(terms)
            if hi - lo < 2 * eps:
                return (lo + hi) / 2
            terms *= 2

    def __repr__(self):
        return "Elt(%d, %s)" % (self.k.n, [str(a) for a in self.c])


def show(q, names):
    """A rational combination as text: q_i names_i, zeros dropped."""
    parts = []
    for a, nm in zip(q, names):
        if a:
            parts.append(str(a) if nm == "1" else ("%s*%s" % (a, nm)))
    return " + ".join(parts).replace("+ -", "- ") if parts else "0"
