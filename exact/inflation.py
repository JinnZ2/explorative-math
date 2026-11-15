"""inflation -- substitution tilings carried exactly. CC0. Co-authored. Stdlib only.

No floats. Lengths are ZPhi, points are Cyc5.

  Fibonacci chain   tiles L (length phi) and S (length 1); L -> LS, S -> L.
                    Inflation by phi: an L of length phi becomes length
                    phi^2 = phi + 1 = L + S; an S of length 1 becomes
                    length phi = L.
  pentagon          vertices zeta^k; the pentagram's inner pentagon is the
                    outer one times -1/phi^2.
  Penrose           Robinson triangles (colour 0 acute, colour 1 obtuse),
                    the standard subdivision with ratio 1/phi = phi - 1.
"""
from fractions import Fraction

from .cyc5 import Cyc5, ZETA, cross, norm2
from .zphi import ZPhi, PHI, ONE, ZERO

INV_PHI = PHI - 1                 # 1/phi, exactly
TILE = {"L": PHI, "S": ONE}
RULE = {"L": "LS", "S": "L"}


# ---- Fibonacci chain ---------------------------------------------------------

def fib_substitute(word):
    return "".join(RULE[c] for c in word)


def fib_positions(word):
    """Left end of every tile, then the right end of the last: len(word) + 1 ZPhi."""
    out, x = [ZERO], ZERO
    for c in word:
        x = x + TILE[c]
        out.append(x)
    return out


def fib_inflate(word, positions):
    """Geometric inflation: scale every position by phi, then cut each scaled L
    into L + S. Returns (new word, new positions). The cut point of a scaled
    L starting at x is x + phi."""
    new_word, new_pos = [], []
    for c, x in zip(word, positions):
        x = x * PHI
        new_pos.append(x)
        if c == "L":
            new_word.append("LS")
            new_pos.append(x + PHI)
        else:
            new_word.append("L")
    new_pos.append(positions[-1] * PHI)
    return "".join(new_word), new_pos


# ---- pentagon ------------------------------------------------------------------

def pentagon(scale=ONE):
    """Vertices scale * zeta^k, k = 0..4 (scale a ZPhi, may be negative)."""
    out, z = [], Cyc5(1, 0)
    for _ in range(5):
        out.append(z * scale)
        z = z * ZETA
    return out


def on_line(q, a, b):
    """True exactly when q lies on the line through a and b."""
    return cross(q - a, b - a) == ZERO


# ---- Penrose (Robinson triangles) ----------------------------------------------

def penrose_seed():
    """One acute triangle: apex 0, legs 1, at 0 and 36 degrees."""
    z10 = -(ZETA * ZETA * ZETA)            # exp(i pi / 5)
    return [(0, Cyc5(0, 0), Cyc5(1, 0), z10)]


def penrose_subdivide(tris):
    out = []
    for colour, a, b, c in tris:
        if colour == 0:
            p = a + (b - a) * INV_PHI
            out.append((0, c, p, b))
            out.append((1, p, c, a))
        else:
            q = b + (a - b) * INV_PHI
            r = b + (c - b) * INV_PHI
            out.append((1, r, c, a))
            out.append((1, q, r, b))
            out.append((0, r, q, a))
    return out


def area_over_w(tri):
    """Signed area of a triangle divided by |w|, a ZPhi (see cyc5.cross)."""
    _, a, b, c = tri
    return cross(b - a, c - a) * Fraction(1, 2)


def edge_norms(tri):
    _, a, b, c = tri
    return [norm2(b - a), norm2(c - b), norm2(a - c)]


def phi_power(x, k_range=range(-80, 81)):
    """The integer k with x == phi^k, or None."""
    for k in k_range:
        if x == PHI ** k:
            return k
    return None
