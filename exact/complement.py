"""complement -- the complement reading of an exact pattern, with no rounding.

CC0. Co-authored. Stdlib only. No floats.

Patterns and samples are Cyc5 points. The field at a sample is the SQUARED
distance to the nearest pattern point, a ZPhi; the smallest is chosen by
exact comparison. A symmetry of the field is exact equality at every sample.
"""
from fractions import Fraction

from .cyc5 import Cyc5, ZETA, norm2
from .zphi import PHI, ZPhi

Z10 = -(ZETA * ZETA * ZETA)        # a turn of 36 degrees


def samples():
    """A sample set closed under turns of 36 degrees and under the mirror in
    the real axis (so it favours no pattern that has fewer symmetries than
    the decagon). Orbits of four seed points."""
    seeds = [Cyc5(Fraction(1, 2), 0),
             Cyc5(Fraction(1, 2), 0) + ZETA * Fraction(1, 3),
             ZETA * (PHI - 1) + Cyc5(Fraction(1, 5), 0),
             Cyc5(Fraction(2, 3), 0) + Z10 * Fraction(1, 4)]
    out = set()
    for z in seeds:
        for _ in range(10):
            out.add(z)
            out.add(z.conjugate())
            z = z * Z10
    return sorted(out, key=lambda z: (z.u.a, z.u.b, z.v.a, z.v.b))


def field2(sample, pattern):
    """Squared distance from a sample to the nearest pattern point (a ZPhi)."""
    best = None
    for p in pattern:
        d = norm2(sample - p)
        if best is None or d < best:
            best = d
    return best


def field_invariant(pattern, pts, move):
    """True exactly when the field takes the same value at every sample and at
    its image under `move` (a function on Cyc5 points)."""
    return all(field2(s, pattern) == field2(move(s), pattern) for s in pts)


def read(pattern):
    pts = samples()
    turns = {}
    z = ZETA * ZETA.conjugate()          # 1, exactly
    for name, step in (("36", Z10), ("72", ZETA), ("144", ZETA * ZETA)):
        turns[name] = field_invariant(pattern, pts, lambda s, m=step: s * m)
    mirror = field_invariant(pattern, pts, lambda s: s.conjugate())
    values = {field2(s, pattern) for s in pts}
    return {"n_pattern_points": len(pattern), "n_samples": len(pts),
            "field_unchanged_by_turn_deg": turns,
            "field_unchanged_by_mirror_in_real_axis": mirror,
            "n_distinct_field_values": len(values),
            "arithmetic": "exact a + b*phi forms; equality is ==, no tolerance",
            "check_one_is_one": z == Cyc5(1, 0)}
