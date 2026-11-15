"""fold -- reflect a point set into the fundamental wedge of its mirrors, exactly.

CC0. Co-authored. Stdlib only. No floats, no tolerance: two folded points are
the same point only if they are equal.

Points are Cyc5 (see cyc5.py). A mirror is a line through the origin, given
by a direction u (a Cyc5). The wedge of a set of mirrors is bounded by two
of them: directions lo and hi, with the wedge on the anticlockwise side of
lo and the clockwise side of hi.
"""
from .cyc5 import Cyc5, cross, norm2
from .zphi import ZERO


def reflect(z, u):
    """z reflected in the line through 0 with direction u: u^2 * conj(z) / |u|^2."""
    return (u * u * z.conjugate()) * norm2(u).inverse()


def side(z, u):
    """1 if z is anticlockwise of the direction u, -1 if clockwise, 0 if on the line."""
    return cross(u, z).sign()


def fold_one(z, lo, hi, max_steps=64):
    """Reflect z across lo and hi in turn until it lies in the wedge between
    them (anticlockwise of lo, clockwise of hi, edges included). Returns the
    folded point and how many reflections it took.

    Not checked here: that lo and hi are mirrors of anything. A point can be
    reflected into the space between any two lines. Test the mirrors with
    is_mirror before folding a set by them."""
    for k in range(max_steps):
        if side(z, lo) < 0:
            z = reflect(z, lo)
        elif side(z, hi) > 0:
            z = reflect(z, hi)
        else:
            return z, k
    raise ValueError("did not reach the wedge within max_steps")


def fold(points, lo, hi):
    """Fold every point into the wedge. Returns [(folded point, multiplicity)],
    in the order the folded points first appear; multiplicity is how many of
    the given points land there. Exact equality, no tolerance."""
    out, index = [], {}
    for z in points:
        f, _ = fold_one(z, lo, hi)
        if f in index:
            out[index[f]][1] += 1
        else:
            index[f] = len(out)
            out.append([f, 1])
    return [(f, m) for f, m in out]


def half_plane_fold(points, u):
    """Fold by ONE mirror: every point clockwise of u is reflected across it."""
    out, index = [], {}
    for z in points:
        f = reflect(z, u) if side(z, u) < 0 else z
        if f in index:
            out[index[f]][1] += 1
        else:
            index[f] = len(out)
            out.append([f, 1])
    return [(f, m) for f, m in out]


def is_mirror(points, u):
    """True exactly when reflecting in the line with direction u maps the set
    of points onto itself (as a set)."""
    s = set(points)
    return all(reflect(z, u) in s for z in s)


def on_ray(z, u):
    """True when z is a non-negative multiple of u (or the origin)."""
    if z == Cyc5(0, 0):
        return True
    from .cyc5 import dot
    return cross(u, z) == ZERO and dot(u, z).sign() > 0
