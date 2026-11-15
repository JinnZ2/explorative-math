"""primitives -- vectors, plane rotations, projection, angles, polytopes, lattices.

CC0. Co-authored. Stdlib only. Points are lists of floats; a shape is a list
of points.
"""
import itertools
import math

PHI = (1 + 5 ** 0.5) / 2


# ---- vectors ---------------------------------------------------------------

def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    n = norm(a)
    if n == 0:
        raise ValueError("zero vector has no direction")
    return [x / n for x in a]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def scale(a, k):
    return [k * x for x in a]


def centroid(points):
    n = len(points)
    return [sum(p[i] for p in points) / n for i in range(len(points[0]))]


def random_unit(rng, d):
    while True:
        v = [rng.gauss(0.0, 1.0) for _ in range(d)]
        if norm(v) > 0:
            return unit(v)


def angle_deg(a, b):
    c = dot(a, b) / (norm(a) * norm(b))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def gram_schmidt(vs, tol=1e-12):
    basis = []
    for v in vs:
        w = list(v)
        for b in basis:
            w = sub(w, scale(b, dot(w, b)))
        if norm(w) > tol:
            basis.append(unit(w))
    return basis


# ---- rotations and projection ------------------------------------------------

def rotate(points, i, j, theta):
    """Rotate every point by theta in the (i, j) coordinate plane. Any dimension."""
    c, s = math.cos(theta), math.sin(theta)
    out = []
    for p in points:
        q = list(p)
        q[i], q[j] = c * p[i] - s * p[j], s * p[i] + c * p[j]
        out.append(q)
    return out


def rotate_many(points, plane_angles):
    """Apply rotations in order; plane_angles is [((i, j), theta), ...]."""
    for (i, j), th in plane_angles:
        points = rotate(points, i, j, th)
    return points


def all_planes(d):
    return list(itertools.combinations(range(d), 2))


def project(points, basis=None, coords=(0, 1)):
    """Project onto an orthonormal basis (list of vectors), or else keep `coords`."""
    if basis is not None:
        return [[dot(p, b) for b in basis] for p in points]
    return [[p[c] for c in coords] for p in points]


def reflect(v, root):
    """Reflect v in the hyperplane orthogonal to root."""
    return sub(v, scale(root, 2.0 * dot(v, root) / dot(root, root)))


def coxeter_plane(simple_roots, h):
    """Orthonormal basis (two vectors) of the Coxeter plane of a reflection group.

    simple_roots: the simple roots. h: the Coxeter number (the order of the
    Coxeter element c = product of the simple reflections). The plane is the
    one c rotates by 360/h degrees. It is taken as the image of
    (2/h) * sum_k cos(2 pi k / h) c^k, which projects onto that plane when the
    eigenvalues exp(+-2 pi i / h) of c are simple. Raises ValueError if c^h
    is not the identity or the image is not two-dimensional (both at 1e-9)."""
    d = len(simple_roots[0])

    def c(v):
        for r in simple_roots:
            v = reflect(v, r)
        return v
    cols = []
    for i in range(d):
        e = [1.0 if j == i else 0.0 for j in range(d)]
        acc, v = [0.0] * d, e
        for k in range(h):
            acc = [a + b * math.cos(2 * math.pi * k / h) for a, b in zip(acc, v)]
            v = c(v)
        if max(abs(a - b) for a, b in zip(v, e)) > 1e-9:
            raise ValueError("Coxeter element does not have order h")
        cols.append(scale(acc, 2.0 / h))
    basis = gram_schmidt(cols, tol=1e-9)
    if len(basis) != 2:
        raise ValueError(f"expected a plane, got dimension {len(basis)}")
    return basis


# ---- polytope vertex generators ---------------------------------------------

def simplex(n):
    """Regular n-simplex in R^n: n+1 vertices, centred, edge sqrt(2)."""
    e = [[1.0 if i == k else 0.0 for i in range(n + 1)] for k in range(n + 1)]
    c = 1.0 / (n + 1)
    centred = [[x - c for x in v] for v in e]
    basis = gram_schmidt(centred[:n])
    return [[dot(v, b) for b in basis] for v in centred]


def hypercube(n):
    """2^n vertices (+-1, ..., +-1)."""
    return [list(v) for v in itertools.product((-1.0, 1.0), repeat=n)]


def cross_polytope(n):
    """2n vertices +-e_i."""
    out = []
    for i in range(n):
        for s in (1.0, -1.0):
            v = [0.0] * n
            v[i] = s
            out.append(v)
    return out


def cell24():
    """24-cell: all permutations of (+-1, +-1, 0, 0)."""
    out = set()
    for i, j in itertools.combinations(range(4), 2):
        for si, sj in itertools.product((-1.0, 1.0), repeat=2):
            v = [0.0] * 4
            v[i], v[j] = si, sj
            out.add(tuple(v))
    return [list(v) for v in sorted(out)]


def icosahedron():
    """12 vertices: cyclic permutations of (0, +-1, +-phi)."""
    out = []
    for a, b in itertools.product((-1.0, 1.0), repeat=2):
        base = (0.0, a, b * PHI)
        for k in range(3):
            out.append([base[(0 - k) % 3], base[(1 - k) % 3], base[(2 - k) % 3]])
    return out


def regular_polygon(k, r=1.0, phase=0.0):
    return [[r * math.cos(phase + 2 * math.pi * m / k), r * math.sin(phase + 2 * math.pi * m / k)]
            for m in range(k)]


POLYTOPES = {"simplex": simplex, "hypercube": hypercube, "cross_polytope": cross_polytope}
FIXED_SHAPES = {"cell24": cell24, "icosahedron": icosahedron}


# ---- lattices ------------------------------------------------------------------

def lattice_zn(n, extent):
    """Integer points of Z^n with every |coordinate| <= extent."""
    return [list(map(float, v)) for v in itertools.product(range(-extent, extent + 1), repeat=n)]


def lattice_hex(extent):
    """Triangular (hexagonal) lattice in the plane, basis (1, 0), (1/2, sqrt3/2)."""
    h = math.sqrt(3) / 2
    return [[m + 0.5 * k, h * k] for m in range(-extent, extent + 1)
            for k in range(-extent, extent + 1)]


def edge_lengths(points):
    return [math.dist(a, b) for a, b in itertools.combinations(points, 2)]


def min_edge(points):
    return min(edge_lengths(points))
