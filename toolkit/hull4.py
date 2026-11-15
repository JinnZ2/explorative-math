"""hull4 -- facets of the convex hull of points on a sphere in R^4, by brute force.

CC0. Co-authored. Stdlib only.

For points that all lie on one sphere about the origin (every one a vertex
of their hull). Every 4 of them span a hyperplane n . x = 1; it is a facet
plane when no point lies beyond it. The facet is all points on it, so a
facet with more than four vertices (a cube, an octahedron) is found whole.
"""
import itertools
import math

from . import primitives as P


def _solve(mat, rhs):
    n = len(rhs)
    m = [row[:] + [b] for row, b in zip(mat, rhs)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(m[r][i]))
        if abs(m[piv][i]) < 1e-12:
            return None
        m[i], m[piv] = m[piv], m[i]
        for r in range(i + 1, n):
            f = m[r][i] / m[i][i]
            for k in range(i, n + 1):
                m[r][k] -= f * m[i][k]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (m[i][n] - sum(m[i][k] * x[k] for k in range(i + 1, n))) / m[i][i]
    return x


def facets(points, max_edge_deg=None, tol=1e-7):
    """Facets of the hull, each a sorted tuple of point indices.

    max_edge_deg: only 4-subsets whose points are pairwise within this angle
    are tried (None tries every 4-subset). A facet with a longer edge would
    then be missed; closed() says whether the facets found close up."""
    n = len(points)
    near = None
    if max_edge_deg is not None:
        c = math.cos(math.radians(max_edge_deg))
        norms = [P.norm(p) for p in points]
        near = [[P.dot(points[i], points[j]) / (norms[i] * norms[j]) >= c - 1e-12
                 for j in range(n)] for i in range(n)]
    found = set()
    for quad in itertools.combinations(range(n), 4):
        if near is not None and not all(near[a][b] for a, b in itertools.combinations(quad, 2)):
            continue
        nvec = _solve([list(points[i]) for i in quad], [1.0] * 4)
        if nvec is None:
            continue
        vals = [P.dot(nvec, p) for p in points]
        if max(vals) > 1.0 + tol:
            continue
        found.add(tuple(i for i, v in enumerate(vals) if v > 1.0 - tol))
    return sorted(found)


def closed(fs):
    """For facets that are all tetrahedra: True when every triangle belongs to
    exactly two of them (the facets close up with no hole). None otherwise."""
    if not fs or any(len(f) != 4 for f in fs):
        return None
    count = {}
    for f in fs:
        for tri in itertools.combinations(f, 3):
            count[tri] = count.get(tri, 0) + 1
    return all(v == 2 for v in count.values())


def edge_pattern(points, facet):
    """Sorted pairwise angles (degrees, 2 decimals) between a facet's vertices."""
    out = []
    for a, b in itertools.combinations(facet, 2):
        c = P.dot(points[a], points[b]) / (P.norm(points[a]) * P.norm(points[b]))
        out.append(round(math.degrees(math.acos(max(-1.0, min(1.0, c)))), 2))
    return sorted(out)


def summary(points, max_edge_deg=None):
    fs = facets(points, max_edge_deg)
    sizes = {}
    for f in fs:
        sizes[len(f)] = sizes.get(len(f), 0) + 1
    patterns = {}
    for f in fs:
        key = tuple(edge_pattern(points, f))
        patterns[key] = patterns.get(key, 0) + 1
    edges = set()
    for f in fs:
        if len(f) == 4:
            edges.update(itertools.combinations(f, 2))
    return {"n_points": len(points), "n_facets": len(fs),
            "facets_by_vertex_count": {str(k): v for k, v in sorted(sizes.items())},
            "tetrahedral_facets_close_up": closed(fs),
            "n_edges_if_all_tetrahedra": len(edges) if all(len(f) == 4 for f in fs) else None,
            "edge_angle_patterns": [{"angles_deg": list(k), "facets": v}
                                    for k, v in sorted(patterns.items())],
            "bound": {"max_edge_deg": max_edge_deg, "tol": 1e-7,
                      "method": "every 4-subset (within max_edge_deg when given); a plane is a "
                                "facet plane when no point lies beyond it"}}
