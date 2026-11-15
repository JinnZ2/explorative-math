"""complement -- read a shape through the space around it.

CC0. Co-authored. Stdlib only.

A shape (a finite point set) is collapsed to one point, the collapse point
(its centroid). What is left is the complement. The complement is sampled,
and the shape's pattern is read from the samples only: how far each sample
is from the pattern (the distance field), what turns and mirrors leave that
field unchanged, where the field has local maxima (voids), and the angular
gaps between pattern points seen from the collapse point.

  plane   a polar grid of samples around the collapse point
  sphere  for a shape in R^n: directions on the unit sphere; the field is
          the angle to the nearest vertex direction, and voids are found by
          a seeded hill climb

Every reading is bounded by its sampling; the bound is returned with it.
"""
import math
import random

from . import primitives as P


# ---- plane ---------------------------------------------------------------------

def centred_unit(points):
    """The pattern about its centroid, scaled so its farthest point is at 1."""
    c = P.centroid(points)
    q = [P.sub(p, c) for p in points]
    s = max((math.hypot(x, y) for x, y in q), default=0.0) or 1.0
    return [[x / s, y / s] for x, y in q]


def turned_to_zero(pattern):
    """The pattern turned about the collapse point so that its first farthest
    point lies at angle 0. Mirror axes of the field are tested on a grid of
    angles; without this a mirror would fall between grid angles."""
    far = max(pattern, key=lambda p: (round(math.hypot(*p), 9), -pattern.index(p)))
    a = -math.atan2(far[1], far[0])
    c, s_ = math.cos(a), math.sin(a)
    return [[c * x - s_ * y, s_ * x + c * y] for x, y in pattern]


def polar_grid(n_rings, n_theta, r_max=1.0):
    """Ring radii r_max * (i + 1/2) / n_rings, angles 360 * j / n_theta."""
    return [r_max * (i + 0.5) / n_rings for i in range(n_rings)], n_theta


def distance_field(pattern, radii, n_theta):
    """field[i][j] = distance from the sample on ring i at angle j to the
    nearest pattern point."""
    out = []
    for r in radii:
        row = []
        for j in range(n_theta):
            a = 2 * math.pi * j / n_theta
            x, y = r * math.cos(a), r * math.sin(a)
            row.append(math.sqrt(min((x - u) ** 2 + (y - v) ** 2 for u, v in pattern)))
        out.append(row)
    return out


def field_rotational_orders(field, n_theta, k_max=48, tol=1e-9):
    """Orders k (k divides n_theta, 2 <= k <= k_max) such that turning the
    field by 360/k leaves every sample unchanged within tol."""
    out = []
    for k in range(2, k_max + 1):
        if n_theta % k:
            continue
        s = n_theta // k
        if all(abs(row[j] - row[(j + s) % n_theta]) <= tol
               for row in field for j in range(n_theta)):
            out.append(k)
    return out


def field_mirror_axes(field, n_theta, tol=1e-9):
    """Mirror axes of the field through the collapse point, in degrees in
    [0, 180). Axes are tested at multiples of 180 / n_theta."""
    out = []
    for m in range(n_theta):
        ok = True
        for row in field:
            for j in range(n_theta):
                if abs(row[j] - row[(m - j) % n_theta]) > tol:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            out.append(180.0 * m / n_theta)
    return out


def angular_gaps(pattern, tol=1e-9):
    """Gaps (degrees) between distinct directions of pattern points seen from
    the collapse point, in order round the circle. Points at the centre have
    no direction and are left out."""
    th = sorted(math.degrees(math.atan2(y, x)) % 360.0 for x, y in pattern
                if math.hypot(x, y) > 1e-9)
    dirs = []
    for t in th:
        if not dirs or t - dirs[-1] > tol:
            dirs.append(t)
    if len(dirs) > 1 and dirs[0] + 360.0 - dirs[-1] <= tol:
        dirs.pop()
    if not dirs:
        return []
    return [(dirs[(i + 1) % len(dirs)] - dirs[i]) % 360.0 or 360.0 for i in range(len(dirs))]


def distinct_values(values, tol=1e-6):
    out = []
    for v in sorted(values):
        if not out or v - out[-1][0] > tol:
            out.append([v, 1])
        else:
            out[-1][1] += 1
    return [{"value": v, "count": c} for v, c in out]


def voids(field, radii, n_theta):
    """Grid points where the field is not below any of its four grid
    neighbours (the ring inside, the ring outside, the two angle steps), with
    strictly more than at least one. Returns them deepest first. The inner
    and outer rings have one radial neighbour fewer."""
    out = []
    n_r = len(radii)
    for i in range(n_r):
        for j in range(n_theta):
            v = field[i][j]
            nb = [field[i][(j - 1) % n_theta], field[i][(j + 1) % n_theta]]
            if i > 0:
                nb.append(field[i - 1][j])
            if i < n_r - 1:
                nb.append(field[i + 1][j])
            if all(v >= x for x in nb) and any(v > x for x in nb):
                out.append({"radius": radii[i], "angle_deg": 360.0 * j / n_theta, "depth": v,
                            "on_outer_ring": i == n_r - 1, "on_inner_ring": i == 0})
    out.sort(key=lambda r: -r["depth"])
    return out


def read_plane(points, n_rings=6, n_theta=1680, k_max=48, tol=1e-9):
    """The complement reading of a plane pattern."""
    pat = turned_to_zero(centred_unit(points))
    radii, n = polar_grid(n_rings, n_theta)
    f = distance_field(pat, radii, n)
    vs = voids(f, radii, n)
    inner = [v for v in vs if not v["on_outer_ring"] and not v["on_inner_ring"]]
    centre = math.sqrt(min(u * u + v * v for u, v in pat))
    gaps = angular_gaps(pat)
    flat = [x for row in f for x in row]
    return {"n_points": len(pat), "pattern": pat,
            "bound": {"n_rings": n_rings, "n_theta": n_theta, "angle_step_deg": 360.0 / n_theta,
                      "ring_radii": radii, "disc_radius": 1.0,
                      "orders_tested": [k for k in range(2, k_max + 1) if n_theta % k == 0],
                      "mirror_axis_step_deg": 180.0 / n_theta, "tol": tol,
                      "orientation": "pattern turned so its first farthest point is at 0",
                      "voids": "grid local maxima; those on the inner and outer rings are "
                               "marked and left out of the inside count, since each lacks "
                               "a radial neighbour"},
            "centre_depth": centre,
            "centre_is_a_void": all(centre >= x for x in f[0]),
            "field": f,
            "field_rotational_orders": field_rotational_orders(f, n, k_max, tol),
            "field_mirror_axes_deg": field_mirror_axes(f, n, tol),
            "gaps_deg": gaps, "distinct_gaps": distinct_values(gaps),
            "n_voids": len(vs), "n_voids_inside": len(inner),
            "voids": vs[:64],
            "void_depths_inside": distinct_values([v["depth"] for v in inner]),
            "deepest_void_inside": inner[0] if inner else None,
            "field_mean": sum(flat) / len(flat), "field_max": max(flat),
            "radial_profile": [sum(row) / len(row) for row in f]}


def field_difference(fa, fb, n_theta, theta_stride=8, shift_step=4):
    """Root-mean-square difference between two fields on the same grid, at the
    turn of the second that makes it smallest. Angles are subsampled by
    theta_stride and turns are tried every shift_step grid steps."""
    idx = range(0, n_theta, theta_stride)
    best, best_shift, zero = None, None, None
    for s in range(0, n_theta, shift_step):
        tot = cnt = 0
        for ra, rb in zip(fa, fb):
            for j in idx:
                d = ra[j] - rb[(j + s) % n_theta]
                tot += d * d
                cnt += 1
        rms = math.sqrt(tot / cnt)
        if s == 0:
            zero = rms
        if best is None or rms < best:
            best, best_shift = rms, s
    return {"rms_min_over_turns": best, "turn_deg": 360.0 * best_shift / n_theta,
            "rms_at_zero_turn": zero,
            "bound": {"theta_stride": theta_stride, "turn_step_deg": 360.0 * shift_step / n_theta}}


# ---- sphere (shape in R^n) ---------------------------------------------------------

def unit_dirs(vertices):
    return [P.unit(v) for v in vertices]


def angle_to_nearest(x, dirs):
    return math.acos(max(-1.0, min(1.0, max(P.dot(x, d) for d in dirs))))


def sphere_voids(vertices, seed=0, starts=60, steps=400, merge_deg=2.0, keep_deg=1.0):
    """Directions far from every vertex direction: seeded hill climb on the
    unit sphere (random moves that shrink, kept when the angle to the nearest
    vertex grows). Returns the distinct end points whose depth is within
    keep_deg of the deepest, and the bound."""
    dirs = unit_dirs(vertices)
    dim = len(dirs[0])
    rng = random.Random(seed)
    ends = []
    for _ in range(starts):
        x = P.unit([rng.gauss(0.0, 1.0) for _ in range(dim)])
        best = angle_to_nearest(x, dirs)
        step = 0.5
        for _ in range(steps):
            y = P.unit([a + step * rng.gauss(0.0, 1.0) for a in x])
            d = angle_to_nearest(y, dirs)
            if d > best:
                x, best = y, d
            else:
                step = max(step * 0.985, 1e-6)
        ends.append((best, x))
    ends.sort(key=lambda e: -e[0])
    top = ends[0][0]
    kept = []
    for depth, x in ends:
        if math.degrees(top - depth) > keep_deg:
            continue
        if all(math.degrees(math.acos(max(-1.0, min(1.0, P.dot(x, k["direction"])))))
               > merge_deg for k in kept):
            kept.append({"direction": x, "depth_deg": math.degrees(depth)})
    return {"deepest_depth_deg": math.degrees(top), "voids": kept, "n_voids_found": len(kept),
            "all_end_depths_deg": distinct_values([math.degrees(e[0]) for e in ends], tol=0.5),
            "bound": {"seed": seed, "starts": starts, "steps": steps, "merge_deg": merge_deg,
                      "keep_within_deg_of_deepest": keep_deg,
                      "found_is_a_lower_bound": "a void no start climbed to is not listed"}}


def match_directions(found, candidates):
    """How a set of found directions sits on a set of candidate directions:
    the angle from each found direction to its nearest candidate, and how
    many candidates have a found direction within 2 degrees."""
    cand = unit_dirs(candidates)
    if not found:
        return {"n_found": 0, "max_angle_deg": None, "mean_angle_deg": None,
                "candidates_hit": 0, "n_candidates": len(cand)}
    ang = [math.degrees(angle_to_nearest(f, cand)) for f in found]
    hit = sum(1 for c in cand
              if any(math.degrees(math.acos(max(-1.0, min(1.0, P.dot(c, f))))) <= 2.0
                     for f in found))
    return {"n_found": len(found), "max_angle_deg": max(ang),
            "mean_angle_deg": sum(ang) / len(ang),
            "candidates_hit": hit, "n_candidates": len(cand)}
