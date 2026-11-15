"""probes -- measurements taken on a shape or a sequence.

CC0. Co-authored. Stdlib only.

  rotational_symmetry  k-fold symmetries of a 2D point set about its centroid
  mirror_symmetry      whether the 2D set has a mirror line through the centroid
                       (axes on a fixed angle grid)
  mirror_symmetry_from_pairs  the same with no grid: candidate axes from
                       every pair of points
  smallest_period      smallest p with s[i] == s[i+p] throughout
  period_profile       smallest period of growing prefixes (aperiodicity shows
                       as growth; a finite sequence always has some period)
  zlib_ratio           compressed length / raw length of the byte encoding
  angle_histogram      counts of angles in fixed-width bins
  three_gap_sum        for three gap lengths: is the longest the sum of the
                       other two (three-distance theorem)? pass / fail, or
                       not_applicable for any other count
  min_stretch_tag      is a gap reading taken on fewer points than about
                       3 / (frequency of the third gap)? "under-resolved"
                       or no tag; frequency unknown -> not assessed
  rotational_symmetry_bound, mirror_symmetry_bound, period_bound
                       the instrument bound a reading was taken at (points,
                       precision, probe step), logged beside the reading
                       (docs/resolution_leak.md)
"""
import math
import zlib

from . import primitives as P


def _matches(a, b, tol):
    used = [False] * len(b)
    for p in a:
        hit = None
        for j, q in enumerate(b):
            if not used[j] and math.dist(p, q) <= tol:
                hit = j
                break
        if hit is None:
            return False
        used[hit] = True
    return True


def _scale_of(points):
    c = P.centroid(points)
    return max((P.norm(P.sub(p, c)) for p in points), default=0.0) or 1.0


def rotational_symmetry(points, k_max=12, rel_tol=1e-6):
    """List of k in 2..k_max such that rotating by 360/k maps the set to itself."""
    if len(points) < 2:
        return []
    c = P.centroid(points)
    pts = [P.sub(p, c) for p in points]
    tol = rel_tol * _scale_of(points)
    return [k for k in range(2, k_max + 1)
            if _matches(P.rotate(pts, 0, 1, 2 * math.pi / k), pts, tol)]


def mirror_symmetry(points, n_axes=360, rel_tol=1e-6):
    """Axis angles (degrees, in [0, 180)) of mirror lines through the centroid."""
    c = P.centroid(points)
    pts = [P.sub(p, c) for p in points]
    tol = rel_tol * _scale_of(points)
    out = []
    for m in range(n_axes):
        a = math.pi * m / n_axes
        c2, s2 = math.cos(2 * a), math.sin(2 * a)
        refl = [[c2 * x + s2 * y, s2 * x - c2 * y] for x, y in pts]
        if _matches(refl, pts, tol):
            out.append(180.0 * m / n_axes)
    return out


def mirror_symmetry_from_pairs(points, rel_tol=1e-6):
    """Mirror axes through the centroid with no angle grid.

    A mirror line through the centroid maps each point to a point of the set,
    so its angle is the mean polar angle of some pair (a point with itself
    included). Every such candidate is tested. Returns axis angles in degrees
    in [0, 180), merged within 1e-9 degrees. Bounded by rel_tol and by the
    centroid as the centre, not by an axis step."""
    c = P.centroid(points)
    pts = [P.sub(p, c) for p in points]
    tol = rel_tol * _scale_of(points)
    th = [math.atan2(y, x) for x, y in pts if math.hypot(x, y) > tol]
    cands = sorted({round(math.degrees((a + b) / 2.0) % 180.0, 9)
                    for i, a in enumerate(th) for b in th[i:]})
    out = []
    for deg in cands:
        a = math.radians(deg)
        c2, s2 = math.cos(2 * a), math.sin(2 * a)
        refl = [[c2 * x + s2 * y, s2 * x - c2 * y] for x, y in pts]
        if _matches(refl, pts, tol) and not (out and deg - out[-1] <= 1e-9):
            out.append(deg)
    if len(out) > 1 and out[-1] - out[0] >= 180.0 - 1e-9:
        out.pop()
    return out


def smallest_period(s):
    n = len(s)
    for p in range(1, n):
        if all(s[i] == s[i + p] for i in range(n - p)):
            return p
    return n


def period_profile(s, lengths=None):
    lengths = lengths or [k for k in (8, 16, 32, 64, 128, 256, 512) if k <= len(s)]
    return {str(k): smallest_period(s[:k]) for k in lengths}


def zlib_ratio(data):
    """Compressed size / raw size at level 9. Bytes, or a str (UTF-8 encoded)."""
    if isinstance(data, str):
        data = data.encode("utf-8")
    if not data:
        return None
    return len(zlib.compress(data, 9)) / len(data)


def angle_histogram(angles_deg, width=10.0, lo=0.0, hi=180.0):
    nb = int(round((hi - lo) / width))
    counts = [0] * nb
    for a in angles_deg:
        b = int((a - lo) // width)
        if 0 <= b < nb:
            counts[b] += 1
        elif a == hi:
            counts[-1] += 1
    return {"lo": lo, "hi": hi, "width": width, "counts": counts}


def polar_angles(points2d):
    c = P.centroid(points2d)
    return sorted(math.degrees(math.atan2(p[1] - c[1], p[0] - c[0])) % 360.0 for p in points2d)



def three_gap_sum(lengths, tol=1e-6):
    """Consistency check on a set of distinct gap lengths.

    With exactly three lengths L1 <= L2 <= L3, the three-distance theorem
    gives L3 = L1 + L2. Returns the verdict "pass" or "fail" with the
    residual L3 - (L1 + L2). Any other count is "not_applicable" (residual
    None): it is recorded, not treated as a pass or a fail."""
    ls = sorted(float(x) for x in lengths)
    if len(ls) != 3:
        return {"n_lengths": len(ls), "verdict": "not_applicable", "residual": None,
                "tol": tol}
    r = ls[2] - (ls[0] + ls[1])
    return {"n_lengths": 3, "verdict": "pass" if abs(r) <= tol else "fail",
            "residual": r, "tol": tol}


def min_stretch_tag(n_points, third_gap_freq, factor=3.0):
    """Bound check on a gap-length reading.

    third_gap_freq: share of gaps that have the rarest of three lengths, from
    whatever reading is the reference for this window (declare which). A
    stretch shorter than factor / third_gap_freq points is tagged
    "under-resolved". None or 0 for the frequency means it is not known (or
    no third gap has been seen): the reading is not assessed and gets no tag."""
    if not third_gap_freq:
        return {"n_points": n_points, "third_gap_freq": third_gap_freq or None,
                "min_points": None, "factor": factor, "assessed": False, "tag": None}
    need = factor / third_gap_freq
    return {"n_points": n_points, "third_gap_freq": third_gap_freq, "min_points": need,
            "factor": factor, "assessed": True,
            "tag": "under-resolved" if n_points < need else None}


# ---- instrument bounds: logged beside a reading, never instead of it --------

def rotational_symmetry_bound(points, k_max=12, rel_tol=1e-6):
    """The bound rotational_symmetry reads at, for the same arguments."""
    return {"n_points": len(points), "k_tested": [2, k_max], "rel_tol": rel_tol,
            "abs_tol": rel_tol * _scale_of(points) if points else None,
            "about": "centroid"}


def mirror_symmetry_bound(points, n_axes=360, rel_tol=1e-6):
    """The bound mirror_symmetry reads at, for the same arguments."""
    return {"n_points": len(points), "n_axes": n_axes, "axis_step_deg": 180.0 / n_axes,
            "rel_tol": rel_tol,
            "abs_tol": rel_tol * _scale_of(points) if points else None,
            "about": "centroid"}


def period_bound(s, lengths=None):
    """The bound smallest_period / period_profile read at: the length of the
    word, and the prefix lengths the profile looks at."""
    lengths = lengths or [k for k in (8, 16, 32, 64, 128, 256, 512) if k <= len(s)]
    return {"n": len(s), "prefix_lengths": list(lengths)}
