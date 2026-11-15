"""hexsheet -- hexagonal cells relaxing under a field. CC0. Co-authored. Stdlib only.

CONSTRUCTED toy. Nothing here is a claim about a real sheet.

  ideal reference   cell centres on the triangular lattice, held exactly as
                    axial integers (q, r). The ideal position is
                    x = q + r/2, y = r * (sqrt(3)/2): x is an exact
                    Fraction, y is an exact Fraction times the fixed
                    irrational sqrt(3)/2. The ideal cell is the regular
                    hexagon whose six neighbour offsets have length 1.
  field             a scalar value at every cell, from its ideal position.
                    gradient: value = slope * (x - x_centre), mean zero
                    over a symmetric sheet. uniform: the same value
                    everywhere.
  local rule        each cell has a preferred spacing set by the field
                    where it sits: 1 + alpha * value. The rest length
                    between two neighbours is the mean of their two
                    preferred spacings. Every cell moves a fraction of the
                    net pull of its neighbours, all cells together, a fixed
                    number of sweeps: it relaxes toward balance with its
                    neighbours under its local field value.
  response          how a cell differs from the ideal after relaxing. It is
                    the cell's answer to the field where it sits, and is
                    named response everywhere, not error.

Positions and responses are binary floats; the ideal reference is exact.
"""
import math
from fractions import Fraction

ROOT3_2 = math.sqrt(3.0) / 2.0
NEIGHBOURS = ((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1))    # axial, anticlockwise


def cells(radius):
    """Axial coordinates of a hexagon-shaped sheet: every (q, r) within
    `radius` steps of the centre. 3*radius*(radius+1) + 1 cells."""
    out = []
    for q in range(-radius, radius + 1):
        for r in range(max(-radius, -q - radius), min(radius, -q + radius) + 1):
            out.append((q, r))
    return out


def ideal_exact(cell):
    """Exact ideal position: (x, y_units) with y = y_units * sqrt(3)/2."""
    q, r = cell
    return Fraction(2 * q + r, 2), Fraction(r)


def ideal_xy(cell):
    x, yu = ideal_exact(cell)
    return [float(x), float(yu) * ROOT3_2]


def ideal_offsets():
    """The six neighbour offsets of the ideal cell (length 1 each)."""
    return [ideal_xy(n) for n in NEIGHBOURS]


def gradient_field(sheet, slope):
    """value = slope * x, x the ideal position; the sheet is centred on x = 0,
    so the mean over a symmetric sheet is zero."""
    return {c: slope * ideal_xy(c)[0] for c in sheet}


def uniform_field(sheet, value=0.0):
    return {c: value for c in sheet}


def edges(sheet):
    s = set(sheet)
    out = []
    for q, r in sheet:
        for dq, dr in NEIGHBOURS[:3]:
            n = (q + dq, r + dr)
            if n in s:
                out.append(((q, r), n))
    return out


def relax(sheet, field, alpha, sweeps, rate=0.15):
    """The local rule, applied `sweeps` times to all cells together.

    Returns (positions, bound): the bound holds what the result was read at
    (sweeps, rate, the largest net pull left on any cell, the largest move
    in the last sweep)."""
    pos = {c: ideal_xy(c) for c in sheet}
    es = [(a, b, 1.0 + alpha * (field[a] + field[b]) / 2.0) for a, b in edges(sheet)]
    max_pull = last_move = 0.0
    for _ in range(sweeps):
        pull = {c: [0.0, 0.0] for c in sheet}
        for a, b, rest in es:
            dx, dy = pos[b][0] - pos[a][0], pos[b][1] - pos[a][1]
            d = math.hypot(dx, dy)
            f = (d - rest) / d
            pull[a][0] += f * dx
            pull[a][1] += f * dy
            pull[b][0] -= f * dx
            pull[b][1] -= f * dy
        max_pull = max(math.hypot(*p) for p in pull.values())
        last_move = rate * max_pull
        for c in sheet:
            pos[c][0] += rate * pull[c][0]
            pos[c][1] += rate * pull[c][1]
    # after the last move, the pull that is left
    pull = {c: [0.0, 0.0] for c in sheet}
    for a, b, rest in es:
        dx, dy = pos[b][0] - pos[a][0], pos[b][1] - pos[a][1]
        d = math.hypot(dx, dy)
        f = (d - rest) / d
        pull[a][0] += f * dx
        pull[a][1] += f * dy
        pull[b][0] -= f * dx
        pull[b][1] -= f * dy
    left = max(math.hypot(*p) for p in pull.values()) if sheet else 0.0
    return pos, {"sweeps": sweeps, "rate": rate, "n_cells": len(sheet), "n_edges": len(es),
                 "max_net_pull_left": left, "max_move_in_last_sweep": last_move,
                 "arithmetic": "binary float (64-bit)"}


def responses(sheet, field, pos):
    """One record per cell: its field value and its response.

    size_response    mean distance to its neighbours, minus 1 (the ideal).
    shape_response   spread (population standard deviation) of those
                     distances: 0 for a regular hexagon of any size.
    offset_response  each neighbour offset minus the ideal offset, in the
                     NEIGHBOURS order; null where the neighbour is off the
                     sheet.
    A cell with fewer than six neighbours is on the edge of the sheet; it is
    logged with what it has and marked interior: false."""
    s = set(sheet)
    ideal = ideal_offsets()
    out = []
    for c in sheet:
        offs, dists = [], []
        for (dq, dr), io in zip(NEIGHBOURS, ideal):
            n = (c[0] + dq, c[1] + dr)
            if n in s:
                v = [pos[n][0] - pos[c][0], pos[n][1] - pos[c][1]]
                offs.append([v[0] - io[0], v[1] - io[1]])
                dists.append(math.hypot(*v))
            else:
                offs.append(None)
        mean = sum(dists) / len(dists) if dists else None
        spread = (math.sqrt(sum((d - mean) ** 2 for d in dists) / len(dists))
                  if dists else None)
        x, yu = ideal_exact(c)
        out.append({"cell": list(c),
                    "ideal_exact": [[x.numerator, x.denominator], [yu.numerator, yu.denominator]],
                    "field_value": field[c],
                    "interior": len(dists) == 6, "n_neighbours": len(dists),
                    "position": list(pos[c]),
                    "size_response": (mean - 1.0) if mean is not None else None,
                    "shape_response": spread,
                    "offset_response": offs})
    return out


def pearson(xs, ys):
    """Pearson correlation; None when either side has no spread."""
    n = len(xs)
    if n < 2:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0 or syy == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(sxx * syy)


def slope(xs, ys):
    """Least-squares slope of y on x; None when x has no spread."""
    n = len(xs)
    if n < 2:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx


def mean_cell(records):
    """The mean hexagon over interior cells: the mean of each of the six
    neighbour offsets, set beside the ideal offsets."""
    inner = [r for r in records if r["interior"]]
    ideal = ideal_offsets()
    if not inner:
        return None
    mean_dev = [[sum(r["offset_response"][k][i] for r in inner) / len(inner) for i in (0, 1)]
                for k in range(6)]
    mean_off = [[io[0] + d[0], io[1] + d[1]] for io, d in zip(ideal, mean_dev)]
    return {"n_interior_cells": len(inner), "ideal_offsets": ideal,
            "mean_offsets": mean_off, "mean_offset_response": mean_dev,
            "max_abs_mean_offset_response": max(abs(v) for d in mean_dev for v in d),
            "mean_side_lengths": [math.hypot(*o) for o in mean_off],
            # an offset and its opposite sit on either side of the cell; their mean
            # length is the cell's own spacing along that direction
            "mean_opposite_pair_lengths": [
                (math.hypot(*mean_off[k]) + math.hypot(*mean_off[k + 3])) / 2.0
                for k in range(3)],
            "largest_single_cell_offset_response":
                max(abs(v) for r in inner for d in r["offset_response"] for v in d)}
