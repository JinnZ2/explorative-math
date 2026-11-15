"""s. exact_inflation -- substitution tilings in exact a + b*phi arithmetic.

CC0. Co-authored. Stdlib only.

Status of the claim: PROPOSED before the run -- carried in exact forms
(exact/: a + b*phi with rational a, b; phi^2 = phi + 1; no floats) an
inflation can be repeated without drift. Three checks, each over LEVELS
levels:

  fibonacci   the chain built two ways at every level: by substituting the
              word and summing exact tile lengths, and by scaling the
              previous level's positions by phi and cutting. The two lists
              of positions must be EQUAL, not close.
  pentagon    nested pentagrams: the inner pentagon is the outer one times
              -1/phi^2. Each inner vertex must lie exactly on two
              diagonals of the outer pentagon; diagonal^2 / side^2 must be
              exactly phi^2.
  penrose     Robinson-triangle subdivision from one triangle: total area
              conserved exactly; every squared edge length exactly a power
              of phi; triangle counts on the Fibonacci recurrence; vertices
              de-duplicated by exact equality, no tolerance.

Beside the Fibonacci check the same two constructions are run in binary
floats, and the largest difference between them is logged as a contrast.
The floats are in this script only; exact/ has none.

  python3 experiments/exact_inflation.py [seed]
  python3 experiments/exact_inflation.py --small [seed]   # phone bound: fewer levels,
                                                          # every output kept, own run file
"""
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from exact import PHI, ONE, ZERO, ZPhi, norm2  # noqa: E402
from exact import inflation as I  # noqa: E402

NAME = "exact_inflation"
LEVELS = {"fibonacci": 16, "pentagon": 12, "penrose": 10}
SMALL = {"fibonacci": 10, "pentagon": 10, "penrose": 6}
CLAIM = ("PROPOSED before the run: in exact a + b*phi forms an inflation repeats with zero "
         "drift; the two constructions of each level agree exactly, at every level.")


def fibonacci(levels):
    word = "L"
    pos = I.fib_positions(word)
    fword, fpos = "L", [0.0, 1.618033988749895]
    fphi = 1.618033988749895
    rows = []
    for n in range(1, levels + 1):
        gw, gp = I.fib_inflate(word, pos)             # geometric: scale and cut
        sw = I.fib_substitute(word)                   # combinatorial: substitute, then sum
        sp = I.fib_positions(sw)
        mism = sum(1 for a, b in zip(gp, sp) if a != b)
        # float contrast: the same two constructions in binary floats
        fg = []
        for c, x in zip(fword, fpos):
            x = x * fphi
            fg.append(x)
            if c == "L":
                fg.append(x + fphi)
        fg.append(fpos[-1] * fphi)
        fs, x = [0.0], 0.0
        for c in sw:
            x += fphi if c == "L" else 1.0
            fs.append(x)
        fdiff = max(abs(a - b) for a, b in zip(fg, fs))
        rows.append({"level": n, "n_tiles": len(sw), "count_L": sw.count("L"),
                     "count_S": sw.count("S"),
                     "words_equal": gw == sw,
                     "positions_equal": gp == sp and len(gp) == len(sp),
                     "n_position_mismatches": mism,
                     "total_length": sp[-1].pair(),
                     "total_length_equals_phi_power": sp[-1] == PHI ** (n + 1),
                     "float_contrast_max_abs_difference": fdiff,
                     "float_contrast_n_unequal": sum(1 for a, b in zip(fg, fs) if a != b)})
        word, pos = gw, gp
        fword, fpos = sw, fg
    return {"levels": levels, "tile_lengths": {"L": PHI.pair(), "S": ONE.pair()},
            "rows": rows,
            "zero_drift": all(r["positions_equal"] and r["words_equal"] for r in rows),
            "final_word_head": word[:64],
            "final_positions_head": [p.pair() for p in pos[:8]]}


def pentagon(levels):
    rows = []
    scale = ONE
    for n in range(levels + 1):
        outer = I.pentagon(scale)
        inner_scale = scale * (-(I.INV_PHI * I.INV_PHI))
        inner = I.pentagon(inner_scale)
        diagonals = [(outer[i], outer[(i + 2) % 5]) for i in range(5)]
        on = [sum(1 for a, b in diagonals if I.on_line(q, a, b)) for q in inner]
        side2 = norm2(outer[1] - outer[0])
        diag2 = norm2(outer[2] - outer[0])
        rows.append({"level": n, "scale": scale.pair(),
                     "diagonals_through_each_inner_vertex": on,
                     "every_inner_vertex_on_two_diagonals": all(k == 2 for k in on),
                     "diag2_over_side2_is_phi2": diag2 == side2 * PHI * PHI,
                     "side2": side2.pair(),
                     "side2_ratio_to_next_is_phi4":
                         side2 == norm2(inner[1] - inner[0]) * PHI ** 4,
                     "vertex_sum_is_zero": sum(outer[1:], outer[0]) == outer[0] - outer[0]})
        scale = inner_scale
    return {"levels": levels, "rows": rows,
            "zero_drift": all(r["every_inner_vertex_on_two_diagonals"]
                              and r["diag2_over_side2_is_phi2"]
                              and r["side2_ratio_to_next_is_phi4"]
                              and r["vertex_sum_is_zero"] for r in rows)}


def penrose(levels):
    tris = I.penrose_seed()
    area0 = I.area_over_w(tris[0])
    rows = []
    n0, n1 = 1, 0
    for n in range(levels + 1):
        if n:
            tris = I.penrose_subdivide(tris)
            n0, n1 = n0 + n1, n0 + 2 * n1
        c0 = sum(1 for t in tris if t[0] == 0)
        c1 = len(tris) - c0
        area = sum((abs(I.area_over_w(t)) for t in tris), ZERO)
        norms = set()
        for t in tris:
            norms.update(I.edge_norms(t))
        powers = sorted(I.phi_power(x) if I.phi_power(x) is not None else 10 ** 9
                        for x in norms)
        verts = set()
        for _, a, b, c in tris:
            verts.update((a, b, c))
        rows.append({"level": n, "n_triangles": len(tris), "n_acute": c0, "n_obtuse": c1,
                     "counts_follow_recurrence": (c0, c1) == (n0, n1),
                     "area_over_w": area.pair(),
                     "area_conserved_exactly": area == abs(area0),
                     "n_distinct_squared_edge_lengths": len(norms),
                     "squared_edge_lengths_as_phi_powers":
                         [p if p != 10 ** 9 else None for p in powers],
                     "every_squared_edge_is_a_phi_power": all(p != 10 ** 9 for p in powers),
                     "n_distinct_vertices_exact": len(verts)})
    return {"levels": levels, "seed_triangle": [p.pair() for p in tris[0][1:]] if levels == 0
            else [p.pair() for p in I.penrose_seed()[0][1:]],
            "rows": rows,
            "zero_drift": all(r["area_conserved_exactly"] and r["counts_follow_recurrence"]
                              and r["every_squared_edge_is_a_phi_power"] for r in rows)}


def run(seed=0, levels=None):
    lv = dict(LEVELS if levels is None else levels)
    return {"claim": CLAIM,
            "bound": {"levels": lv, "arithmetic": "fractions.Fraction pairs (a, b) for a + b*phi; "
                      "points u + v*w with w^2 = -(2 + phi); equality is ==, no tolerance",
                      "float_contrast": "binary float (64-bit), Fibonacci check only"},
            "fibonacci": fibonacci(lv["fibonacci"]),
            "pentagon": pentagon(lv["pentagon"]),
            "penrose": penrose(lv["penrose"])}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed, SMALL if small else LEVELS)
    f, g, t = out["fibonacci"], out["pentagon"], out["penrose"]
    comp = emlib.compressibility(list(_final_word(f["levels"])),
                                 "Fibonacci word at the last level, letters L / S")
    fl = [r["float_contrast_max_abs_difference"] for r in f["rows"]]
    ai = (f"OBSERVED, Fibonacci chain, {f['levels']} levels (last level "
          f"{f['rows'][-1]['n_tiles']} tiles): positions from the two constructions equal at "
          f"every level: {f['zero_drift']}; position mismatches in total "
          f"{sum(r['n_position_mismatches'] for r in f['rows'])}. Float contrast, largest "
          f"difference between the same two constructions per level: first {fl[0]:.3e}, last "
          f"{fl[-1]:.3e}, unequal positions at the last level "
          f"{f['rows'][-1]['float_contrast_n_unequal']} of {f['rows'][-1]['n_tiles'] + 1}. "
          f"OBSERVED, pentagon, {g['levels']} nested levels: every inner vertex on exactly two "
          f"diagonals, diagonal^2 / side^2 = phi^2, side^2 ratio phi^4, at every level: "
          f"{g['zero_drift']}. OBSERVED, Penrose subdivision, {t['levels']} levels (last level "
          f"{t['rows'][-1]['n_triangles']} triangles, {t['rows'][-1]['n_distinct_vertices_exact']} "
          f"distinct vertices by exact equality): area conserved exactly, counts on the "
          f"recurrence, every squared edge a power of phi, at every level: {t['zero_drift']}; "
          f"distinct squared edge lengths at the last level "
          f"{t['rows'][-1]['n_distinct_squared_edge_lengths']}. Claim logged before the "
          f"reading: {CLAIM} AI reading: not interpretable atm beyond these values.")
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "levels": out["bound"]["levels"], "claim": CLAIM,
                         "status_before_run": "PROPOSED"},
                        out, None, ai, comp, None)
    print(p)
    print(ai)


def _final_word(levels):
    w = "L"
    for _ in range(levels):
        w = I.fib_substitute(w)
    return w


if __name__ == "__main__":
    main(sys.argv[1:])
