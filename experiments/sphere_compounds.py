"""ae. sphere_compounds -- pairs of 4-polytopes on the sphere, checked closely.

CC0. Co-authored. Stdlib only. Follows experiments/complement_fit_v2.py.

  1 two 16-cells   the rotation search again with more refinement, target
                   separation 60 degrees; and the placement built outright
                   in exact fractions: the tesseract's sixteen vertices
                   split by sign parity into two sets of eight.
  2 catalogue      two names CARRIED from a co-author, with what could be
                   checked here: "compound of two 5-cells" (a 5-cell and
                   its negative, 10 directions) and "disphenoidal 288-cell"
                   (48 directions: a 24-cell and its dual). Checked by
                   computing each set's convex hull: vertices, edges, cells,
                   and the shape of the cells. The names themselves were not
                   checked against a source: the pages were out of reach.
  3 two 24-cells   the dual placement built outright, set beside the
                   searched placement of complement_fit_v2 and beside the
                   48-direction catalogue entry.

Observer slots null.

  python3 experiments/sphere_compounds.py [seed]
  python3 experiments/sphere_compounds.py --small [seed]   # phone bound; own run file
"""
import itertools
import math
import random
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import complement_fit_v2 as V  # noqa: E402
import emlib  # noqa: E402
from toolkit import complement as CP  # noqa: E402
from toolkit import hull4 as H4  # noqa: E402
from toolkit import primitives as P  # noqa: E402

NAME = "sphere_compounds"
STARTS, JIGGLES, POLISH = 16, 2500, 20000
SMALL = {"starts": 4, "jiggles": 600, "polish": 3000}
FITV2 = "runs/complement_fit_v2_seed{seed}.json"
ASKED = ("PROPOSED: two 16-cells reach a separation of 60 degrees and join into a tesseract; "
         "the 48 directions of a 24-cell and its dual are the disphenoidal 288-cell (CARRIED); "
         "a 5-cell and its negative are the compound of two 5-cells (CARRIED).")


def polish(base, probe, angles, rng, moves):
    """Many small random moves of all six angles, shrinking from 0.5 degrees."""
    ang = list(angles)
    val = V.separation(base, V.rotate4(probe, ang))
    size = 0.5
    for _ in range(moves):
        t = [a + rng.gauss(0.0, size) for a in ang]
        v = V.separation(base, V.rotate4(probe, t))
        if v > val + 1e-13:
            ang, val = t, v
        else:
            size = max(size * 0.999, 1e-8)
    return ang, val


def against(sig, cat_sigs):
    rows = []
    for k, sg in cat_sigs.items():
        if len(sg) == len(sig):
            rows.append({"name": k,
                         "largest_angle_difference_deg": max(abs(x - y) for x, y in zip(sig, sg)),
                         "match_at_0.5_deg": V.same_signature(sig, sg)})
    return sorted(rows, key=lambda r: r["largest_angle_difference_deg"])


def two_16_cells(seed, starts, jiggles, polish_moves, cat, sigs):
    base = cat["16-cell"]
    rng = random.Random(f"{seed}/two-16-cells")
    s = V.search_rotation(base, base, rng, starts, jiggles=jiggles)
    ang, val = polish(base, base, s["angles_deg"], rng, polish_moves)
    union = base + V.rotate4(base, ang)
    sig = V.signature(union)
    # the placement built outright, in exact fractions
    h = Fraction(1, 2)
    verts = [tuple(h * sg for sg in signs) for signs in itertools.product((1, -1), repeat=4)]
    even = [v for v in verts if sum(1 for x in v if x < 0) % 2 == 0]
    odd = [v for v in verts if sum(1 for x in v if x < 0) % 2 == 1]

    def dots(a, b, same):
        return sorted({sum(x * y for x, y in zip(p, q)) for i, p in enumerate(a)
                       for j, q in enumerate(b) if not (same and i >= j)})
    exact = {"arithmetic": "exact fractions",
             "set_sizes": [len(even), len(odd)],
             "dots_within_even_set": [str(x) for x in dots(even, even, True)],
             "dots_within_odd_set": [str(x) for x in dots(odd, odd, True)],
             "each_set_is_a_16_cell": dots(even, even, True) == [Fraction(-1), Fraction(0)]
             and dots(odd, odd, True) == [Fraction(-1), Fraction(0)],
             "dots_between_the_sets": [str(x) for x in dots(even, odd, False)],
             "smallest_angle_between_the_sets": "60 degrees exactly (largest dot 1/2)"
             if max(sum(x * y for x, y in zip(p, q)) for p in even for q in odd) == h else None,
             "union_is_the_tesseract_vertex_set": set(even) | set(odd) == set(verts)}
    built = [[float(x) for x in v] for v in even + odd]
    return {"search": {"starts": starts, "random_moves": jiggles, "polish_moves": polish_moves,
                       "separation_after_search_deg": s["separation_deg"],
                       "separation_after_polish_deg": val, "target_deg": 60.0,
                       "short_of_target_deg": 60.0 - val, "rotation_angles_deg": ang,
                       "against_catalogue": against(sig, sigs),
                       "matches_tesseract_at_0.5_deg":
                           V.same_signature(sig, sigs["tesseract"]),
                       "distinct_pairwise_angles_deg": CP.distinct_values(sig, tol=0.05)},
            "built_outright": dict(exact,
                                   separation_deg=V.separation(CP.unit_dirs(built[:8]),
                                                               CP.unit_dirs(built[8:])),
                                   against_catalogue=against(V.signature(CP.unit_dirs(built)),
                                                             sigs),
                                   hull=H4.summary(CP.unit_dirs(built)))}


def hull_entry(name, dirs, carried, max_edge_deg=None, hull_name=None):
    s = H4.summary(dirs, max_edge_deg)
    v, c = s["n_points"], s["n_facets"]
    e = s["n_edges_if_all_tetrahedra"]
    faces = (e - v + c) if e is not None else None          # Euler: V - E + F - C = 0
    pats = s["edge_angle_patterns"]
    disph = None
    if len(pats) == 1 and len(pats[0]["angles_deg"]) == 6:
        a = pats[0]["angles_deg"]
        vals = sorted(set(a))
        disph = (len(vals) == 2 and sorted(a.count(x) for x in vals) == [2, 4])
    profile = set()
    for i, p in enumerate(dirs):
        profile.add(tuple(round(math.degrees(math.acos(max(-1.0, min(1.0, P.dot(p, q))))), 3)
                          for q in sorted((q for j, q in enumerate(dirs) if j != i),
                                          key=lambda q: -P.dot(p, q))))
    return {"name_carried": name, "name_checked_against_a_source": False,
            # the name carried for the convex hull of these directions, when it is
            # not the name of the set itself (a compound is not its hull)
            "hull_name_carried": hull_name,
            "hull_name_checked_against_a_source": False if hull_name else None,
            "why_not": "the reference pages were outside the allowed network",
            "carried_statement": carried,
            "computed": {"vertices": v, "edges": e, "faces_by_euler": faces, "cells": c,
                         "cells_by_vertex_count": s["facets_by_vertex_count"],
                         "cell_edge_angle_patterns": pats,
                         "all_cells_one_shape": len(pats) == 1,
                         "cells_are_tetragonal_disphenoids": disph,
                         "cells_close_up": s["tetrahedral_facets_close_up"],
                         "every_vertex_sees_the_same_angles": len(profile) == 1,
                         "distinct_pairwise_angles_deg":
                             CP.distinct_values(V.signature(dirs), tol=0.05),
                         "bound": s["bound"]}}


def two_24_cells(seed, cat, sigs):
    c24 = cat["24-cell"]
    dual = cat["16-cell + tesseract (as built)"]
    union = c24 + dual
    sig = V.signature(union)
    out = {"built_outright": {
        "placement": "the 24-cell as built with the 16-cell and tesseract directions as built",
        "dual_is_a_24_cell": V.same_signature(V.signature(dual), sigs["24-cell"]),
        "separation_deg": V.separation(c24, dual),
        "deepest_void": CP.sphere_voids(union, seed, 60, 300)["deepest_depth_deg"],
        "against_catalogue": against(sig, sigs),
        "distinct_pairwise_angles_deg": CP.distinct_values(sig, tol=0.05)}}
    try:
        rec = emlib.load_run(emlib.ROOT / FITV2.format(seed=seed))
        row = next(r for r in rec["outputs"]["sphere"]["pairs"]
                   if r["baseline"] == "24-cell" and r["probe"] == "24-cell")
        ssig = V.signature(row["joined_directions"])
        out["searched_in_complement_fit_v2"] = {
            "separation_deg": row["separation_after_search_deg"],
            "residual_gap_deg": row["residual_gap_after_search_deg"],
            "largest_angle_difference_from_the_built_placement_deg":
                max(abs(x - y) for x, y in zip(ssig, sig)),
            "same_as_the_built_placement_at_0.5_deg": V.same_signature(ssig, sig),
            "hull": H4.summary(row["joined_directions"], 61.0)}
    except FileNotFoundError:
        out["searched_in_complement_fit_v2"] = None
    return out


def run(seed=0, starts=STARTS, jiggles=JIGGLES, polish_moves=POLISH):
    cat = V.catalogue()
    sigs = {k: V.signature(v) for k, v in cat.items()}
    known = {n: H4.summary(cat[n]) for n in ("5-cell", "16-cell", "tesseract", "24-cell")}
    return {"asked": ASKED,
            "bound": {"seed": seed, "starts": starts, "random_moves": jiggles,
                      "polish_moves": polish_moves, "match_tol_deg": V.MATCH_TOL_DEG,
                      "hull": "brute force over 4-subsets; for the 48-direction set only "
                              "subsets with every pair within 61 degrees, checked to close up",
                      "names": "carried; not checked against a source here",
                      "search_is_a_lower_bound": True},
            "hull_self_check_on_known_polytopes":
                {n: {"cells": s["n_facets"], "cells_by_vertex_count": s["facets_by_vertex_count"]}
                 for n, s in known.items()},
            "two_16_cells": two_16_cells(seed, starts, jiggles, polish_moves, cat, sigs),
            "catalogue_additions": {
                "compound of two 5-cells": hull_entry(
                    "compound of two 5-cells", cat["5-cell + its negative (as built)"],
                    "a 5-cell together with its negative (dual position); 10 vertices",
                    hull_name="disphenoidal 30-cell (dual of bitruncated 5-cell)"),
                "disphenoidal 288-cell": hull_entry(
                    "disphenoidal 288-cell", cat["24-cell + its dual (as built)"],
                    "48 vertices = a 24-cell and its dual 24-cell; 288 cells", 61.0)},
            "two_24_cells": two_24_cells(seed, cat, sigs)}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["starts"], SMALL["jiggles"], SMALL["polish"]) if small else run(seed))
    t = out["two_16_cells"]
    comp = emlib.compressibility(
        [str(len(str(c["computed"]["cells"]))) for c in out["catalogue_additions"].values()],
        "digits in the cell count of each catalogue addition, in order")
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "starts": out["bound"]["starts"],
                         "random_moves": out["bound"]["random_moves"],
                         "polish_moves": out["bound"]["polish_moves"],
                         "fit_v2_reference": FITV2.format(seed=seed), "asked": ASKED,
                         "status_before_run": "PROPOSED"},
                        out, None, None, comp, None, gz=True)
    print(p)
    print("hull self-check:", out["hull_self_check_on_known_polytopes"])
    s = t["search"]
    print(f"two 16-cells, searched: {s['separation_after_search_deg']:.4f} -> "
          f"{s['separation_after_polish_deg']:.6f} (target 60), matches tesseract "
          f"{s['matches_tesseract_at_0.5_deg']}, against {s['against_catalogue']}")
    b = t["built_outright"]
    print("two 16-cells, built:", {k: b[k] for k in ("each_set_is_a_16_cell",
                                                     "dots_between_the_sets",
                                                     "union_is_the_tesseract_vertex_set",
                                                     "separation_deg")},
          "against", b["against_catalogue"], "hull", b["hull"]["facets_by_vertex_count"])
    for k, c in out["catalogue_additions"].items():
        print(k, "| hull name carried:", c["hull_name_carried"], {x: c["computed"][x] for x in ("vertices", "edges", "faces_by_euler", "cells",
                                                "all_cells_one_shape",
                                                "cells_are_tetragonal_disphenoids",
                                                "cells_close_up",
                                                "every_vertex_sees_the_same_angles")},
              c["computed"]["cell_edge_angle_patterns"])
    w = out["two_24_cells"]
    print("two 24-cells, built:", {k: w["built_outright"][k] for k in
                                   ("dual_is_a_24_cell", "separation_deg", "deepest_void")},
          w["built_outright"]["against_catalogue"])
    if w["searched_in_complement_fit_v2"]:
        x = w["searched_in_complement_fit_v2"]
        print("two 24-cells, searched:", {k: x[k] for k in x if k != "hull"},
              "hull cells", x["hull"]["n_facets"], x["hull"]["facets_by_vertex_count"])


if __name__ == "__main__":
    main(sys.argv[1:])
