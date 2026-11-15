"""af. dual_partner -- is a shape's best-separated partner its dual?

CC0. Co-authored. Stdlib only. Follows experiments/sphere_compounds.py.

Pattern, PROPOSED: among a set of probe shapes, the one that can be placed
farthest from a baseline (largest smallest-angle between the two) is the
baseline's dual; it sits in the baseline's voids.

  sphere   baselines 5-cell, 16-cell, tesseract, 24-cell; probes the same
           four. For every pair the probe is turned in R^4 to sit as far
           from the baseline as the search can put it. Logged per baseline:
           every probe's separation, the rank of the dual (ties shown), the
           deepest void left, and whether the probe's directions sit on the
           baseline's deepest voids and cover them.
  plane    regular polygons on one circle; a k-gon's dual is the k-gon
           turned by 180/k. Separation over all turns, by scanning.
  control  baselines whose dual is NOT in the probe set: a random set, the
           rectified 5-cell, a 5-cell with its negative (sphere); a 9-gon
           and a 10-gon (plane, probes without them).

Duals used (standard, not derived here): the 5-cell and the 24-cell are
self-dual; the 16-cell and the tesseract are each other's dual; a regular
polygon is self-dual.

Observer slots null.

  python3 experiments/dual_partner.py [seed]
  python3 experiments/dual_partner.py --small [seed]   # phone bound; own run file
"""
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import complement_fit_v2 as V  # noqa: E402
import emlib  # noqa: E402
import sphere_compounds as SC  # noqa: E402
from toolkit import complement as CP  # noqa: E402
from toolkit import primitives as P  # noqa: E402

NAME = "dual_partner"
PROBES = ("5-cell", "16-cell", "tesseract", "24-cell")
DUAL = {"5-cell": "5-cell", "16-cell": "tesseract", "tesseract": "16-cell", "24-cell": "24-cell"}
CONTROLS = ("random_8", "rectified 5-cell", "5-cell + its negative (as built)")
POLYGONS = (3, 4, 5, 6, 7, 8, 12)
POLYGON_CONTROLS = (9, 10)
STARTS, POLISH = 12, 4000
VOID_STARTS, VOID_STEPS = 200, 400
TIE_DEG = 0.05               # separations closer than this are a tie
SMALL = {"starts": 3, "polish": 800, "void_starts": 40, "void_steps": 250,
         "probes": ("5-cell", "16-cell", "tesseract"), "controls": ("random_8",)}
PATTERN = ("PROPOSED: the probe that can be placed farthest from a baseline is the baseline's "
           "dual, and it sits in the baseline's voids.")


def rank_with_ties(rows, key, reverse=True):
    """Competition ranks: values within TIE_DEG share a rank."""
    order = sorted(rows, key=lambda r: -r[key] if reverse else r[key])
    rank, last = 0, None
    for i, r in enumerate(order):
        if last is None or abs(r[key] - last) > TIE_DEG:
            rank, last = i + 1, r[key]
        r[key + "_rank"] = rank
    return order


def sphere_part(seed, probes, controls, starts, polish_moves, void_starts, void_steps):
    cat = V.catalogue()
    rng0 = random.Random(f"{seed}/dual-partner-control")
    shapes = {k: cat[k] for k in cat}
    shapes["random_8"] = CP.unit_dirs([[rng0.gauss(0, 1) for _ in range(4)] for _ in range(8)])
    cache = {}

    def best(a, b):
        key = tuple(sorted((a, b)))
        if key not in cache:
            x, y = shapes[key[0]], shapes[key[1]]
            rng = random.Random(f"{seed}/{key[0]}/{key[1]}")
            n_starts, n_jig = V.search_effort(len(x), len(y), starts)
            s = V.search_rotation(x, y, rng, n_starts, jiggles=n_jig)
            ang, val = SC.polish(x, y, s["angles_deg"], rng, polish_moves)
            cache[key] = {"first": key[0], "angles_deg": ang, "separation_deg": val,
                          "starts": n_starts, "random_moves": n_jig, "polish_moves": polish_moves}
        return cache[key]

    out = []
    for b in list(probes) + list(controls):
        base = shapes[b]
        voids = CP.sphere_voids(base, seed, void_starts, void_steps)
        vdirs = [v["direction"] for v in voids["voids"]]
        rows = []
        for p in probes:
            c = best(b, p)
            # the search turns the second of the sorted pair; put the baseline first
            if c["first"] == b:
                turned = V.rotate4(shapes[p], c["angles_deg"])
                base_used = base
            else:
                turned = shapes[p]
                base_used = V.rotate4(base, c["angles_deg"])
            v_used = CP.sphere_voids(base_used, seed, void_starts, void_steps)
            vd = [v["direction"] for v in v_used["voids"]]
            m = CP.match_directions(turned, vd) if vd else None
            cover = CP.match_directions(vd, turned)
            union = CP.sphere_voids(base_used + turned, seed, max(40, void_starts // 3),
                                    void_steps)
            rows.append({"probe": p, "n_points": len(shapes[p]),
                         "separation_deg": c["separation_deg"],
                         "residual_gap_deg": union["deepest_depth_deg"],
                         "probe_points_on_baseline_voids":
                             (m["max_angle_deg"] <= 2.0) if m else None,
                         "largest_angle_probe_point_to_void_deg": m["max_angle_deg"] if m else None,
                         "baseline_voids_covered": cover["candidates_hit"] if False else
                         sum(1 for vv in vd if math.degrees(CP.angle_to_nearest(vv, turned)) <= 2.0),
                         "n_baseline_voids": len(vd),
                         "sits_in_and_covers_the_voids":
                             bool(m) and m["max_angle_deg"] <= 2.0
                             and all(math.degrees(CP.angle_to_nearest(vv, turned)) <= 2.0
                                     for vv in vd),
                         "search": {k: c[k] for k in ("starts", "random_moves", "polish_moves")}})
        rank_with_ties(rows, "separation_deg")
        rank_with_ties(rows, "residual_gap_deg", reverse=False)
        rows.sort(key=lambda r: (r["separation_deg_rank"], r["probe"]))
        dual = DUAL.get(b)
        drow = next((r for r in rows if r["probe"] == dual), None)
        top = [r["probe"] for r in rows if r["separation_deg_rank"] == 1]
        out.append({"baseline": b, "n_points": len(base), "control": b in controls,
                    "dual": dual, "dual_in_probe_set": drow is not None,
                    "baseline_deepest_void_deg": voids["deepest_depth_deg"],
                    "n_deepest_voids_found": len(vdirs),
                    "probes": rows,
                    "best_separated": top, "best_separated_is_a_tie": len(top) > 1,
                    "dual_is_best_separated": (dual in top) if drow else None,
                    "dual_is_the_only_best": (top == [dual]) if drow else None,
                    "rank_of_dual_by_separation": drow["separation_deg_rank"] if drow else None,
                    "rank_of_dual_by_residual_gap": drow["residual_gap_deg_rank"] if drow else None,
                    "probes_sitting_in_and_covering_the_voids":
                        [r["probe"] for r in rows if r["sits_in_and_covers_the_voids"]],
                    "dual_sits_in_and_covers_the_voids":
                        drow["sits_in_and_covers_the_voids"] if drow else None})
    return {"probes": list(probes), "controls": list(controls), "tie_deg": TIE_DEG,
            "baselines": out}


def polygon_separation(k, m, step_deg=0.05):
    """Largest, over turns of the m-gon, of the smallest angle between a vertex
    of the k-gon and a vertex of the m-gon (same circle). Scanned in steps."""
    a = [360.0 * i / k for i in range(k)]
    best, best_turn = -1.0, 0.0
    period = 360.0 / m
    steps = int(round(period / step_deg))
    for s in range(steps + 1):
        t = s * step_deg
        b = [(t + 360.0 * j / m) % 360.0 for j in range(m)]
        sep = min(min(abs(x - y), 360.0 - abs(x - y)) for x in a for y in b)
        if sep > best + 1e-12:
            best, best_turn = sep, t
    return best, best_turn


def plane_part(polygons, controls):
    out = []
    for k in list(polygons) + list(controls):
        rows = []
        for m in polygons:
            sep, turn = polygon_separation(k, m)
            rows.append({"probe": f"{m}-gon", "n_points": m, "separation_deg": sep,
                         "at_turn_deg": turn})
        rank_with_ties(rows, "separation_deg")
        rows.sort(key=lambda r: (r["separation_deg_rank"], r["n_points"]))
        dual = f"{k}-gon"
        drow = next((r for r in rows if r["probe"] == dual), None)
        top = [r["probe"] for r in rows if r["separation_deg_rank"] == 1]
        out.append({"baseline": f"{k}-gon", "control": k in controls,
                    "dual": dual + " turned by 180/k", "dual_in_probe_set": drow is not None,
                    "half_gap_deg": 180.0 / k,
                    "probes": rows, "best_separated": top,
                    "best_separated_is_a_tie": len(top) > 1,
                    "dual_is_best_separated": (dual in top) if drow else None,
                    "dual_is_the_only_best": (top == [dual]) if drow else None,
                    "rank_of_dual_by_separation": drow["separation_deg_rank"] if drow else None,
                    "dual_turn_deg": drow["at_turn_deg"] if drow else None})
    return {"probes": [f"{m}-gon" for m in polygons], "controls": [f"{k}-gon" for k in controls],
            "scan_step_deg": 0.05, "tie_deg": TIE_DEG, "baselines": out}


def run(seed=0, probes=PROBES, controls=CONTROLS, starts=STARTS, polish_moves=POLISH,
        void_starts=VOID_STARTS, void_steps=VOID_STEPS):
    sph = sphere_part(seed, probes, controls, starts, polish_moves, void_starts, void_steps)
    pl = plane_part(POLYGONS, POLYGON_CONTROLS)
    real = [b for b in sph["baselines"] if not b["control"]]
    preal = [b for b in pl["baselines"] if not b["control"]]
    return {"pattern": PATTERN,
            "bound": {"seed": seed, "starts": starts, "polish_moves": polish_moves,
                      "void_starts": void_starts, "void_steps": void_steps,
                      "tie_deg": TIE_DEG, "void_match_deg": 2.0,
                      "search_is_a_lower_bound": "a probe's separation may be larger than "
                                                 "found; ranks can change with more search",
                      "more_points_sit_closer": "a probe with more points cannot be placed as "
                                                "far; separation is not scaled for that",
                      "duals_are": "standard, carried; not derived here",
                      "probe_sets_are_short": "four polytopes; seven polygons"},
            "sphere": sph, "plane": pl,
            "summary": {
                "sphere_baselines": len(real),
                "sphere_dual_best_separated": sum(1 for b in real if b["dual_is_best_separated"]),
                "sphere_dual_the_only_best": sum(1 for b in real if b["dual_is_the_only_best"]),
                "sphere_dual_sits_in_and_covers_voids":
                    sum(1 for b in real if b["dual_sits_in_and_covers_the_voids"]),
                "sphere_non_dual_probes_that_also_cover":
                    {b["baseline"]: [p for p in b["probes_sitting_in_and_covering_the_voids"]
                                     if p != b["dual"]] for b in real},
                "sphere_controls_with_a_covering_probe":
                    {b["baseline"]: b["probes_sitting_in_and_covering_the_voids"]
                     for b in sph["baselines"] if b["control"]},
                "plane_baselines": len(preal),
                "plane_dual_best_separated": sum(1 for b in preal if b["dual_is_best_separated"]),
                "plane_dual_the_only_best": sum(1 for b in preal if b["dual_is_the_only_best"])}}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["probes"], SMALL["controls"], SMALL["starts"], SMALL["polish"],
               SMALL["void_starts"], SMALL["void_steps"]) if small else run(seed))
    symbols = [str(b["rank_of_dual_by_separation"] or 0)
               for b in out["sphere"]["baselines"] + out["plane"]["baselines"]]
    comp = emlib.compressibility(
        symbols, "rank of the dual by separation, per baseline (0 = dual not in the probe "
                 "set): sphere baselines, then plane")
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "probes": out["sphere"]["probes"],
                         "controls": out["sphere"]["controls"], "polygons": list(POLYGONS),
                         "polygon_controls": list(POLYGON_CONTROLS),
                         "starts": out["bound"]["starts"],
                         "polish_moves": out["bound"]["polish_moves"],
                         "pattern": PATTERN, "status_before_run": "PROPOSED"},
                        out, None, None, comp, None, gz=True)
    print(p)
    print("SPHERE")
    for b in out["sphere"]["baselines"]:
        print(f"  {b['baseline'][:28]:28s} dual {str(b['dual']):9s} | "
              + "  ".join(f"{r['probe']} {r['separation_deg']:.2f} (gap {r['residual_gap_deg']:.1f}"
                          f"{', covers' if r['sits_in_and_covers_the_voids'] else ''})"
                          for r in b["probes"])
              + f" | dual rank {b['rank_of_dual_by_separation']}, best {b['best_separated']}")
    print("PLANE")
    for b in out["plane"]["baselines"]:
        print(f"  {b['baseline']:7s} | "
              + "  ".join(f"{r['probe']} {r['separation_deg']:.2f}" for r in b["probes"])
              + f" | dual rank {b['rank_of_dual_by_separation']}, best {b['best_separated']}")
    print("SUMMARY", out["summary"])


if __name__ == "__main__":
    main(sys.argv[1:])
