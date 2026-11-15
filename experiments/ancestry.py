"""ah. ancestry -- do a shape's voids sit on its parent, or on its parent's dual?

CC0. Co-authored. Stdlib only. Follows experiments/dual_partner.py, where
one control (the rectified 5-cell) had its voids covered by the 5-cell.

Pattern, PROPOSED: the probe that covers a shape's voids is its dual OR its
parent.

  baselines   the rectified and the truncated forms of the 5-cell, 16-cell,
              tesseract and 24-cell, built from the parent's vertices and
              edges (rectified: edge midpoints; truncated: two points on
              every edge, at the fraction that makes the cut faces regular).
  probes      the parent, and the parent's dual, each in the frame it is
              built in (the dual of the 5-cell is its negative; the 16-cell
              and tesseract are each other's dual as built; the dual of the
              24-cell is the 16-cell and tesseract directions together).
  reading     every void the hill climb finds, grouped by depth. For each
              depth class: do the parent's directions sit on it and cover
              it? do the dual's? do the two together?

No rotation is searched: the forms share their parent's frame. Observer
slots null.

  python3 experiments/ancestry.py [seed]
  python3 experiments/ancestry.py --small [seed]   # phone bound; own run file
"""
import itertools
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from toolkit import complement as CP  # noqa: E402
from toolkit import primitives as P  # noqa: E402

NAME = "ancestry"
STARTS, STEPS = 240, 400
MATCH_DEG = 2.0
DEPTH_CLASS_DEG = 1.0
SMALL = {"starts": 60, "steps": 250, "parents": ("5-cell", "16-cell")}
PATTERN = "PROPOSED: the probe that covers a shape's voids is its dual OR its parent."


def parents():
    sim = P.simplex(4)
    cross, cube, c24 = P.cross_polytope(4), P.hypercube(4), P.cell24()
    return {"5-cell": {"vertices": sim, "dual": [[-x for x in v] for v in sim],
                       "dual_name": "5-cell, negated", "face": 3},
            "16-cell": {"vertices": cross, "dual": cube, "dual_name": "tesseract", "face": 3},
            "tesseract": {"vertices": cube, "dual": cross, "dual_name": "16-cell", "face": 4},
            "24-cell": {"vertices": c24,
                        "dual": CP.unit_dirs(cross) + CP.unit_dirs(cube),
                        "dual_name": "24-cell in dual position", "face": 3}}


def edges(verts):
    d = {(i, j): math.dist(verts[i], verts[j])
         for i, j in itertools.combinations(range(len(verts)), 2)}
    m = min(d.values())
    return [k for k, v in d.items() if v <= m * (1 + 1e-9)], m


def rectified(verts):
    es, _ = edges(verts)
    return [[(a + b) / 2 for a, b in zip(verts[i], verts[j])] for i, j in es]


def truncated(verts, face):
    """Two points on every edge, at fraction t from each end; t makes a face
    with `face` sides into a regular polygon with twice as many."""
    t = 1.0 / (2.0 + 2.0 * math.cos(math.pi / face))
    es, _ = edges(verts)
    out = []
    for i, j in es:
        a, b = verts[i], verts[j]
        out.append([x + t * (y - x) for x, y in zip(a, b)])
        out.append([y + t * (x - y) for x, y in zip(a, b)])
    return out, t


def on_one_sphere(points):
    n = [P.norm(p) for p in points]
    return (max(n) - min(n)) / max(n)


def depth_classes(voids):
    classes = []
    for v in sorted(voids, key=lambda v: -v["depth_deg"]):
        if classes and classes[-1]["depth_deg"] - v["depth_deg"] <= DEPTH_CLASS_DEG:
            classes[-1]["directions"].append(v["direction"])
        else:
            classes.append({"depth_deg": v["depth_deg"], "directions": [v["direction"]]})
    return classes


def cover(void_dirs, probe):
    """Do the probe's directions sit on these voids, and are all of them covered?"""
    pd = CP.unit_dirs(probe)
    probe_on = [math.degrees(CP.angle_to_nearest(p, void_dirs)) <= MATCH_DEG for p in pd]
    void_hit = [math.degrees(CP.angle_to_nearest(v, pd)) <= MATCH_DEG for v in void_dirs]
    return {"n_probe_directions": len(pd), "probe_directions_on_a_void": sum(probe_on),
            "n_voids": len(void_dirs), "voids_with_a_probe_direction": sum(void_hit),
            "every_void_has_a_probe_direction": all(void_hit),
            "every_probe_direction_is_on_a_void": all(probe_on),
            "sits_in_and_covers": all(void_hit) and all(probe_on)}


def read(name, form, verts, parent, seed, starts, steps):
    dirs = CP.unit_dirs(verts)
    v = CP.sphere_voids(dirs, seed, starts, steps, merge_deg=2.0, keep_deg=360.0)
    classes = depth_classes(v["voids"])
    rows = []
    for c in classes:
        cp, cd = cover(c["directions"], parent["vertices"]), cover(c["directions"], parent["dual"])
        both = cover(c["directions"],
                     CP.unit_dirs(parent["vertices"]) + CP.unit_dirs(parent["dual"]))
        rows.append({"depth_deg": c["depth_deg"], "n_voids_found": len(c["directions"]),
                     "parent": cp, "dual": cd, "parent_and_dual_together": both,
                     "covered_by": [k for k, x in (("parent", cp), ("dual", cd),
                                                   ("parent+dual", both))
                                    if x["sits_in_and_covers"]]})
    top = rows[0]
    return {"baseline": f"{form} {name}", "parent": name, "form": form,
            "n_vertices": len(dirs), "vertices_on_one_sphere_spread": on_one_sphere(verts),
            "dual_of_parent": parent["dual_name"],
            "n_void_depth_classes": len(rows), "void_classes": rows,
            "deepest_voids_covered_by": top["covered_by"],
            "parent_covers_the_deepest_voids": top["parent"]["sits_in_and_covers"],
            "dual_covers_the_deepest_voids": top["dual"]["sits_in_and_covers"],
            "classes_covered_by_parent": [i for i, r in enumerate(rows)
                                          if "parent" in r["covered_by"]],
            "classes_covered_by_dual": [i for i, r in enumerate(rows) if "dual" in r["covered_by"]],
            "classes_covered_by_none": [i for i, r in enumerate(rows) if not r["covered_by"]],
            "bound": v["bound"]}


def run(seed=0, starts=STARTS, steps=STEPS, which=None):
    ps = parents()
    rows = []
    for name, par in ps.items():
        if which and name not in which:
            continue
        # the parent itself, as a reference: its deepest voids against its dual
        rows.append(read(name, "plain", par["vertices"], par, seed, starts, steps))
        rows.append(read(name, "rectified", rectified(par["vertices"]), par, seed, starts, steps))
        tv, t = truncated(par["vertices"], par["face"])
        r = read(name, "truncated", tv, par, seed, starts, steps)
        r["truncation_fraction"] = t
        rows.append(r)
    forms = [r for r in rows if r["form"] != "plain"]
    return {"pattern": PATTERN,
            "bound": {"seed": seed, "starts": starts, "steps": steps, "match_deg": MATCH_DEG,
                      "depth_class_deg": DEPTH_CLASS_DEG,
                      "frame": "every form in its parent's frame; no rotation searched",
                      "voids_found_is_a_lower_bound": "a void no start climbed to is not "
                                                      "listed; a class may be incomplete",
                      "covers_means": "every void of the class has a probe direction within "
                                      f"{MATCH_DEG} degrees AND every probe direction is on a "
                                      "void of the class",
                      "duals_are": "standard, carried; not derived here"},
            "baselines": rows,
            "summary": {
                "forms": len(forms),
                "deepest_covered_by_parent": [r["baseline"] for r in forms
                                              if r["parent_covers_the_deepest_voids"]],
                "deepest_covered_by_dual": [r["baseline"] for r in forms
                                            if r["dual_covers_the_deepest_voids"]],
                "deepest_covered_by_parent_and_dual_together":
                    [r["baseline"] for r in forms
                     if "parent+dual" in r["deepest_voids_covered_by"]],
                "deepest_covered_by_none_of_the_three": [r["baseline"] for r in forms
                                                         if not r["deepest_voids_covered_by"]],
                "some_class_covered_by_parent": [r["baseline"] for r in forms
                                                 if r["classes_covered_by_parent"]],
                "some_class_covered_by_dual": [r["baseline"] for r in forms
                                               if r["classes_covered_by_dual"]],
                "plain_deepest_covered_by": {r["baseline"]: r["deepest_voids_covered_by"]
                                             for r in rows if r["form"] == "plain"}}}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["starts"], SMALL["steps"], SMALL["parents"]) if small else run(seed))
    symbols = [("p" if r["parent_covers_the_deepest_voids"] else "")
               + ("d" if r["dual_covers_the_deepest_voids"] else "") or "n"
               for r in out["baselines"]]
    comp = emlib.compressibility(
        symbols, "per baseline, who covers its deepest voids: p parent, d dual, pd both, n "
                 "neither")
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "starts": out["bound"]["starts"],
                         "steps": out["bound"]["steps"],
                         "parents": sorted({r["parent"] for r in out["baselines"]}),
                         "forms": ["plain", "rectified", "truncated"],
                         "pattern": PATTERN, "status_before_run": "PROPOSED"},
                        out, None, None, comp, None, gz=True)
    print(p)
    for r in out["baselines"]:
        print(f"{r['baseline']:22s} n={r['n_vertices']:3d} | "
              + "  ".join(f"[{c['depth_deg']:.1f} deg x{c['n_voids_found']}: "
                          f"{'+'.join(c['covered_by']) or 'neither'}]"
                          for c in r["void_classes"][:4]))
    print("SUMMARY", out["summary"])


if __name__ == "__main__":
    main(sys.argv[1:])
