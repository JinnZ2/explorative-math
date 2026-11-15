"""q. degeneracy_check -- coincident or collapsed points in every final projection.

CC0. Co-authored. Stdlib only.

A projection in which two points land on each other, or all points land on
a line, has mirrors the polytope's own symmetry did not give it. This run
reads every final projection logged by the pursuit runs, the Coxeter-plane
run and frame 0 of the tour with probes.degeneracy (smallest pairwise
distance, collinearity), tags the ones where points coincide or collapse,
and recounts the starts reaching h mirrors with those excluded.

No run file is changed; the rows go to this experiment's own run file. A
source run that is absent is listed as absent, not skipped silently.

  python3 experiments/degeneracy_check.py [seed]
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "degeneracy_check"
TOLS = (1e-6, 1e-3)          # degeneracy tags are read at both
TAG_TOL = "1e-06"            # the recount excludes projections degenerate at this one
COUNT_TOL = 1e-6             # mirror counts are the ones logged at this tolerance
SOURCES = ("mirror_pursuit_seed0.json", "layered_pursuit_seed0.json",
           "anchor_generate_seed0.json", "constrained_generate_seed0.json",
           "coxeter_plane_seed0.json", "mirror_sweep_seed0.json",
           "fold_loop_seed0.json", "fold_shape_seed0.json", "combine_seed0.json")
MAKE = {"5-cell": lambda: P.simplex(4), "16-cell": lambda: P.cross_polytope(4)}


def row(run, shape, where, pts, n_axes, h, method=None, start=None, recomputed=False):
    d = Q.degeneracy(pts, TOLS)
    return {"run": run, "shape": shape, "where": where, "method": method, "start_index": start,
            "h": h, "n_mirror_axes": n_axes, "reaches_h": (n_axes == h) if h else None,
            "projection_recomputed_from_angles": recomputed,
            "min_pair_distance": d["min_pair_distance"], "closest_pair": d["closest_pair"],
            "collinearity": d["collinearity"], "tags": d["tags"],
            "degenerate": d["tags"][TAG_TOL]["degenerate"],
            # a point at the centre is fixed by every mirror; logged beside the
            # degenerate tag, not folded into it
            "min_radius": d["min_radius"],
            "central": d["tags"][TAG_TOL]["central"],
            "n_at_centre": d["tags"][TAG_TOL]["n_at_centre"],
            "n_points": d["n_points"],
            "more_axes_than_off_centre_points":
                (n_axes is not None
                 and n_axes > d["n_points"] - d["tags"][TAG_TOL]["n_at_centre"])}


def axes_at(counts):
    return next(c["n_mirror_axes"] for c in counts if c["rel_tol"] == COUNT_TOL)


def from_angles(shape, angles_deg):
    verts = MAKE[shape]()
    planes = P.all_planes(len(verts[0]))
    return C.rotate_project(verts, [(pl, math.radians(a)) for pl, a in zip(planes, angles_deg)])


def rows_of(name, rec):
    out = rec["outputs"]
    rows = []
    if name.startswith("mirror_pursuit"):
        for s in out["shapes"]:
            h = s["coxeter_plane_mirror_count"]
            for r in s["starts"]:
                rows.append(row(name, s["shape"], f"starts[{r['start_index']}]",
                                r["final_projected"], axes_at(r["final_counts"]), h,
                                "stage_1_smallest_residual", r["start_index"]))
    elif name.startswith("layered_pursuit"):
        for s in out["shapes"]:
            for r in s["starts"]:
                i = r["start_index"]
                for st in r["layered_stages"]:
                    last = st["k"] == s["h"]
                    rows.append(row(name, s["shape"], f"starts[{i}].layered_stages[k={st['k']}]",
                                    st["final_projected"], axes_at(st["counts"]), s["h"],
                                    "layered_unanchored" if last else
                                    "layered_unanchored_earlier_stage", i))
                rows.append(row(name, s["shape"], f"starts[{i}].alt", r["alt"]["final_projected"],
                                axes_at(r["alt"]["counts"]), s["h"], "alt_k_equals_h", i))
    elif name.startswith(("anchor_generate", "constrained_generate")):
        base = "anchor_generate" if name.startswith("anchor") else "constrained_generate"
        for s in out["shapes"]:
            for r in s["starts"]:
                i = r["start_index"]
                for label, g in r["rules"].items():
                    rows.append(row(name, s["shape"], f"starts[{i}].rules.{label}",
                                    g["final_projected"], axes_at(g["counts"]), s["h"],
                                    f"{base}_{label}", i))
                    b = g.get("bound_check")
                    if b:
                        rows.append(row(name, s["shape"], f"starts[{i}].rules.{label}.bound_check",
                                        from_angles(s["shape"], b["final_angles_deg"]),
                                        axes_at(b["counts"]), s["h"],
                                        f"{base}_{label}_bound_check", i, recomputed=True))
    elif name.startswith(("fold_loop", "fold_shape", "combine")):
        base = name.split("_seed")[0]
        for s in out["shapes"]:
            for r in s["starts"]:
                i = r["start_index"]
                for v, g in r["variants"].items():
                    rows.append(row(name, s["shape"], f"starts[{i}].variants.{v}",
                                    g["final_projected"], g["final_count"], s["h"],
                                    f"{base}_{v}", i))
    elif name.startswith("coxeter_plane"):
        for c in out["cases"]:
            rows.append(row(name, c["name"], f"cases[{c['name']}]", c["projected"],
                            c["readings"][0]["n_mirror_axes"], c["coxeter_number_h"],
                            "coxeter_plane_by_construction"))
    elif name.startswith("mirror_sweep"):
        maxes = {s["shape"]: s["max_axis_count"]["shape"] for s in out["sweep"]}
        for f in out["frame_zero"]:
            rows.append(row(name, f["shape"], f"frame_zero[{f['shape']}]", f["projected_points"],
                            maxes.get(f["shape"]), None, "tour_frame_0"))
    return rows


def recount(rows):
    """Per (run, shape, method): starts reaching h, how many of those are
    degenerate, and the count with them excluded."""
    groups = {}
    for r in rows:
        if r["h"] is None:
            continue
        groups.setdefault((r["run"], r["shape"], r["method"]), []).append(r)
    out = []
    for (run, shape, method), rs in groups.items():
        reach = [r for r in rs if r["reaches_h"]]
        out.append({"run": run, "shape": shape, "method": method, "h": rs[0]["h"],
                    "n_projections": len(rs),
                    "n_degenerate": sum(1 for r in rs if r["degenerate"]),
                    "final_counts": [r["n_mirror_axes"] for r in rs],
                    "degenerate_flags": [r["degenerate"] for r in rs],
                    "final_counts_excluding_degenerate":
                        [None if r["degenerate"] else r["n_mirror_axes"] for r in rs],
                    "reaching_h": len(reach),
                    "reaching_h_degenerate": sum(1 for r in reach if r["degenerate"]),
                    "reaching_h_excluding_degenerate":
                        sum(1 for r in reach if not r["degenerate"]),
                    "n_central": sum(1 for r in rs if r["central"]),
                    "central_flags": [r["central"] for r in rs],
                    "reaching_h_central": sum(1 for r in reach if r["central"]),
                    "largest_count_excluding_degenerate":
                        max((r["n_mirror_axes"] for r in rs if not r["degenerate"]),
                            default=None)})
    return out


def run(seed=0, runs_dir=None, sources=SOURCES):
    d = Path(runs_dir) if runs_dir else emlib.RUNS
    rows, absent = [], []
    for name in sources:
        try:
            rec = emlib.load_run(d / name)
        except FileNotFoundError:
            absent.append(name)
            continue
        rows.extend(rows_of(name, rec))
    tagged = [r for r in rows if r["degenerate"]]
    near = [r for r in rows if r["tags"]["0.001"]["degenerate"] and not r["degenerate"]]
    return {"bound": {"tag_tol": float(TAG_TOL), "tols": list(TOLS), "count_rel_tol": COUNT_TOL,
                      "normalised_by": "largest distance from the centroid",
                      "collinearity": "sqrt(smaller / larger eigenvalue of the scatter matrix)"},
            "sources_read": [s for s in sources if s not in absent], "sources_absent": absent,
            "n_rows": len(rows), "n_degenerate": len(tagged),
            "n_degenerate_only_at_1e-3": len(near),
            "n_central": sum(1 for r in rows if r["central"]),
            "n_central_not_degenerate": sum(1 for r in rows if r["central"]
                                            and not r["degenerate"]),
            "central_not_degenerate": [
                {k: r[k] for k in ("run", "shape", "where", "method", "n_mirror_axes",
                                   "n_at_centre", "n_points", "min_radius")}
                for r in rows if r["central"] and not r["degenerate"]],
            "degenerate": [{k: r[k] for k in ("run", "shape", "where", "method", "n_mirror_axes",
                                              "min_pair_distance", "collinearity")}
                           for r in tagged],
            "recount": recount(rows), "rows": rows}


def main(argv):
    seed = int(argv[0]) if argv else 0
    out = run(seed)
    symbols = ["d" if r["degenerate"] else ("n" if r["tags"]["0.001"]["degenerate"] else "c")
               for r in out["rows"]]
    comp = emlib.compressibility(
        symbols, "one letter per final projection: d degenerate at 1e-6, n degenerate only "
                 "at 1e-3, c neither")
    final = [g for g in out["recount"] if "earlier_stage" not in (g["method"] or "")]
    changed = [g for g in final if g["reaching_h_degenerate"]]
    lines = "; ".join(f"{g['run'].split('_seed')[0]} {g['shape']} {g['method']}: reaching h "
                      f"{g['reaching_h']} -> {g['reaching_h_excluding_degenerate']} "
                      f"({g['n_degenerate']} of {g['n_projections']} degenerate)"
                      for g in final if g["n_degenerate"])
    ai = (f"OBSERVED: {out['n_rows']} final projections read from {out['sources_read']}; "
          f"absent: {out['sources_absent'] or 'none'}. Degenerate at 1e-6 (points coincide or "
          f"lie on a line): {out['n_degenerate']}; degenerate only at 1e-3: "
          f"{out['n_degenerate_only_at_1e-3']}. Groups with a degenerate final projection, "
          f"starts reaching h before -> after excluding them: {lines or 'none'}. Groups where "
          f"the count of starts reaching h changes: {len(changed)}. A point at the centre "
          f"(fixed by every mirror; not part of the degenerate tag): {out['n_central']} "
          f"projections, {out['n_central_not_degenerate']} of them not degenerate; starts "
          f"reaching h with a point at the centre: "
          f"{sum(g['reaching_h_central'] for g in out['recount'])}. "
          "AI reading: not interpretable atm beyond the rows.")
    p = emlib.write_run(NAME, seed, {"sources": list(SOURCES), "tols": list(TOLS),
                                     "tag_tol": float(TAG_TOL), "count_rel_tol": COUNT_TOL,
                                     "runs_dir": "runs"},
                        out, None, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
