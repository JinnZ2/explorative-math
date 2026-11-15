"""ab. collapse_gap -- collapse a shape to one point and read it from the complement.

CC0. Co-authored. Stdlib only.

Carried method; status PROPOSED. This file is one concrete reading of it;
the choices that make it concrete are in "What is done here" and in the
run's bound.

  collapse   the shape's points are identified with one point, the collapse
             point (the centroid). What remains is the complement, sampled.
  baseline   the ORIGINAL shape's pattern is placed around the collapse point
             and read from the samples only: distance field, the turns and
             mirrors that leave the field unchanged, voids, angular gaps.
  probes     the same operation with other shapes; each reading is set
             beside the baseline's, and the difference logged per probe.
  checks     does the complement reading recover the shape's own symmetry?
             control: a random point set, collapsed the same way.

What is done here
  plane      shapes as plane patterns: 4-polytopes on their Coxeter planes,
             and regular polygons. Samples: a polar grid inside the unit
             disc. The pattern is scaled to radius 1 and turned so a
             farthest point is at angle 0.
  sphere     shapes in full R^4: samples are directions; the field is the
             angle to the nearest vertex direction; voids are found by a
             seeded hill climb, and set beside the vertex directions of each
             probe shape, in the frame the shapes are built in.
  exact      the pentagon and decagon patterns in exact a + b*phi forms: the
             field is a squared distance, compared with no tolerance.

Pointer, structural only: Alexander duality relates a set to its complement.
Nothing here computes homology; no claim is made beyond the readings.

Observer slots null.

  python3 experiments/collapse_gap.py [seed]
  python3 experiments/collapse_gap.py --small [seed]   # phone bound; own run file
"""
import math
import random
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from exact import Cyc5, ZETA  # noqa: E402
from exact import complement as XC  # noqa: E402
from exact import inflation as XI  # noqa: E402
from toolkit import complement as CP  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "collapse_gap"
N_RINGS, N_THETA = 6, 1680
K_MAX = 48
SPHERE_STARTS, SPHERE_STEPS = 200, 400
BASELINES = ("5-cell", "16-cell", "24-cell")
SMALL = {"n_rings": 4, "n_theta": 840, "sphere_starts": 40, "sphere_steps": 250,
         "baselines": ("5-cell", "16-cell")}
METHOD = ("PROPOSED, carried: collapse the shape to one point; the remaining space is the "
          "complement; read the shape's pattern from the complement, then read other shapes "
          "the same way and log each against the baseline.")

A4 = [[1.0, -1.0, 0, 0, 0], [0, 1.0, -1.0, 0, 0], [0, 0, 1.0, -1.0, 0], [0, 0, 0, 1.0, -1.0]]
B4 = [[1.0, -1.0, 0, 0], [0, 1.0, -1.0, 0], [0, 0, 1.0, -1.0], [0, 0, 0, 1.0]]
F4 = [[0, 1.0, -1.0, 0], [0, 0, 1.0, -1.0], [0, 0, 0, 1.0], [0.5, -0.5, -0.5, -0.5]]


def simplex_r5():
    return [[(1.0 if i == k else 0.0) - 0.2 for i in range(5)] for k in range(5)]


def plane_shapes(seed):
    rng = random.Random(f"{seed}/plane-controls")
    out = [("5-cell", "Coxeter plane A4", P.project(simplex_r5(), basis=P.coxeter_plane(A4, 5))),
           ("16-cell", "Coxeter plane B4",
            P.project(P.cross_polytope(4), basis=P.coxeter_plane(B4, 8))),
           ("tesseract", "Coxeter plane B4",
            P.project(P.hypercube(4), basis=P.coxeter_plane(B4, 8))),
           ("24-cell", "Coxeter plane F4", P.project(P.cell24(), basis=P.coxeter_plane(F4, 12)))]
    out += [(f"{k}-gon", "regular polygon", P.regular_polygon(k)) for k in (3, 5, 6, 7, 8, 12)]
    out += [(f"random_{n}", "control: random points in the unit square",
             [[rng.uniform(-1, 1), rng.uniform(-1, 1)] for _ in range(n)]) for n in (5, 8)]
    return out


def own_symmetry(pattern):
    """The pattern's symmetry read on its POINTS (not through the complement)."""
    return {"rotational_orders": Q.rotational_symmetry(pattern, K_MAX, 1e-6),
            "n_mirror_axes": len(Q.mirror_symmetry_from_pairs(pattern, 1e-6)),
            "bound": Q.rotational_symmetry_bound(pattern, K_MAX, 1e-6)}


def plane_reading(name, how, pts, n_rings, n_theta, keep_field):
    r = CP.read_plane(pts, n_rings, n_theta, K_MAX)
    own = own_symmetry(r["pattern"])
    tested = set(r["bound"]["orders_tested"])
    own_on_grid = [k for k in own["rotational_orders"] if k in tested]
    field = r.pop("field")
    r.update({"shape": name, "construction": how, "own_symmetry_from_points": own,
              "recovers_rotational_orders": r["field_rotational_orders"] == own_on_grid,
              "own_orders_the_grid_cannot_test":
                  [k for k in own["rotational_orders"] if k not in tested],
              "recovers_mirror_count":
                  len(r["field_mirror_axes_deg"]) == own["n_mirror_axes"],
              "n_field_mirror_axes": len(r["field_mirror_axes_deg"]),
              "any_symmetry_read_from_the_complement":
                  bool(r["field_rotational_orders"] or r["field_mirror_axes_deg"]),
              "degeneracy": Q.degeneracy(r["pattern"]),
              "middle_ring": field[len(field) // 2][::max(1, n_theta // 120)]})
    if keep_field:
        r["field"] = field
    return r, field


def delta(base, bfield, probe, pfield, n_theta):
    d = CP.field_difference(bfield, pfield, n_theta)
    bg = [g["value"] for g in base["distinct_gaps"]]
    pg = [g["value"] for g in probe["distinct_gaps"]]
    return {"probe": probe["shape"],
            "field_rms_min_over_turns": d["rms_min_over_turns"], "at_turn_deg": d["turn_deg"],
            "field_rms_at_zero_turn": d["rms_at_zero_turn"], "field_bound": d["bound"],
            "same_field_rotational_orders":
                base["field_rotational_orders"] == probe["field_rotational_orders"],
            "orders_only_in_baseline": [k for k in base["field_rotational_orders"]
                                        if k not in probe["field_rotational_orders"]],
            "orders_only_in_probe": [k for k in probe["field_rotational_orders"]
                                     if k not in base["field_rotational_orders"]],
            "mirror_axes_baseline_minus_probe":
                base["n_field_mirror_axes"] - probe["n_field_mirror_axes"],
            "n_points_baseline_minus_probe": base["n_points"] - probe["n_points"],
            "n_gaps_baseline_minus_probe": len(base["gaps_deg"]) - len(probe["gaps_deg"]),
            "same_distinct_gaps": (len(bg) == len(pg)
                                   and all(abs(a - b) < 1e-6 for a, b in zip(bg, pg))),
            "n_voids_inside_baseline_minus_probe":
                base["n_voids_inside"] - probe["n_voids_inside"],
            "centre_depth_baseline_minus_probe": base["centre_depth"] - probe["centre_depth"],
            "field_mean_baseline_minus_probe": base["field_mean"] - probe["field_mean"],
            "radial_profile_baseline_minus_probe":
                [a - b for a, b in zip(base["radial_profile"], probe["radial_profile"])]}


def sphere_shapes(seed):
    rng = random.Random(f"{seed}/sphere-controls")
    sim = P.simplex(4)
    cross, cube, c24 = P.cross_polytope(4), P.hypercube(4), P.cell24()
    return [("5-cell", sim), ("5-cell negated", [[-x for x in v] for v in sim]),
            ("16-cell", cross), ("tesseract", cube), ("24-cell", c24),
            ("16-cell + tesseract directions", cross + cube),
            ("random_8", [[rng.gauss(0, 1) for _ in range(4)] for _ in range(8)])]


def sphere_part(seed, starts, steps):
    shapes = sphere_shapes(seed)
    rows = []
    for name, verts in shapes:
        v = CP.sphere_voids(verts, seed, starts, steps)
        dirs = [x["direction"] for x in v["voids"]]
        matches = []
        for pname, pverts in shapes:
            m = CP.match_directions(dirs, pverts)
            m.update({"probe": pname,
                      "voids_sit_on_probe_directions":
                          m["max_angle_deg"] is not None and m["max_angle_deg"] <= 2.0,
                      "every_probe_direction_has_a_void":
                          m["candidates_hit"] == m["n_candidates"]})
            matches.append(m)
        rows.append({"shape": name, "n_vertices": len(verts),
                     "deepest_void_depth_deg": v["deepest_depth_deg"],
                     "n_deepest_voids_found": v["n_voids_found"],
                     "end_depths_deg": v["all_end_depths_deg"],
                     "voids": v["voids"], "bound": v["bound"],
                     "matches": matches,
                     "probes_the_voids_sit_on":
                         [m["probe"] for m in matches if m["voids_sit_on_probe_directions"]],
                     "probes_fully_covered":
                         [m["probe"] for m in matches if m["voids_sit_on_probe_directions"]
                          and m["every_probe_direction_has_a_void"]]})
    return {"frame": "each shape as toolkit.primitives builds it; no rotation is searched",
            "match_threshold_deg": 2.0, "shapes": rows}


def exact_part():
    pent = XI.pentagon()
    dec = pent + [z * XC.Z10 for z in pent]
    lopsided = [Cyc5(1, 0), ZETA, ZETA * ZETA * 2 + Cyc5(Fraction(1, 3), 0)]
    out = {}
    for name, pat, own in (("pentagon", pent, "turns of 72 and 144, mirror; not 36"),
                           ("decagon", dec, "turns of 36, 72 and 144, mirror"),
                           ("control: three points with no symmetry", lopsided, "none")):
        r = XC.read(pat)
        r["own_symmetry"] = own
        out[name] = r
    return out


def run(seed=0, n_rings=N_RINGS, n_theta=N_THETA, sphere_starts=SPHERE_STARTS,
        sphere_steps=SPHERE_STEPS, baselines=BASELINES):
    readings, fields = {}, {}
    for name, how, pts in plane_shapes(seed):
        readings[name], fields[name] = plane_reading(name, how, pts, n_rings, n_theta,
                                                     keep_field=name in baselines)
    per_baseline = []
    for b in baselines:
        deltas = [delta(readings[b], fields[b], readings[p], fields[p], n_theta)
                  for p in readings if p != b]
        deltas.sort(key=lambda d: d["field_rms_min_over_turns"])
        per_baseline.append({"baseline": b, "deltas": deltas,
                             "probes_by_field_difference": [d["probe"] for d in deltas],
                             "probes_with_the_same_field_orders":
                                 [d["probe"] for d in deltas
                                  if d["same_field_rotational_orders"]]})
    checks = [{"shape": n, "own_rotational_orders": r["own_symmetry_from_points"]["rotational_orders"],
               "field_rotational_orders": r["field_rotational_orders"],
               "own_mirror_axes": r["own_symmetry_from_points"]["n_mirror_axes"],
               "field_mirror_axes": r["n_field_mirror_axes"],
               "recovers_rotational_orders": r["recovers_rotational_orders"],
               "recovers_mirror_count": r["recovers_mirror_count"],
               "any_symmetry_read_from_the_complement":
                   r["any_symmetry_read_from_the_complement"],
               "control": n.startswith("random")} for n, r in readings.items()]
    return {"method": METHOD,
            "bound": {"seed": seed, "n_rings": n_rings, "n_theta": n_theta, "k_max": K_MAX,
                      "sphere_starts": sphere_starts, "sphere_steps": sphere_steps,
                      "baselines": list(baselines),
                      "collapse_point": "the centroid of the pattern",
                      "complement_sampled": "plane: polar grid inside the unit disc; "
                                            "sphere: directions reached by a seeded hill climb",
                      "full_field_logged_for": list(baselines),
                      "not_done": "no homology is computed; the duality pointer is structural"},
            "plane": {"readings": readings, "per_baseline": per_baseline, "checks": checks},
            "sphere": sphere_part(seed, sphere_starts, sphere_steps),
            "exact": exact_part()}


def main(argv):
    small = "--small" in argv
    argv = [a for a in argv if a != "--small"]
    seed = int(argv[0]) if argv else 0
    out = (run(seed, SMALL["n_rings"], SMALL["n_theta"], SMALL["sphere_starts"],
               SMALL["sphere_steps"], SMALL["baselines"]) if small else run(seed))
    symbols = ["y" if c["recovers_rotational_orders"] and c["recovers_mirror_count"] else "n"
               for c in out["plane"]["checks"]]
    comp = emlib.compressibility(
        symbols, "one letter per plane shape: y the complement reading gives the shape's own "
                 "rotational orders and mirror count, n it does not")
    angles = {n: r["field_mirror_axes_deg"] for n, r in out["plane"]["readings"].items()}
    b = out["bound"]
    # observer slots null: ai_interest is void
    p = emlib.write_run(NAME + "_small" if small else NAME, seed,
                        {"small": small, "n_rings": b["n_rings"], "n_theta": b["n_theta"],
                         "k_max": K_MAX, "sphere_starts": b["sphere_starts"],
                         "sphere_steps": b["sphere_steps"], "baselines": b["baselines"],
                         "method": METHOD, "status_before_run": "PROPOSED"},
                        out, angles, None, comp, None, gz=True)
    print(p)
    print("PLANE checks  (own orders / field orders | own mirrors / field mirrors)")
    for c in out["plane"]["checks"]:
        print(f"  {c['shape']:10s} {c['own_rotational_orders']} / {c['field_rotational_orders']} | "
              f"{c['own_mirror_axes']} / {c['field_mirror_axes']}  recovers "
              f"{c['recovers_rotational_orders']}, {c['recovers_mirror_count']}")
    for pb in out["plane"]["per_baseline"]:
        print("PLANE baseline", pb["baseline"])
        for d in pb["deltas"]:
            print(f"  {d['probe']:10s} field rms {d['field_rms_min_over_turns']:.4f} at "
                  f"{d['at_turn_deg']:.1f} deg  same orders {d['same_field_rotational_orders']}  "
                  f"gaps same {d['same_distinct_gaps']}  centre depth diff "
                  f"{d['centre_depth_baseline_minus_probe']:+.3f}  voids diff "
                  f"{d['n_voids_inside_baseline_minus_probe']:+d}")
    print("SPHERE")
    for r in out["sphere"]["shapes"]:
        print(f"  {r['shape']:30s} deepest void {r['deepest_void_depth_deg']:.2f} deg, found "
              f"{r['n_deepest_voids_found']}; sit on {r['probes_the_voids_sit_on']}; fully cover "
              f"{r['probes_fully_covered']}")
    print("EXACT")
    for n, r in out["exact"].items():
        print(f"  {n}: turns {r['field_unchanged_by_turn_deg']} mirror "
              f"{r['field_unchanged_by_mirror_in_real_axis']}")


if __name__ == "__main__":
    main(sys.argv[1:])
