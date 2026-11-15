"""m. coxeter_plane -- Petrie / Coxeter-plane projections, read for mirrors.

CC0. Co-authored. Stdlib only.

Four 4-polytopes are projected onto the Coxeter plane of their symmetry
group (primitives.coxeter_plane, from the simple roots and the Coxeter
number h), and the projection is read with
probes.mirror_symmetry_from_pairs (no angle grid), with the grid probe and
with rotational_symmetry at k_max 24.

  5-cell         A4, h = 5    (vertices in R^5, in the hyperplane sum = 0)
  16-cell        B4, h = 8    (cross-polytope)
  tesseract      B4, h = 8    (hypercube)
  24-cell        F4, h = 12

Prediction, PROPOSED before the run by the operator: the Coxeter-plane
projection has the maximal mirror count. "Maximal" is not something this
run can show: it logs the count beside h and beside the largest count the
360-frame tour found (runs/mirror_sweep_seed0.json), and nothing here
searches all orientations.

Beside open_questions.md Q1, Q2, Q3: those ask for rotation angles. This
run gives a projection basis, not a path of rotation angles.

  python3 experiments/coxeter_plane.py [seed]
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from toolkit import primitives as P  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "coxeter_plane"
REL_TOLS = (1e-6, 1e-9)
K_MAX = 24
N_AXES = 360
SWEEP = "runs/mirror_sweep_seed0.json"
PREDICTION = ("PROPOSED before the run, by the operator: the Petrie / Coxeter-plane "
              "projection has the maximal mirror count.")

A4 = [[1.0, -1.0, 0, 0, 0], [0, 1.0, -1.0, 0, 0], [0, 0, 1.0, -1.0, 0], [0, 0, 0, 1.0, -1.0]]
B4 = [[1.0, -1.0, 0, 0], [0, 1.0, -1.0, 0], [0, 0, 1.0, -1.0], [0, 0, 0, 1.0]]
F4 = [[0, 1.0, -1.0, 0], [0, 0, 1.0, -1.0], [0, 0, 0, 1.0], [0.5, -0.5, -0.5, -0.5]]


def simplex_r5():
    """5-cell as the centred unit vectors of R^5 (sum of coordinates 0)."""
    return [[(1.0 if i == k else 0.0) - 0.2 for i in range(5)] for k in range(5)]


CASES = (
    ("5-cell", "A4", 5, A4, simplex_r5, "simplex"),
    ("16-cell", "B4", 8, B4, lambda: P.cross_polytope(4), "crosspolytope"),
    ("tesseract", "B4", 8, B4, lambda: P.hypercube(4), None),
    ("24-cell", "F4", 12, F4, P.cell24, None),
)


def sweep_max(path=None):
    f = Path(path) if path else emlib.ROOT / SWEEP
    if not f.exists():
        return {}
    rec = json.loads(f.read_text())
    return {s["shape"]: {"max_axis_count_shape_alone": s["max_axis_count"]["shape"],
                         "frames": s["frames"]} for s in rec["outputs"]["sweep"]}


def case(name, group, h, roots, make, sweep_name, sweeps):
    verts = make()
    basis = P.coxeter_plane(roots, h)
    pts = P.project(verts, basis=basis)
    readings = []
    for tol in REL_TOLS:
        axes = Q.mirror_symmetry_from_pairs(pts, tol)
        grid = Q.mirror_symmetry(pts, N_AXES, tol)
        rot = Q.rotational_symmetry(pts, K_MAX, tol)
        readings.append({"rel_tol": tol,
                         "mirror_axes_from_pairs_deg": axes, "n_mirror_axes": len(axes),
                         "mirror_axes_grid_deg": grid, "n_mirror_axes_grid": len(grid),
                         "rotational_symmetry": rot,
                         "mirror_bound": {"axes": "every pair of points; no angle grid",
                                          "rel_tol": tol, "about": "centroid",
                                          "n_points": len(pts)},
                         "mirror_grid_bound": Q.mirror_symmetry_bound(pts, N_AXES, tol),
                         "rotational_symmetry_bound":
                             Q.rotational_symmetry_bound(pts, K_MAX, tol)})
    distinct = len({(round(x, 9), round(y, 9)) for x, y in pts})
    sw = sweeps.get(sweep_name) if sweep_name else None
    return {"name": name, "group": group, "coxeter_number_h": h,
            "simple_roots": roots, "plane_basis": basis,
            "n_vertices": len(verts), "n_distinct_projected_points": distinct,
            "projected": pts, "readings": readings,
            "n_mirror_axes_equals_h": readings[0]["n_mirror_axes"] == h,
            "sweep_reference": sw,
            "more_axes_than_sweep_max": (readings[0]["n_mirror_axes"]
                                         > sw["max_axis_count_shape_alone"]) if sw else None}


def run(seed=0, sweep_path=None):
    sweeps = sweep_max(sweep_path)
    return {"prediction": PREDICTION,
            "cases": [case(*c, sweeps) for c in CASES]}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    angles = {c["name"]: c["readings"][0]["mirror_axes_from_pairs_deg"] for c in out["cases"]}
    symbols = [chr(ord("a") + c["readings"][0]["n_mirror_axes"]) for c in out["cases"]]
    comp = emlib.compressibility(
        symbols, "mirror-axis count per case as a letter (a = 0), cases in order")
    rows = "; ".join(
        f"{c['name']} ({c['group']}, h = {c['coxeter_number_h']}, {c['n_vertices']} vertices, "
        f"{c['n_distinct_projected_points']} distinct projected points): "
        f"{c['readings'][0]['n_mirror_axes']} mirror axes with no grid, "
        f"{c['readings'][0]['n_mirror_axes_grid']} on the {N_AXES}-axis grid, rotational "
        f"symmetry {c['readings'][0]['rotational_symmetry']} at k_max {K_MAX}"
        + (f", tour maximum for the shape alone "
           f"{c['sweep_reference']['max_axis_count_shape_alone']} over "
           f"{c['sweep_reference']['frames']} frames" if c["sweep_reference"] else
           ", no tour reference")
        for c in out["cases"])
    same = all(c["readings"][0]["n_mirror_axes"] == c["readings"][1]["n_mirror_axes"]
               and c["readings"][0]["rotational_symmetry"] == c["readings"][1]["rotational_symmetry"]
               for c in out["cases"])
    ai = (f"OBSERVED, tolerance 1e-6: {rows}. Same counts at tolerance 1e-9: {same}. "
          f"Prediction logged before the reading: {PREDICTION} "
          "DERIVED: a plane figure with h-fold rotational symmetry and one mirror has "
          "exactly h mirror axes through its centre; whether any other orientation of "
          "the same polytope shows more is not tested by this run, which reads one "
          "plane per shape. AI reading: not interpretable atm beyond the counts.")
    p = emlib.write_run(NAME, seed,
                        {"cases": [[c[0], c[1], c[2]] for c in CASES],
                         "rel_tols": list(REL_TOLS), "k_max": K_MAX, "n_axes_grid": N_AXES,
                         "sweep_reference": SWEEP, "prediction": PREDICTION,
                         "about": "centroid"},
                        out, angles, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
