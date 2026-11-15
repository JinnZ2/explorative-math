"""l. mirror_sweep -- what frame 0 is, and mirrors along a 360-frame tour.

CC0. Co-authored. Stdlib only.

runs/symmetry_resolution_seed0.json read a mirror only on frame 0 of each
shape (4 frames). Two things are logged here.

  frame_zero   what frame 0 of Build A's tour is: its rotation angles, and
               whether its projected points equal the first two coordinates
               of the unrotated vertices.
  sweep        the same tour at 360 frames per shape (frame k is at k/360 of
               each plane's sweep; frames 0, 90, 180, 270 are the four frames
               read before). On every frame, mirror axes with no angle grid
               (probes.mirror_symmetry_from_pairs) of three point sets: the
               projected shape alone, Build A's union (phi scale + rotate)
               and Build B's union (phi scale).

  residual     on the same frames, probes.mirror_residual of each set: a
               continuous reading (0 = exact mirror), so a frame that comes
               near a mirror without having one shows as a dip.

Prediction, PROPOSED before the run by the operator: frames with a mirror
are rare and clustered at special orientations.

Build A (explorative_math) is imported read-only, as compare/ does.

  python3 experiments/mirror_sweep.py [seed]
"""
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from explorative_math import projections as A_proj  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "mirror_sweep"
SHAPES = ("simplex", "crosspolytope")
DIMENSION = 4
FRAMES = 360
REL_TOL = 1e-6
SETS = ("shape", "a_union", "b_union")
PREDICTION = ("PROPOSED before the run, by the operator: frames with any mirror are rare, "
              "and clustered at special orientations.")


def frames_of(shape, frames, seed=0):
    run = next(A_proj.iter_projection_runs(
        {"dimension": DIMENSION, "shape": shape, "frames": frames, "seed": seed,
         "composition": "phi"}))
    original = next(o for o in run["outputs"] if o["kind"] == "geometry_definition")
    return original["original_vertices"], [o for o in run["outputs"]
                                           if o["kind"] == "projection_frame"]


def frame_zero(shape):
    original, frames = frames_of(shape, 4)
    f = frames[0]
    rots = f["rotation_metadata"]["rotations"]
    pts = f["projected_points"]
    first_two = [v[:2] for v in original]
    return {"shape": shape, "frame_index": f["frame_index"],
            "rotation_angles_degrees": [[r["plane"], r["angle_degrees"]] for r in rots],
            "all_rotation_angles_zero": all(r["angle_degrees"] == 0.0 for r in rots),
            "original_vertices": original, "projected_points": pts,
            "projected_equals_first_two_coordinates": pts == first_two,
            "max_abs_difference": max(abs(a - b) for p, q in zip(pts, first_two)
                                      for a, b in zip(p, q)),
            "degeneracy": Q.degeneracy(pts)}


def sweep(shape, frames=FRAMES, rel_tol=REL_TOL):
    _, fs = frames_of(shape, frames)
    counts = {k: [] for k in SETS}
    resid = {k: [] for k in SETS}
    with_mirror = []
    for f in fs:
        pts = f["projected_points"]
        sets = {"shape": pts, "a_union": f["composed_points"],
                "b_union": [r["point"] for r in C.phi_stack(pts, (0, 1))]}
        axes = {k: Q.mirror_symmetry_from_pairs(v, rel_tol) for k, v in sets.items()}
        for k in SETS:
            counts[k].append(len(axes[k]))
            resid[k].append(Q.mirror_residual(sets[k])["residual"])
        if any(axes.values()):
            rots = _rotations(f)
            with_mirror.append({"frame_index": f["frame_index"],
                                "degeneracy": {k: Q.degeneracy(v) for k, v in sets.items()},
                                "degenerate": {k: Q.degeneracy(v)["tags"]["1e-06"]["degenerate"]
                                               for k, v in sets.items()},
                                "progress": f["frame_index"] / frames,
                                "rotation_angles_degrees": rots,
                                "axes_deg": axes})
    n = len(fs)
    result = _sweep_result(shape, fs, n, rel_tol, counts, with_mirror, resid)
    deg = {}
    for k in SETS:
        top = max(counts[k])
        hits = [m for m in with_mirror if len(m["axes_deg"][k]) == top]
        deg[k] = (all(m["degenerate"][k] for m in hits) if top and hits else None)
    result["max_axis_count_degenerate"] = deg
    result["max_axis_count_excluding_degenerate"] = {
        k: max([len(m["axes_deg"][k]) for m in with_mirror if not m["degenerate"][k]] + [0])
        for k in SETS}
    return result


def _sweep_result(shape, fs, n, rel_tol, counts, with_mirror, resid):
    return {"shape": shape, "frames": n,
            "bound": {"frames": n, "frame_step": f"1/{n} of each plane's sweep",
                      "axes": "every pair of points; no angle grid", "rel_tol": rel_tol,
                      "about": "centroid",
                      "n_points": {"shape": len(fs[0]["projected_points"]),
                                   "a_union": len(fs[0]["composed_points"]),
                                   "b_union": 2 * len(fs[0]["projected_points"])}},
            "axis_count_per_frame": counts,
            "frames_with_any_mirror": {k: sum(1 for c in counts[k] if c) for k in SETS},
            "fraction_with_any_mirror": {k: sum(1 for c in counts[k] if c) / n for k in SETS},
            "max_axis_count": {k: max(counts[k]) for k in SETS},
            "max_axis_count_degenerate": None,      # filled below: is that maximum read on
                                                    # a projection with coincident points?
            "mirror_frames": with_mirror,
            # continuous reading, added beside the hit / miss counts
            "residual_per_frame": resid,
            "residual_bound": {"axes": "mean polar angle of every pair of points",
                               "matching": "nearest point, not one-to-one",
                               "normalised_by": "largest distance from the centroid",
                               "frames": n},
            "residual_min_excluding_frame_0": {
                k: min((r, i) for i, r in enumerate(resid[k]) if i != 0)[::-1]
                for k in SETS} if n > 1 else None,
            "residual_median": {k: statistics.median(resid[k]) for k in SETS},
            "residual_local_minima": {k: local_minima(resid[k]) for k in SETS}}


def local_minima(values, keep=8):
    """Frames whose residual is below both neighbours (the tour is treated as
    open: frame 0 and the last frame count if below their one neighbour).
    The `keep` lowest, lowest first, as [frame_index, residual]."""
    n = len(values)
    out = []
    for i, v in enumerate(values):
        left = values[i - 1] if i > 0 else float("inf")
        right = values[i + 1] if i < n - 1 else float("inf")
        if v < left and v < right:
            out.append([i, v])
    out.sort(key=lambda r: r[1])
    return out[:keep]


def _rotations(frame):
    return [[r["plane"], r["angle_degrees"]] for r in frame["rotation_metadata"]["rotations"]]


def _fmt(pair):
    return f"[{pair[0]}, {pair[1]:.3e}]"


def run(seed=0, frames=FRAMES):
    return {"prediction": PREDICTION,
            "frame_zero": [frame_zero(s) for s in SHAPES],
            "sweep": [sweep(s, frames) for s in SHAPES]}


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    symbols = [str(min(c, 9)) for s in out["sweep"] for c in s["axis_count_per_frame"]["b_union"]]
    comp = emlib.compressibility(
        symbols, "mirror-axis count of Build B's union per frame (digit, capped at 9), "
                 "simplex frames then cross-polytope frames")
    z = "; ".join(f"{f['shape']}: all rotation angles zero {f['all_rotation_angles_zero']}, "
                  f"projected points equal the first two coordinates "
                  f"{f['projected_equals_first_two_coordinates']}" for f in out["frame_zero"])
    sw = "; ".join(
        f"{s['shape']}: frames with a mirror of {s['frames']} -- shape alone "
        f"{s['frames_with_any_mirror']['shape']}, A union {s['frames_with_any_mirror']['a_union']}, "
        f"B union {s['frames_with_any_mirror']['b_union']}; frame indices "
        f"{[m['frame_index'] for m in s['mirror_frames']]}; largest axis count "
        f"{s['max_axis_count']}" for s in out["sweep"])
    rs = "; ".join(
        f"{s['shape']}: lowest residual off frame 0 [frame, value] shape alone "
        f"{_fmt(s['residual_min_excluding_frame_0']['shape'])}, A union "
        f"{_fmt(s['residual_min_excluding_frame_0']['a_union'])}, B union "
        f"{_fmt(s['residual_min_excluding_frame_0']['b_union'])}; median over frames, shape alone "
        f"{s['residual_median']['shape']:.3e}; lowest local minima of "
        f"the shape alone {[_fmt(m) for m in s['residual_local_minima']['shape'][:4]]}"
        for s in out["sweep"])
    ai = (f"OBSERVED, frame 0: {z}. OBSERVED, sweep at {FRAMES} frames, tolerance {REL_TOL:g}, "
          f"no axis grid: {sw}. Prediction logged before the reading: {PREDICTION} "
          f"OBSERVED, residual (0 = exact mirror) -- {rs}. "
          "DERIVED: the sweep reads one tour (seed 0) at steps of 1/360 of each plane's "
          "sweep; an orientation with a mirror that falls between two frames is not seen, "
          "so the counts are bounded by the frame step and by the path. "
          "AI reading: not interpretable atm beyond the counts.")
    p = emlib.write_run(NAME, seed,
                        {"shapes": list(SHAPES), "dimension": DIMENSION, "frames": FRAMES,
                         "rel_tol": REL_TOL, "tour_seed": 0, "prediction": PREDICTION,
                         "source": "explorative_math.projections.iter_projection_runs, "
                                   "composition phi"},
                        out, None, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
