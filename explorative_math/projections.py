"""Finite sampled Givens tours, planar composition, and unfiltered comparisons.

This is not an exhaustive grand tour, projection pursuit, a Coxeter construction,
or a Penrose/quasicrystal generator. Every original vertex and every requested
frame is retained. Raw geometry is OBSERVED; angles and comparisons are DERIVED.
Catalog agreement means angular coincidence only, never a natural mechanism.

Public APIs accept a configuration mapping. ``max_vertices`` bounds the total
number of vertices generated for both shapes together, not just the first shape.
``iter_projection_runs`` yields exactly one run with one geometry definition and
``frames`` frame outputs. Supply a catalog list explicitly for filesystem-free
use, or package ``nature_catalog.json`` beside this module (zipapps supported).
"""
from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from importlib import resources
from itertools import product
import hashlib
import json
import math
import random
from typing import Any

PHI = (1.0 + math.sqrt(5.0)) / 2.0
GOLDEN_ANGLE_DEGREES = 360.0 / (PHI * PHI)
DISCLAIMER = "Angular coincidence only; no natural-mechanism inference."
SHAPES = ("crosspolytope", "hypercube", "simplex")
COMPOSITIONS = ("none", "phi", "second")


def _integer(value: Any, name: str, minimum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return value


def _finite(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    try:
        value = float(value)
    except OverflowError as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    return value


def _canonical(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"), allow_nan=False)


def validate_projection_config(config: Mapping[str, Any]) -> dict[str, Any]:
    """Normalize bounds and refuse oversized geometry before generating vertices.

    Defaults: dimension=5, crosspolytope, frames=12, seed=0, composition=phi,
    tolerance_degrees=1, max_vertices=4096. Integers cannot be booleans.
    Extra configuration keys are ignored, consistently with the geometry API.
    Composition ``second`` uses a second copy of the same original shape.
    """
    if not isinstance(config, Mapping):
        raise ValueError("config must be a mapping")
    dimension = _integer(config.get("dimension", 5), "dimension", 2)
    frames = _integer(config.get("frames", 12), "frames", 1)
    seed = _integer(config.get("seed", 0), "seed")
    limit = _integer(config.get("max_vertices", 4096), "max_vertices", 1)
    shape = config.get("shape", "crosspolytope")
    if shape not in SHAPES:
        raise ValueError("shape must be crosspolytope, hypercube, or simplex")
    composition = config.get("composition", "phi")
    if composition not in COMPOSITIONS:
        raise ValueError("composition must be none, phi, or second")
    tolerance = _finite(config.get("tolerance_degrees", 1.0), "tolerance_degrees")
    if not 0.0 <= tolerance <= 180.0:
        raise ValueError("tolerance_degrees must be in [0, 180]")
    copies = 1 if composition == "none" else 2
    # Check exponent against the bit length first; even absurd dimensions cannot
    # force allocation of a giant integer merely to reject a hypercube.
    per_shape_limit = limit // copies
    if shape == "hypercube" and dimension >= per_shape_limit.bit_length():
        raise ValueError(f"generated vertex count (2**{dimension} * {copies}) "
                         f"exceeds max_vertices={limit}; no vertices truncated")
    count = 2 ** dimension if shape == "hypercube" else (
        2 * dimension if shape == "crosspolytope" else dimension + 1)
    if copies * count > limit:
        raise ValueError(f"generated vertex count {copies * count} exceeds "
                         f"max_vertices={limit}; no vertices truncated")
    return {"dimension": dimension, "shape": shape, "frames": frames,
            "seed": seed, "composition": composition,
            "tolerance_degrees": tolerance, "max_vertices": limit}


def count_projection(config: Mapping[str, Any]) -> dict[str, int]:
    """Return exact archive counts without reading the catalog or sampling."""
    cfg = validate_projection_config(config)
    return {"runs": 1, "outputs": cfg["frames"] + 1}


def _validate_catalog(catalog: Any) -> list[dict[str, Any]]:
    # Also accept common packaged JSON envelopes while the explicit API remains
    # a list of records. Envelope metadata does not affect measurements.
    if isinstance(catalog, Mapping):
        keys = [key for key in ("records", "targets", "entries", "catalog")
                if key in catalog]
        if len(keys) != 1:
            raise ValueError("catalog object must contain one records/targets/entries/catalog list")
        catalog = catalog[keys[0]]
    if not isinstance(catalog, (list, tuple)):
        raise ValueError("catalog must be a list of records")
    checked: list[dict[str, Any]] = []
    for index, record in enumerate(catalog):
        if not isinstance(record, Mapping):
            raise ValueError(f"catalog record {index} must be an object")
        copy = dict(record)
        for key in ("id", "label", "source_url", "interpretation"):
            if not isinstance(copy.get(key), str) or not copy[key].strip():
                raise ValueError(f"catalog record {index} requires a nonempty {key}")
        target = _finite(copy.get("target_angle_degrees"), "target_angle_degrees")
        if not 0.0 <= target <= 180.0:
            raise ValueError("target_angle_degrees must be in [0, 180]")
        if copy.get("status") != "DERIVED":
            raise ValueError(f"catalog record {index} status must be DERIVED")
        copy["target_angle_degrees"] = target
        try:
            copy = json.loads(_canonical(copy))
        except (TypeError, ValueError) as exc:
            raise ValueError("catalog metadata must be finite JSON data") from exc
        checked.append(copy)
    return checked


def default_catalog() -> list[dict[str, Any]]:
    """Read the package-relative nature catalog, including from a zipapp.

    No source code is imported dynamically and no current-working-directory
    paths are consulted. Missing resources fail explicitly instead of silently
    omitting targets. Pass ``catalog=[...]`` to avoid resource loading entirely.
    """
    try:
        text = resources.files("explorative_math").joinpath("nature_catalog.json").read_text(encoding="utf-8")
        catalog = json.loads(text)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("cannot load package nature_catalog.json; package the resource "
                         "or pass catalog explicitly") from exc
    return _validate_catalog(catalog)


def _vertices(shape: str, n: int) -> list[list[float]]:
    """All vertices in stable order, centered with unit circumradius."""
    if shape == "crosspolytope":
        vertices = []
        for axis in range(n):
            for sign in (-1.0, 1.0):
                vector = [0.0] * n
                vector[axis] = sign
                vertices.append(vector)
        return vertices
    if shape == "hypercube":
        size = 1.0 / math.sqrt(n)
        return [[size * sign for sign in signs]
                for signs in product((-1.0, 1.0), repeat=n)]
    # Column k of the Helmert basis has k+1 leading +1/sqrt((k+1)(k+2)),
    # then -(k+1)/sqrt((k+1)(k+2)), then zeros. These n columns are an
    # orthonormal basis of the centered hyperplane in R**(n+1). Its rows
    # represent centered regular-simplex vertices in R**n; rescale to radius 1.
    radius_scale = math.sqrt((n + 1.0) / n)
    vertices = [[0.0] * n for _ in range(n + 1)]
    for k in range(n):
        size = radius_scale / math.sqrt((k + 1) * (k + 2))
        for row in range(k + 1):
            vertices[row][k] = size
        vertices[k + 1][k] = -(k + 1) * size
    return vertices


def _tour(n: int, seed: int) -> list[dict[str, Any]]:
    """Seeded ordering of a connected O(n) subset of coordinate planes.

    Neighbor planes plus the closing plane touch every coordinate. This subset
    and its phases describe a sampled path, not all coordinate planes or all
    projections. A linear phase-adjusted sweep makes frame zero the identity.
    """
    planes = [(i, i + 1) for i in range(n - 1)]
    if n > 2:
        planes.append((0, n - 1))
    rng = random.Random(seed)
    rng.shuffle(planes)
    result = []
    for order, plane in enumerate(planes):
        phase = rng.uniform(-180.0, 180.0)
        result.append({"order": order, "plane": list(plane),
                       "seeded_phase_degrees": phase,
                       "sweep_degrees": 360.0 + phase})
    return result


def _rotate(vertices: Sequence[Sequence[float]], rotations: list[dict[str, Any]]) -> list[list[float]]:
    result = [list(vertex) for vertex in vertices]
    for rotation in rotations:
        i, j = rotation["plane"]
        theta = math.radians(rotation["angle_degrees"])
        c, s = math.cos(theta), math.sin(theta)
        for vertex in result:
            x, y = vertex[i], vertex[j]
            vertex[i], vertex[j] = c * x - s * y, s * x + c * y
    return result


def _measure(points: Sequence[Sequence[float]], origins: list[dict[str, int]]) -> dict[str, Any]:
    """Keep one polar record per point and every cyclic neighboring valid angle.

    Exactly zero vectors have null angles and cannot enter the angular order.
    Nonzero floating-point vectors, however small, are not dropped. Equal
    angles are stably ordered by point index; their zero gaps are retained.
    """
    polar = []
    for index, (x, y) in enumerate(points):
        radius = math.hypot(x, y)
        angle = None if radius == 0.0 else math.degrees(math.atan2(y, x)) % 360.0
        polar.append({"point_index": index, **origins[index], "radius": radius,
                      "angle_degrees": angle, "status": "DERIVED",
                      "degenerate": angle is None,
                      "degeneracy_reason": "zero_projection" if angle is None else None})
    ordered = sorted((p for p in polar if p["angle_degrees"] is not None),
                     key=lambda p: (p["angle_degrees"], p["point_index"]))
    gaps = []
    for index, left in enumerate(ordered):
        right = ordered[(index + 1) % len(ordered)]
        # Add the complete wrap, including when every direction is identical.
        gap = right["angle_degrees"] - left["angle_degrees"]
        if index == len(ordered) - 1:
            gap += 360.0
        reason = "duplicate_direction" if gap == 0.0 else (
            "single_defined_direction" if len(ordered) == 1 else None)
        gaps.append({"measurement_index": index, "left_point_index": left["point_index"],
                     "right_point_index": right["point_index"],
                     "left_origin": dict(origins[left["point_index"]]),
                     "right_origin": dict(origins[right["point_index"]]),
                     "gap_degrees": gap, "wraparound": index == len(ordered) - 1,
                     "degenerate": reason is not None, "degeneracy_reason": reason,
                     "status": "DERIVED"})
    groups: dict[float, list[int]] = {}
    for point in ordered:
        groups.setdefault(point["angle_degrees"], []).append(point["point_index"])
    return {"status": "DERIVED", "angle_convention": "polar degrees in [0,360); cyclic CCW gaps in [0,360]",
            "polar": polar, "polar_angles_degrees": [p["angle_degrees"] for p in polar],
            "sorted_point_indices": [p["point_index"] for p in ordered],
            "angular_gaps": gaps,
            "degeneracies": {"undefined_point_indices": [p["point_index"] for p in polar
                                                         if p["angle_degrees"] is None],
                             "duplicate_direction_groups": [v for v in groups.values() if len(v) > 1],
                             "undefined_gap_reason": "no_defined_angles" if not ordered else None}}


def _compare(gaps: list[dict[str, Any]], catalog: list[dict[str, Any]],
             tolerance: float, frame_index: int) -> list[dict[str, Any]]:
    comparisons = []
    for catalog_index, target in enumerate(catalog):
        angle = target["target_angle_degrees"]
        errors = [min(abs(g["gap_degrees"] - angle),
                      360.0 - abs(g["gap_degrees"] - angle)) for g in gaps]
        nearest = min(range(len(errors)), key=lambda i: (errors[i], i)) if errors else None
        error = errors[nearest] if nearest is not None else None
        comparisons.append({"catalog_index": catalog_index, "catalog_id": target["id"],
                            "label": target["label"], "frame_index": frame_index,
                            "target_angle_degrees": angle,
                            "observed_nearest_degrees": gaps[nearest]["gap_degrees"] if nearest is not None else None,
                            "absolute_error_degrees": error, "tolerance_degrees": tolerance,
                            "matched": error is not None and error <= tolerance,
                            "result": "match" if error is not None and error <= tolerance else "miss",
                            "measurement_index": nearest,
                            "nearest_measurement_indices": [i for i, e in enumerate(errors) if e == error],
                            "matching_measurement_indices": [i for i, e in enumerate(errors) if e <= tolerance],
                            "measurement_errors_degrees": errors,
                            "nearest_measurement": dict(gaps[nearest]) if nearest is not None else None,
                            "reason": "no_defined_angular_gaps" if nearest is None else None,
                            "comparison_metric": "min(abs(gap-target), 360-abs(gap-target))",
                            "catalog_record": json.loads(_canonical(target)),
                            "source_url": target["source_url"], "interpretation": target["interpretation"],
                            "status": "DERIVED", "disclaimer": DISCLAIMER})
    return comparisons


def iter_projection_runs(config: Mapping[str, Any],
                         catalog: list[dict[str, Any]] | None = None) -> Iterator[dict[str, Any]]:
    """Yield one complete, deterministic run, with all frames and catalog misses.

    Each frame independently applies ordered Givens rotations to the original
    vertices. Frame k has evenly spaced progress k/frames (the closing endpoint
    is excluded) and plane angle progress*(360+seeded_phase) degrees. Frame zero
    is unrotated. Projection is onto coordinates 0 and 1.

    Phi composition retains the base and adds a copy scaled by phi, rotated
    in-plane by the golden angle. Second composition adds the same original
    shape at unit scale with that rotation. Uniform positive scaling alone
    preserves every defined angle; it does NOT create golden-angle signatures.
    Comparisons use all composed cyclic gaps, with exact measurement indices,
    tied nearest indices, every tolerance match, and all errors retained.
    """
    cfg = validate_projection_config(config)
    targets = default_catalog() if catalog is None else _validate_catalog(catalog)
    original = _vertices(cfg["shape"], cfg["dimension"])
    schedule = _tour(cfg["dimension"], cfg["seed"])
    composition = cfg["composition"]
    scale = PHI if composition == "phi" else 1.0
    shape_count = 1 if composition == "none" else 2
    params = dict(cfg)
    params.update({"catalog": targets, "normalization": "centered; unit circumradius per original shape",
                   "projection_coordinates": [0, 1],
                   "tour": "finite sampled Givens tour; not exhaustive; no projection pursuit",
                   "plane_selection": "adjacent coordinate planes plus closing plane; seeded ordering",
                   "angle_schedule": "angle_degrees = (frame_index / frames) * (360 + seeded_phase_degrees)",
                   "rotation_schedule": schedule, "shape_count": shape_count,
                   "generated_vertex_count": shape_count * len(original),
                   "max_vertices_scope": "total across all shapes; never truncate",
                   "composition_scale": scale, "composition_rotation_degrees": GOLDEN_ANGLE_DEGREES if shape_count == 2 else 0.0,
                   "scale_angle_invariance": "positive uniform scaling alone leaves polar angles and angular gaps unchanged",
                   "disclaimer": DISCLAIMER})
    identity = "projection-" + hashlib.sha256(_canonical({"config": cfg, "catalog": targets}).encode("utf-8")).hexdigest()
    shapes = [{"shape_index": index, "shape": cfg["shape"], "original_vertices": [v[:] for v in original],
               "vertex_count": len(original), "original_circumradius": 1.0,
               "composition_scale": 1.0 if index == 0 else scale,
               "composition_rotation_degrees": 0.0 if index == 0 else GOLDEN_ANGLE_DEGREES,
               "status": "OBSERVED"} for index in range(shape_count)]
    outputs: list[dict[str, Any]] = [{"kind": "geometry_definition", "status": "OBSERVED",
                                     "dimension": cfg["dimension"], "shape": cfg["shape"],
                                     "original_vertices": [v[:] for v in original], "shapes": shapes,
                                     "generated_vertex_count": shape_count * len(original),
                                     "normalization": params["normalization"]}]
    primary_origins = [{"shape_index": 0, "vertex_index": i} for i in range(len(original))]
    for frame_index in range(cfg["frames"]):
        progress = frame_index / cfg["frames"]
        rotations = [{**item, "angle_degrees": progress * item["sweep_degrees"]} for item in schedule]
        rotated = _rotate(original, rotations)
        projected = [v[:2] for v in rotated]
        before = _measure(projected, primary_origins)
        # Exact invariant angles are retained rather than allowing atan2 roundoff
        # on uniformly scaled copies to suggest a nonexistent angle change.
        pure_scale = json.loads(_canonical(before))
        for point in pure_scale["polar"]:
            point["radius"] *= scale
        pure_scale["derivation"] = "positive uniform scaling preserves directions; copied angles are exact invariants"
        scaled_points = [[scale * x for x in point] for point in projected]
        combined = [point[:] for point in projected]
        origins = [dict(origin) for origin in primary_origins]
        frame_shapes = [{"shape_index": 0, "rotated_vertices": [v[:] for v in rotated],
                         "projected_points": [p[:] for p in projected],
                         "rotation_metadata": rotations, "composition_scale": 1.0,
                         "composition_rotation_degrees": 0.0, "status": "OBSERVED"}]
        if shape_count == 2:
            copy_rotation = [{"order": 0, "plane": [0, 1], "angle_degrees": GOLDEN_ANGLE_DEGREES}]
            second_vertices = _rotate([[scale * x for x in v] for v in rotated], copy_rotation)
            second_points = [v[:2] for v in second_vertices]
            combined.extend(second_points)
            origins.extend({"shape_index": 1, "vertex_index": i} for i in range(len(original)))
            frame_shapes.append({"shape_index": 1, "rotated_vertices": second_vertices,
                                 "projected_points": second_points,
                                 "rotation_metadata": rotations + copy_rotation,
                                 "transform_order": "tour rotation, uniform scale, golden-angle in-plane rotation",
                                 "composition_scale": scale, "composition_rotation_degrees": GOLDEN_ANGLE_DEGREES,
                                 "status": "OBSERVED"})
        after = _measure(combined, origins)
        outputs.append({"kind": "projection_frame", "status": "OBSERVED", "frame_index": frame_index,
                        "rotation_metadata": {"method": "ordered Givens coordinate-plane rotations",
                                              "sampled_not_exhaustive": True, "seed": cfg["seed"],
                                              "progress": progress, "rotations": rotations},
                        "rotated_vertices": rotated, "projected_points": projected,
                        "projected_polar_angles_degrees": before["polar_angles_degrees"],
                        "precomposition_measurements": before,
                        "pure_scale_points": scaled_points, "pure_scale_measurements": pure_scale,
                        "shapes": frame_shapes, "composed_points": combined,
                        "composed_point_origins": origins,
                        "composed_measurements": after, "angular_gaps": after["angular_gaps"],
                        "catalog_comparisons": _compare(after["angular_gaps"], targets,
                                                        cfg["tolerance_degrees"], frame_index)})
    yield {"id": identity, "run_id": identity, "family": "projections", "params": params,
           "outputs": outputs,
           "human_interest": None, "ai_interest": None, "compressibility": None,
           "divergence_note": None, "status": "OBSERVED"}
