"""Finite, seeded high-dimensional geometry experiments; no selection by interest.

Raw samples and complete paths are OBSERVED execution results. Mathematical
comparisons and aggregate measurements are explicitly DERIVED outputs.
"""
from __future__ import annotations

import hashlib
import json
import math
import random
from collections.abc import Iterator, Mapping
from typing import Any

EXPERIMENTS = ("orthogonality", "shell", "superposition", "trajectories")
DEFAULT_DIMENSIONS = [2, 8, 32, 128]


def _integer(value: Any, name: str, minimum: int | None = None) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return value


def _finite(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    return value


def validate_geometry_config(config: Mapping[str, Any]) -> dict[str, Any]:
    """Normalize finite bounds, rejecting ambiguous or invalid input."""
    experiment = config.get("experiment", "all")
    if experiment not in (*EXPERIMENTS, "all"):
        raise ValueError("unknown geometry experiment")
    dimensions = config.get("dimensions", DEFAULT_DIMENSIONS)
    if isinstance(dimensions, str):
        try:
            dimensions = [int(item.strip()) for item in dimensions.split(",")]
        except ValueError as exc:
            raise ValueError("dimensions must be comma-separated positive integers") from exc
    if not isinstance(dimensions, (list, tuple)) or not dimensions:
        raise ValueError("dimensions must be a nonempty list")
    dimensions = [_integer(n, "dimension", 1) for n in dimensions]
    if len(set(dimensions)) != len(dimensions):
        raise ValueError("dimensions must not contain duplicates")
    samples = _integer(config.get("samples", 32), "samples", 1)
    vectors = _integer(config.get("vectors", 0), "vectors", 0)
    steps = _integer(config.get("steps", 32), "steps", 0)
    seed = _integer(config.get("seed", 0), "seed")
    thickness = _finite(config.get("shell_thickness", 0.1), "shell_thickness")
    scale = _finite(config.get("map_scale", 0.8), "map_scale")
    if not 0 <= thickness <= 1:
        raise ValueError("shell_thickness must be in [0, 1]")
    if abs(scale) > 1:
        raise ValueError("map_scale must be in [-1, 1]")
    if experiment in ("all", "superposition") and vectors and vectors <= max(dimensions):
        raise ValueError("superposition requires vectors > every selected dimension; use 0 for 2*n")
    return {"experiment": experiment, "dimensions": dimensions, "samples": samples,
            "vectors": vectors, "steps": steps, "seed": seed,
            "shell_thickness": thickness, "map_scale": scale}


def count_geometry(config: Mapping[str, Any]) -> dict[str, int]:
    """Return exact run and output-entry counts without sampling anything."""
    cfg = validate_geometry_config(config)
    runs = outputs = 0
    selected = EXPERIMENTS if cfg["experiment"] == "all" else (cfg["experiment"],)
    for n in cfg["dimensions"]:
        for experiment in selected:
            if experiment in ("orthogonality", "shell"):
                runs += 1
                outputs += cfg["samples"] + 1
            elif experiment == "superposition":
                k = cfg["vectors"] or 2 * n
                runs += 1
                outputs += k + k * (k - 1) // 2 + 1
            else:
                runs += 2
                outputs += 2 * (cfg["steps"] + 1)
    return {"runs": runs, "outputs": outputs}


def unit_vector(rng: random.Random, dimension: int) -> list[float]:
    """Sample a sphere direction by normalizing independent standard Gaussians."""
    _integer(dimension, "dimension", 1)
    while True:
        vector = [rng.gauss(0.0, 1.0) for _ in range(dimension)]
        norm = math.sqrt(math.fsum(x * x for x in vector))
        if norm > 0:
            return [x / norm for x in vector]


def angle_measurement(left: list[float], right: list[float]) -> dict[str, float]:
    """Measure the oriented angle of unit vectors, clamping rounding at +/-1."""
    if len(left) != len(right) or not left:
        raise ValueError("angle vectors must have equal nonzero dimension")
    raw_dot = math.fsum(a * b for a, b in zip(left, right))
    dot = min(1.0, max(-1.0, raw_dot))
    return {"raw_dot": raw_dot, "dot": dot, "angle_degrees": math.degrees(math.acos(dot))}


def shell_fraction(dimension: int, thickness: float) -> float:
    """DERIVED: unit-ball volume in radius [1-t,1] is 1-(1-t)**n."""
    _integer(dimension, "dimension", 1)
    thickness = _finite(thickness, "thickness")
    if not 0 <= thickness <= 1:
        raise ValueError("thickness must be in [0, 1]")
    if thickness == 1:
        return 1.0
    return -math.expm1(dimension * math.log1p(-thickness))


def _run(params: dict[str, Any], outputs: list[dict[str, Any]]) -> dict[str, Any]:
    canonical = json.dumps(params, sort_keys=True, separators=(",", ":"), allow_nan=False)
    identity = "geom-" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:20]
    return {"id": identity, "run_id": identity, "family": "geometry", "params": params,
            "outputs": outputs, "human_interest": None, "ai_interest": None,
            "compressibility": None, "divergence_note": None, "status": "OBSERVED"}


def _rng(params: dict[str, Any]) -> random.Random:
    """Use an independent deterministic stream for each parameterized run."""
    text = json.dumps(params, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return random.Random(int.from_bytes(hashlib.sha256(text.encode("utf-8")).digest(), "big"))


def _orthogonality(n: int, cfg: dict[str, Any]) -> dict[str, Any]:
    params = {"experiment": "orthogonality", "dimension": n, "samples": cfg["samples"],
              "seed": cfg["seed"], "sampling": "independent normalized standard Gaussian pairs",
              "angle_convention": "oriented degrees in [0,180]"}
    rng = _rng(params)
    outputs: list[dict[str, Any]] = []
    for index in range(cfg["samples"]):
        left, right = unit_vector(rng, n), unit_vector(rng, n)
        outputs.append({"kind": "pair", "index": index, "left": left, "right": right,
                        **angle_measurement(left, right), "status": "OBSERVED"})
    outputs.append({"kind": "summary", "status": "DERIVED", "count": cfg["samples"],
                    "mean_angle_degrees": math.fsum(x["angle_degrees"] for x in outputs) / cfg["samples"],
                    "mean_squared_dot": math.fsum(x["dot"] ** 2 for x in outputs) / cfg["samples"],
                    "expected_squared_dot": 1.0 / n,
                    "expected_angle_degrees": 90.0,
                    "derivation": "sphere symmetry gives E[dot^2]=1/n and E[angle]=90 degrees; empirical means use all pairs"})
    return _run(params, outputs)


def _shell(n: int, cfg: dict[str, Any]) -> dict[str, Any]:
    thickness = cfg["shell_thickness"]
    params = {"experiment": "shell", "dimension": n, "samples": cfg["samples"],
              "seed": cfg["seed"], "shell_thickness": thickness, "outer_radius": 1.0,
              "sampling": "uniform unit-ball volume: direction * U**(1/n)"}
    rng = _rng(params)
    outputs: list[dict[str, Any]] = []
    for index in range(cfg["samples"]):
        direction = unit_vector(rng, n)
        uniform = rng.random()
        radius = uniform ** (1.0 / n)
        outputs.append({"kind": "point", "index": index, "uniform": uniform,
                        "radius": radius, "point": [radius * x for x in direction],
                        "in_shell": radius >= 1.0 - thickness, "status": "OBSERVED"})
    outputs.append({"kind": "volume_comparison", "status": "DERIVED",
                    "exact_fraction": shell_fraction(n, thickness),
                    "empirical_fraction": sum(x["in_shell"] for x in outputs) / cfg["samples"],
                    "derivation": "ball volume scales as radius**n; shell fraction = 1-(1-thickness)**dimension"})
    return _run(params, outputs)


def _superposition(n: int, cfg: dict[str, Any]) -> dict[str, Any]:
    k = cfg["vectors"] or 2 * n
    params = {"experiment": "superposition", "dimension": n, "vectors": k,
              "seed": cfg["seed"], "sampling": "independent normalized standard Gaussian directions",
              "packing_method": "random baseline, not an optimized or proven-optimal packing",
              "angle_convention": "oriented directions, antipodes distinct"}
    rng = _rng(params)
    directions = [unit_vector(rng, n) for _ in range(k)]
    outputs: list[dict[str, Any]] = [
        {"kind": "direction", "index": index, "vector": vector, "status": "OBSERVED"}
        for index, vector in enumerate(directions)
    ]
    minimum = math.inf
    closest: list[int] | None = None
    max_abs_dot = 0.0
    for i in range(k):
        for j in range(i + 1, k):
            measurement = angle_measurement(directions[i], directions[j])
            outputs.append({"kind": "pair_angle", "left_index": i, "right_index": j,
                            **measurement, "status": "OBSERVED"})
            if measurement["angle_degrees"] < minimum:
                minimum, closest = measurement["angle_degrees"], [i, j]
            max_abs_dot = max(max_abs_dot, abs(measurement["dot"]))
    outputs.append({"kind": "packing_summary", "status": "DERIVED",
                    "min_angle_degrees": minimum, "closest_pair": closest,
                    "max_absolute_dot": max_abs_dot, "pair_count": k * (k - 1) // 2,
                    "derivation": "minimum angle and maximum absolute dot over every recorded unordered pair; not an optimality claim"})
    return _run(params, outputs)


def _trajectories(n: int, cfg: dict[str, Any]) -> Iterator[dict[str, Any]]:
    params = {"experiment": "random_walk", "dimension": n, "steps": cfg["steps"],
              "seed": cfg["seed"], "step_length": 1.0,
              "sampling": "independent normalized standard Gaussian increments"}
    rng = _rng(params)
    state = [0.0] * n
    outputs: list[dict[str, Any]] = [{"kind": "path_state", "step": 0, "state": state[:],
                                    "increment": None, "status": "OBSERVED"}]
    for step in range(1, cfg["steps"] + 1):
        increment = unit_vector(rng, n)
        state = [x + dx for x, dx in zip(state, increment)]
        outputs.append({"kind": "path_state", "step": step, "state": state[:],
                        "increment": increment, "status": "OBSERVED"})
    yield _run(params, outputs)
    params = {"experiment": "iterated_map", "dimension": n, "steps": cfg["steps"],
              "seed": cfg["seed"], "scale": cfg["map_scale"],
              "map": "x_next = scale*x + translation, coordinatewise",
              "sampling": "standard Gaussian initial point and fixed translation"}
    rng = _rng(params)
    state = [rng.gauss(0.0, 1.0) for _ in range(n)]
    translation = [rng.gauss(0.0, 1.0) for _ in range(n)]
    params["initial"] = state[:]
    params["translation"] = translation
    outputs = [{"kind": "path_state", "step": 0, "state": state[:], "status": "OBSERVED"}]
    for step in range(1, cfg["steps"] + 1):
        state = [cfg["map_scale"] * x + b for x, b in zip(state, translation)]
        outputs.append({"kind": "path_state", "step": step, "state": state[:], "status": "OBSERVED"})
    yield _run(params, outputs)


def iter_geometry_runs(config: Mapping[str, Any]) -> Iterator[dict[str, Any]]:
    """Yield every requested experiment/dimension run, with all its outputs.

    Finite stochastic samples are not exhaustive over the continuum. No sample,
    duplicate, repeated path point, or pairwise measurement is discarded.
    """
    cfg = validate_geometry_config(config)
    selected = EXPERIMENTS if cfg["experiment"] == "all" else (cfg["experiment"],)
    for n in cfg["dimensions"]:
        for experiment in selected:
            if experiment == "orthogonality":
                yield _orthogonality(n, cfg)
            elif experiment == "shell":
                yield _shell(n, cfg)
            elif experiment == "superposition":
                yield _superposition(n, cfg)
            else:
                yield from _trajectories(n, cfg)
