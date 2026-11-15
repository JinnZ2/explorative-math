"""Finite, deterministic, unfiltered mathematical-space enumeration.

These iterators retain repeated states and duplicate seeds.  They do not judge
interestingness, compress trajectories, or claim to enumerate an infinite space.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from itertools import product
import hashlib
import json
from typing import Any, Iterator, Mapping, Sequence

SCHEMA_VERSION = "1.0"
GENERATOR_VERSION = "1.0.0"
STATUSES = frozenset({"OBSERVED", "DERIVED", "PROPOSED"})
OBSERVER_SLOTS = ("human_interest", "ai_interest", "compressibility", "divergence_note")
FAMILIES = ("ca", "substitution", "integer")


class ValidationError(ValueError):
    """A configuration, dataset, or explicit input violates the contract."""


def require_int(name: str, value: Any, minimum: int | None = None) -> int:
    """Accept integers, never booleans or floating-point lookalikes."""
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValidationError(f"{name} must be an integer (not a boolean)")
    if minimum is not None and value < minimum:
        raise ValidationError(f"{name} must be at least {minimum}")
    return value


def require_binary_word(name: str, word: Any) -> str:
    if not isinstance(word, str) or not word or any(c not in "01" for c in word):
        raise ValidationError(f"{name} must be a nonempty binary string")
    return word


def canonical_json(value: Any) -> str:
    """Canonical UTF-8 JSON text: sorted keys, compact separators, no NaN."""
    return json.dumps(value, ensure_ascii=False, sort_keys=True,
                      separators=(",", ":"), allow_nan=False)


def stable_id(prefix: str, value: Any) -> str:
    digest = hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()
    return f"{prefix}-{digest}"


@dataclass(frozen=True)
class EnumerationConfig:
    """Explicit finite bounds; all seeds are evaluated without deduplication.

    CA single mode is one 1 at index ``width // 2``, with all other cells 0.
    Substitution seeds may repeat intentionally.  Integer seeds are 0..m-1.
    Each run contains its initial state plus exactly ``steps`` updates.
    """

    family: str = "all"
    steps: int = 4
    width: int = 5
    ca_seeds: str = "all"
    max_word_length: int = 2
    substitution_seeds: tuple[str, ...] = ("0", "1")
    a_min: int = -1
    a_max: int = 1
    b_min: int = -1
    b_max: int = 1
    modulus: int = 8

    def __post_init__(self) -> None:
        if self.family not in (*FAMILIES, "all"):
            raise ValidationError("family must be ca, substitution, integer, or all")
        require_int("steps", self.steps, 0)
        require_int("width", self.width, 1)
        require_int("max_word_length", self.max_word_length, 1)
        require_int("modulus", self.modulus, 1)
        for name in ("a_min", "a_max", "b_min", "b_max"):
            require_int(name, getattr(self, name))
        if self.a_min > self.a_max or self.b_min > self.b_max:
            raise ValidationError("inclusive lower bounds must not exceed upper bounds")
        if self.ca_seeds not in ("all", "single"):
            raise ValidationError("ca_seeds must be all or single")
        if not isinstance(self.substitution_seeds, (tuple, list)) or not self.substitution_seeds:
            raise ValidationError("substitution_seeds must be a nonempty sequence")
        for word in self.substitution_seeds:
            require_binary_word("substitution seed", word)
        object.__setattr__(self, "substitution_seeds", tuple(self.substitution_seeds))

    @property
    def families(self) -> tuple[str, ...]:
        return FAMILIES if self.family == "all" else (self.family,)

    def to_dict(self) -> dict[str, Any]:
        result = {field.name: getattr(self, field.name) for field in fields(self)}
        result["substitution_seeds"] = list(self.substitution_seeds)
        return result

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> EnumerationConfig:
        if not isinstance(value, Mapping):
            raise ValidationError("config must be an object")
        names = {field.name for field in fields(cls)}
        if set(value) != names:
            raise ValidationError("config must contain exactly the declared configuration keys")
        return cls(**dict(value))

    @property
    def dataset_id(self) -> str:
        return stable_id("dataset", {"schema_version": SCHEMA_VERSION,
                                     "generator_version": GENERATOR_VERSION,
                                     "config": self.to_dict()})

    def expected_counts(self) -> dict[str, int]:
        runs = 0
        if "ca" in self.families:
            runs += 256 * (2 ** self.width if self.ca_seeds == "all" else 1)
        if "substitution" in self.families:
            choices = 2 ** (self.max_word_length + 1) - 2
            runs += choices * choices * len(self.substitution_seeds)
        if "integer" in self.families:
            runs += (self.a_max - self.a_min + 1) * (self.b_max - self.b_min + 1) * self.modulus
        return {"runs": runs, "outputs": runs * (self.steps + 1)}


def iter_ca_rules() -> Iterator[int]:
    """ECA rule numbers in ascending order, 0 through 255 inclusive."""
    yield from range(256)


def iter_ca_seeds(width: int = 5, mode: str = "all") -> Iterator[list[int]]:
    """Lexicographic binary rows, or a centered single-cell seed."""
    require_int("width", width, 1)
    if mode == "single":
        row = [0] * width
        row[width // 2] = 1
        yield row
    elif mode == "all":
        for bits in product((0, 1), repeat=width):
            yield list(bits)
    else:
        raise ValidationError("CA seed mode must be all or single")


def _binary_row(row: Sequence[int]) -> list[int]:
    if not isinstance(row, (list, tuple)) or not row:
        raise ValidationError("CA row must be a nonempty list or tuple")
    for bit in row:
        require_int("CA cell", bit, 0)
        if bit > 1:
            raise ValidationError("CA cells must be 0 or 1")
    return list(row)


def evolve_ca(row: Sequence[int], rule: int) -> list[int]:
    """One synchronous Wolfram ECA step with circular boundary conditions."""
    require_int("ECA rule", rule, 0)
    if rule > 255:
        raise ValidationError("ECA rule must be at most 255")
    state = _binary_row(row)
    width = len(state)
    return [(rule >> ((state[(i - 1) % width] << 2) |
                      (state[i] << 1) | state[(i + 1) % width])) & 1
            for i in range(width)]


def iter_ca_trajectory(rule: int, seed: Sequence[int], steps: int = 4) -> Iterator[list[int]]:
    """Yield step 0 and every update, including all repeats."""
    require_int("steps", steps, 0)
    require_int("ECA rule", rule, 0)
    if rule > 255:
        raise ValidationError("ECA rule must be at most 255")
    state = _binary_row(seed)
    yield state.copy()
    for _ in range(steps):
        state = evolve_ca(state, rule)
        yield state.copy()


def iter_binary_words(max_word_length: int = 2) -> Iterator[str]:
    """Nonempty binary words, shortest first, then lexicographic."""
    require_int("max_word_length", max_word_length, 1)
    for length in range(1, max_word_length + 1):
        for bits in product("01", repeat=length):
            yield "".join(bits)


def iter_substitution_rules(max_word_length: int = 2) -> Iterator[dict[str, str]]:
    """All ordered RHS pairs; the 0 RHS is the outer loop, 1 the inner."""
    require_int("max_word_length", max_word_length, 1)
    # Reiterate the inner words instead of materializing the RHS space.
    for zero in iter_binary_words(max_word_length):
        for one in iter_binary_words(max_word_length):
            yield {"0": zero, "1": one}


def _substitution_rule(rule: Mapping[str, str]) -> dict[str, str]:
    if not isinstance(rule, Mapping) or set(rule) != {"0", "1"}:
        raise ValidationError("substitution rule must map exactly 0 and 1")
    return {symbol: require_binary_word(f"RHS for {symbol}", rule[symbol])
            for symbol in ("0", "1")}


def substitute(state: str, rule: Mapping[str, str]) -> str:
    """Apply a binary substitution simultaneously to every symbol."""
    require_binary_word("substitution state", state)
    checked = _substitution_rule(rule)
    return "".join(checked[symbol] for symbol in state)


def iter_substitution_trajectory(rule: Mapping[str, str], seed: str,
                                 steps: int = 4) -> Iterator[str]:
    require_int("steps", steps, 0)
    checked = _substitution_rule(rule)
    state = require_binary_word("substitution seed", seed)
    yield state
    for _ in range(steps):
        state = "".join(checked[symbol] for symbol in state)
        yield state


def iter_integer_rules(a_min: int = -1, a_max: int = 1, b_min: int = -1,
                       b_max: int = 1, modulus: int = 8) -> Iterator[dict[str, int]]:
    """Inclusive ascending a then b; syntactically different rules stay distinct."""
    for name, value in (("a_min", a_min), ("a_max", a_max),
                        ("b_min", b_min), ("b_max", b_max)):
        require_int(name, value)
    require_int("modulus", modulus, 1)
    if a_min > a_max or b_min > b_max:
        raise ValidationError("inclusive lower bounds must not exceed upper bounds")
    for a in range(a_min, a_max + 1):
        for b in range(b_min, b_max + 1):
            yield {"a": a, "b": b, "modulus": modulus}


def iter_integer_seeds(modulus: int = 8) -> Iterator[int]:
    require_int("modulus", modulus, 1)
    yield from range(modulus)


def evolve_integer(value: int, a: int, b: int, modulus: int) -> int:
    require_int("value", value, 0)
    require_int("a", a)
    require_int("b", b)
    require_int("modulus", modulus, 1)
    if value >= modulus:
        raise ValidationError("integer state must be less than modulus")
    return (a * value + b) % modulus


def iter_integer_trajectory(a: int, b: int, modulus: int, seed: int,
                            steps: int = 4) -> Iterator[int]:
    require_int("steps", steps, 0)
    # Validate even a step-0-only run.
    evolve_integer(seed, a, b, modulus)
    value = seed
    yield value
    for _ in range(steps):
        value = (a * value + b) % modulus
        yield value


def iter_runs(config: EnumerationConfig) -> Iterator[tuple[str, dict[str, Any], Any, Iterator[Any]]]:
    """Yield (family, rule parameters, seed, trajectory iterator) in fixed order."""
    if not isinstance(config, EnumerationConfig):
        raise ValidationError("config must be an EnumerationConfig")
    if "ca" in config.families:
        for rule in iter_ca_rules():
            for seed in iter_ca_seeds(config.width, config.ca_seeds):
                params = {"rule_number": rule, "boundary": "circular", "width": config.width}
                yield "ca", params, seed, iter_ca_trajectory(rule, seed, config.steps)
    if "substitution" in config.families:
        for rule in iter_substitution_rules(config.max_word_length):
            for seed in config.substitution_seeds:
                yield "substitution", rule, seed, iter_substitution_trajectory(rule, seed, config.steps)
    if "integer" in config.families:
        for rule in iter_integer_rules(config.a_min, config.a_max, config.b_min, config.b_max,
                                       config.modulus):
            for seed in iter_integer_seeds(config.modulus):
                yield "integer", rule, seed, iter_integer_trajectory(
                    rule["a"], rule["b"], rule["modulus"], seed, config.steps)


def iter_records(config: EnumerationConfig | None = None) -> Iterator[dict[str, Any]]:
    """Stream one OBSERVED record per run with every state in ``outputs``.

    Run ordinals distinguish duplicate seeds and modulo-equivalent rules. IDs
    are deterministic across output paths, timestamps, and repeated executions.
    Each yielded object owns its seed/rule/state containers. All explicit
    observer slots are top-level and remain null in the base archive.
    """
    config = EnumerationConfig() if config is None else config
    if not isinstance(config, EnumerationConfig):
        raise ValidationError("config must be an EnumerationConfig")
    namespace = config.dataset_id
    for ordinal, (family, rule, seed, trajectory) in enumerate(iter_runs(config)):
        run_id = stable_id("run", {"dataset_id": namespace, "ordinal": ordinal,
                                   "family": family, "rule": rule, "seed": seed})
        outputs = [{"step": step, "state": list(state) if family == "ca" else state}
                   for step, state in enumerate(trajectory)]
        yield {"id": run_id, "run_id": run_id, "family": family, "kind": "trajectory",
               "params": {"rule": dict(rule),
                          "seed": list(seed) if family == "ca" else seed,
                          "steps": config.steps},
               "outputs": outputs, "status": "OBSERVED", **dict.fromkeys(OBSERVER_SLOTS)}
