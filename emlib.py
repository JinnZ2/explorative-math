"""emlib -- shared helpers for the experiments. CC0. Co-authored. Stdlib only.

- write_run(): writes runs/<name>_seed<k>.json and checks the schema.
- lz76(): Lempel-Ziv 1976 phrase count of a symbol string (compressibility).
- vector helpers: dot, norm, unit, random_unit, angle_deg.
"""
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUNS = ROOT / "runs"

SCHEMA = ("params", "outputs", "angles", "human_interest", "ai_interest",
          "compressibility", "nature_match", "divergence_note")
NOT_INTERPRETABLE = "not interpretable atm"


class SchemaError(ValueError):
    pass


def check_run(rec):
    missing = [k for k in SCHEMA if k not in rec]
    extra = [k for k in rec if k not in SCHEMA]
    if missing or extra:
        raise SchemaError(f"missing {missing} extra {extra}")
    if rec["human_interest"] is not None:
        raise SchemaError("human_interest is the human's to fill; must be null when written")
    if rec["ai_interest"] is None:
        raise SchemaError("ai_interest must be a reading or NOT_INTERPRETABLE, not null")
    c = rec["compressibility"]
    if c is not None and not {"definition", "n", "phrases"} <= set(c):
        raise SchemaError("compressibility needs definition, n, phrases (or null)")
    return rec


def write_run(name, seed, params, outputs, angles, ai_interest,
              compressibility, nature_match, out_dir=None):
    rec = {"params": dict(params, seed=seed), "outputs": outputs, "angles": angles,
           "human_interest": None, "ai_interest": ai_interest,
           "compressibility": compressibility, "nature_match": nature_match,
           "divergence_note": None}
    check_run(rec)
    d = Path(out_dir) if out_dir else RUNS
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"{name}_seed{seed}.json"
    p.write_text(json.dumps(rec, indent=1, sort_keys=False) + "\n", encoding="utf-8")
    return p


def lz76(s):
    """Number of phrases in the Lempel-Ziv (1976) parsing of s.

    Each new phrase is the shortest substring starting at i that has not
    appeared in s[:i+len-1]. Low count relative to n = more compressible."""
    n = len(s)
    if n == 0:
        return 0
    i, c = 0, 0
    while i < n:
        k = 1
        while i + k <= n and s[i:i + k] in s[:i + k - 1]:
            k += 1
        c += 1
        i += k
    return c


def compressibility(symbols, definition):
    s = "".join(symbols)
    return {"definition": definition, "n": len(s), "phrases": lz76(s)}


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    n = norm(a)
    return [x / n for x in a]


def random_unit(rng, d):
    while True:
        v = [rng.gauss(0.0, 1.0) for _ in range(d)]
        if norm(v) > 0:
            return unit(v)


def angle_deg(a, b):
    c = dot(a, b) / (norm(a) * norm(b))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))
