"""emlib -- shared helpers for the experiments. CC0. Co-authored. Stdlib only.

- write_run(): writes runs/<name>_seed<k>.json and checks the schema.
- lz76(), compressibility(): LZ76 phrase count and zlib ratio of a symbol string.
- vector helpers: re-exported from toolkit.primitives (one copy).
"""
import json
import zlib
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
    """LZ76 phrase count and zlib ratio of the same symbol string."""
    s = "".join(symbols)
    return {"definition": definition, "n": len(s), "phrases": lz76(s),
            "zlib_ratio": (len(zlib.compress(s.encode("utf-8"), 9)) / len(s)) if s else None}


# vector helpers live in toolkit.primitives; re-exported so the experiments keep working
from toolkit.primitives import dot, norm, unit, random_unit, angle_deg  # noqa: E402,F401
