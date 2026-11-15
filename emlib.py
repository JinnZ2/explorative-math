"""emlib -- shared helpers for the experiments. CC0. Co-authored. Stdlib only.

- write_run(): writes runs/<name>_seed<k>.json and checks the schema. A run
  whose JSON is larger than GZ_THRESHOLD bytes is written as .json.gz
  instead: the same text, complete, gzip-compressed.
- load_run(), run_files(): read a run from either form; list both forms.
- lz76(), compressibility(): LZ76 phrase count and zlib ratio of a symbol string.
- vector helpers: re-exported from toolkit.primitives (one copy).
"""
import gzip
import json
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUNS = ROOT / "runs"

SCHEMA = ("params", "outputs", "angles", "human_interest", "ai_interest",
          "compressibility", "nature_match", "divergence_note")
NOT_INTERPRETABLE = "not interpretable atm"
GZ_THRESHOLD = 256 * 1024   # bytes of JSON text above which a run is written gzipped


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
              compressibility, nature_match, out_dir=None, gz=None):
    """gz: True / False forces the form; None chooses by GZ_THRESHOLD. Either
    way the full record is written, and a stale file of the other form for
    the same run is removed, so a run has one file."""
    rec = {"params": dict(params, seed=seed), "outputs": outputs, "angles": angles,
           "human_interest": None, "ai_interest": ai_interest,
           "compressibility": compressibility, "nature_match": nature_match,
           "divergence_note": None}
    check_run(rec)
    d = Path(out_dir) if out_dir else RUNS
    d.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(rec, indent=1, sort_keys=False) + "\n").encode("utf-8")
    plain = d / f"{name}_seed{seed}.json"
    packed = d / f"{name}_seed{seed}.json.gz"
    if gz if gz is not None else len(data) > GZ_THRESHOLD:
        with packed.open("wb") as fh:
            # mtime 0 and no stored name: the same record gives the same bytes
            with gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0) as z:
                z.write(data)
        plain.unlink(missing_ok=True)
        return packed
    plain.write_bytes(data)
    packed.unlink(missing_ok=True)
    return plain


def load_run(path):
    """Read a run from .json or .json.gz. If the named form is absent and the
    other form is there, the other form is read."""
    p = Path(path)
    if not p.exists():
        other = (p.with_name(p.name[:-3]) if p.name.endswith(".json.gz")
                 else p.with_name(p.name + ".gz"))
        if not other.exists():
            raise FileNotFoundError(str(p))
        p = other
    if p.name.endswith(".gz"):
        with gzip.open(p, "rt", encoding="utf-8") as fh:
            return json.load(fh)
    return json.loads(p.read_text(encoding="utf-8"))


def run_files(runs_dir=None):
    """Every run file in a directory, both forms, sorted by name."""
    d = Path(runs_dir) if runs_dir else RUNS
    return sorted(list(d.glob("*.json")) + list(d.glob("*.json.gz")))


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
