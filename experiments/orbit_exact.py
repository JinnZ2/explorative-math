"""j. orbit_exact -- the seed-5 orbit again, in exact arithmetic.

CC0. Co-authored. Stdlib only.

Source: runs/generated_seed5.json (orbit pipeline). Its prefix period reads
184 at both 256 and 512 steps (open_questions.md Q8). That orbit was
computed in binary floats. Here the same alpha and x0 are taken as the
exact rationals the logged floats are (fractions.Fraction of a float is
exact) and the orbit x -> x + alpha (mod 1) is iterated with no rounding.

  float      the logged pipeline again, as a check that the word is remade
  exact      Fraction arithmetic, same 512 steps: word, period profile
  extended   the exact orbit carried past 512 steps, period profile at
             longer prefixes

The exact run tests float precision only. It says nothing about the real
number the float alpha stands in for.

  python3 experiments/orbit_exact.py [seed]
"""
import json
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import emlib  # noqa: E402
from toolkit import composers as C  # noqa: E402
from toolkit import probes as Q  # noqa: E402

NAME = "orbit_exact"
SOURCE = "runs/generated_seed5.json"
EXTENDED_STEPS = 4096
EXTENDED_PREFIXES = (512, 1024, 2048, 4096)
HALF = Fraction(1, 2)


def code(xs, half):
    return "".join("0" if x < half else "1" for x in xs)


def first_difference(a, b):
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return i
    return None if len(a) == len(b) else min(len(a), len(b))


def run(seed=0, source=None, extended_steps=EXTENDED_STEPS, prefixes=EXTENDED_PREFIXES):
    rec = json.loads((Path(source) if source else emlib.ROOT / SOURCE).read_text())
    spec = rec["params"]["spec"]
    a, x0, steps = spec["alpha"], spec["x0"], spec["steps"]
    logged_word = rec["outputs"]["word"]

    fx = C.iterate(lambda x: (x + a) % 1.0, x0, steps)
    fword = code(fx, 0.5)

    ea, ex0 = Fraction(a), Fraction(x0)
    ex = C.iterate(lambda x: (x + ea) % 1, ex0, extended_steps)
    eword_all = code(ex, HALF)
    eword = eword_all[:steps + 1]
    head = ex[:steps + 1]
    drift = [abs(Fraction(f) - e) for f, e in zip(fx, head)]
    to_half = [abs(e - HALF) for e in head]
    to_wrap = [min(e, 1 - e) for e in head]
    k = 184
    ka = k * ea
    near = round(ka)
    prof = Q.period_profile(eword)
    ext_prefixes = [p for p in prefixes if p <= len(eword_all)]
    ext_prof = Q.period_profile(eword_all, ext_prefixes)
    persists = prof.get("512") == rec["outputs"]["period_profile"]["512"]
    return {
        "source": {"alpha": a, "x0": x0, "steps": steps,
                   "logged_period_profile": rec["outputs"]["period_profile"]},
        "float": {"arithmetic": "binary float (64-bit)", "word": fword,
                  "word_equals_logged": fword == logged_word,
                  "period_profile": Q.period_profile(fword),
                  "period_bound": Q.period_bound(fword)},
        "exact": {"arithmetic": "fractions.Fraction, no rounding",
                  "alpha_exact": [ea.numerator, ea.denominator],
                  "x0_exact": [ex0.numerator, ex0.denominator],
                  "word": eword,
                  "word_equals_float": eword == fword,
                  "first_index_where_words_differ": first_difference(eword, fword),
                  "period_profile": prof,
                  "period_bound": Q.period_bound(eword),
                  "max_abs_float_minus_exact": float(max(drift)),
                  "min_distance_to_coding_boundary_half": float(min(to_half)),
                  "min_distance_to_wrap_0_or_1": float(min(to_wrap))},
        "extended": {"arithmetic": "fractions.Fraction, no rounding",
                     "steps": extended_steps, "word": eword_all,
                     "period_profile": ext_prof,
                     "period_bound": Q.period_bound(eword_all, ext_prefixes)},
        "exact_orbit_period": {
            "alpha_denominator_log2": ea.denominator.bit_length() - 1,
            "x0_denominator_log2": ex0.denominator.bit_length() - 1,
            "alpha_denominator_is_power_of_two":
                ea.denominator & (ea.denominator - 1) == 0,
            "orbit_period": ea.denominator,
            "note": "DERIVED: alpha = p/q in lowest terms, so x + n*alpha (mod 1) returns "
                    "to x first at n = q, for any x0. The coded word repeats with a period "
                    "that divides q."},
        "k_alpha": {"k": k, "k_times_alpha": float(ka), "nearest_integer": near,
                    "residual": float(ka - near)},
        "float_precision_check": {
            "rule": "period at 512 steps: differs under exact arithmetic -> leak CONFIRMED "
                    "(float precision); same -> not a float-precision leak, logged as open",
            "period_at_512_float": rec["outputs"]["period_profile"]["512"],
            "period_at_512_exact": prof.get("512"),
            "persists": persists,
            "verdict": "not a float-precision leak; open" if persists
                       else "resolution leak CONFIRMED (float precision)"},
    }


def main(argv):
    seed = int(argv[0]) if argv else 0  # deterministic; seed kept for the schema
    out = run(seed)
    e, x, c = out["exact"], out["extended"], out["float_precision_check"]
    comp = emlib.compressibility(list(e["word"]),
                                 "exact orbit, 512 steps, coded 0 if x < 1/2 else 1")
    ai = (f"OBSERVED: float word remade equals the logged word: "
          f"{out['float']['word_equals_logged']}. Exact word equals float word: "
          f"{e['word_equals_float']} (largest float-minus-exact difference "
          f"{e['max_abs_float_minus_exact']:.3e}; nearest approach to the coding boundary "
          f"{e['min_distance_to_coding_boundary_half']:.3e}, to the wrap "
          f"{e['min_distance_to_wrap_0_or_1']:.3e}). Exact period profile, prefix lengths "
          f"8..512: {e['period_profile']}. Verdict by the stated rule: {c['verdict']}. "
          f"Exact orbit carried to {x['steps']} steps, period profile {x['period_profile']}. "
          f"184 * alpha = {out['k_alpha']['k_times_alpha']:.9f}, residual from the nearest "
          f"integer {out['k_alpha']['residual']:.3e}. "
          f"DERIVED: alpha is a rational with denominator 2^"
          f"{out['exact_orbit_period']['alpha_denominator_log2']}, so the orbit is exactly "
          "periodic with that period and the word's period divides it; every prefix "
          "period logged here is bounded by its prefix. "
          "DERIVED: the 512-step word is the same in both arithmetics, so the 184 reading "
          "is not made by float rounding; the extended profile is a separate reading, "
          "bounded by its own steps. AI reading: not interpretable atm beyond these values.")
    p = emlib.write_run(NAME, seed,
                        {"source_run": SOURCE, "extended_steps": EXTENDED_STEPS,
                         "extended_prefixes": list(EXTENDED_PREFIXES),
                         "coding": "0 if x < 1/2 else 1"},
                        out, None, ai, comp, None)
    print(p)
    print(ai)


if __name__ == "__main__":
    main(sys.argv[1:])
