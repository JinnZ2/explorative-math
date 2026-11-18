# Opt-in options on four probes, 2026-11-18

CC0. Co-authored. Written by Manus (an AI agent) at the request of the repository owner, as a proposal.
Status words follow `CLAUDE.md`: OBSERVED = read from a run or a test, DERIVED = follows from shown values,
PROPOSED = not decided.

## What this changes

Four keyword options, all **off by default**. With them off, every function returns what it returned at commit
`3162c99`: the same keys and the same values. `tests/test_probe_optin.py` holds the old bodies as frozen copies and
compares thousands of random inputs against them. The mutation checks (leaking an option key into a default result,
turning an option on by default) each fail those tests. No stored run, experiment or caller is edited, so nothing
already committed changes.

| option | answers | what it adds | what it does not change |
|---|---|---|---|
| `three_gap_sum(..., rel_tol=None)` | A01, A02 | with `rel_tol`, passes when `abs(residual) <= rel_tol * L3`, and records `rel_tol` and `tol_used` | the absolute `tol` test when `rel_tol` is None |
| `angle_histogram(..., account=False)` | A04, A05 | `n_input`, `n_binned`, `below_lo`, `above_hi`, `bins_cover`, `binned_above_hi` | every count; angles outside the range are still not counted, they are now counted aside |
| `rotational_symmetry_bound(..., with_spacing=False)` | A07 | `min_pair_distance` and `n_points_within_tol` | the symmetry list; it is still read at the same tolerance |
| `emlib.write_run(..., strict_json=False)`, `emlib.nonfinite_paths` | A19 | refuses a record with NaN or an infinite float, before any file is written, naming where they are | the default writer; it still writes bare `NaN` / `Infinity` |

## Why each one, and against what

The finding labels are those of `docs/probe_provenance_2026-11-18.json` (PR #3). "Audited" classes are a post-hoc
revision by the author of the predictions; the quotes are there so that you can disagree with them.

**`rel_tol` (A01, A02): LIMIT, independent truth.** The docstring says only that the function "Returns the verdict
"pass" or "fail" with the residual L3 - (L1 + L2)." It does not say whether `tol` is absolute. The three-distance
identity is scale-free, so an absolute 1e-6 is a unit-dependent bound. OBSERVED in the tests: lengths
(1e-9, 2e-9, 7e-9) pass the absolute test with a residual of 57 % of the largest length, and fail with
`rel_tol=1e-6`; the golden triple times 1e12, each value rounded to a float on its own, fails the absolute test on one
float step (4.9e-4) and passes with `rel_tol=1e-9`. Both callers in the repo use lengths of order 1.

**`account` (A04, A05): A05 is a DEFECT on the function's own output; A04 is UNDEFINED.** `angle_histogram` has no
docstring. A05: for width 7 on [0, 180] the record says `hi = 180` and has 26 bins reaching 182, and counts 181;
`bins_cover` now says `[0.0, 182.0]` and `binned_above_hi` says 1. A04: angles outside [lo, hi] vanish without a count;
"keep every requested output" (`CONTRIBUTING.md`) reaches that only if they were requested, so the option records
the count and leaves the decision. Only caller: `toolkit/generator.py:90`, width 15 on [0, 360], where neither
case arises.

**`with_spacing` (A07): LIMIT, repo rule met in part.** `CLAUDE.md`: "any count / period / symmetry reading logs the
instrument bound it was read at". The bound line does log `abs_tol`, and `degeneracy()` already tags the A07 set
`coincident` and `central`. What the bound line lacks is the comparison that makes `abs_tol` readable: on the A07 set
`abs_tol = 1.000000004553418` against a smallest point distance of about 0.45, with 3 points within tolerance of another.
The reading `[2, 3, 4, 6, 12]` is unchanged.

**`strict_json` (A19): DEFECT against RFC 8259, no direct repo rule for this build.** RFC 8259 section 6: "Numeric
values that cannot be represented in the grammar below (such as Infinity and NaN) are not permitted." Section 10:
"The resulting text MUST strictly conform to the JSON grammar." `docs/schema.md` rejects nonfinite numbers, but it
governs the other build. Seven committed files in `runs/` hold 191 bare `+Infinity` values (`combine_seed0`,
`combine_seed1`, `combine_small_seed0`, `fold_loop_seed0`, `fold_shape_seed0`, `layered_pursuit_seed0`,
`layered_pursuit_small_seed0`), written where `fold_loop.py:180`, `fold_shape.py:199` and `layered_pursuit.py:78` map
"no result" to `float("inf")`. Python's `json` reads them back; strict parsers do not. That is why the default stays as it
is: turning it on would stop those experiments from rerunning.

## What is not done, and why

| not done | reason |
|---|---|
| default changes | each would alter output for some input; whether that is wanted is the maintainers' decision. A05 alone would alter no current call |
| `polar_angles` (A06) | no docstring; what a point at the centroid should return is a design choice |
| first-fit matching (A09), mirror grid (A10) | LIMITs that the repo's own bound functions describe; a different algorithm is a larger change than an option |
| the seven committed files with `+Infinity` | recorded values; changing `inf` to `null` rewrites results. `CONTRIBUTING.md`: "Unknown is `null`" points that way, but only for observer slots |
| turning the options on in `generator.py` or any experiment | would add keys to new runs; not asked for |

## Reproduce

```text
python -m unittest discover -s tests -p "test_probe_optin.py" -v
python -m unittest discover -s tests
```
