# Adversarial probe test and Q5 slopes, 2026-11-18

CC0. Co-authored. Written by Manus (an AI agent) at the request of the repository owner. Status words follow
`CLAUDE.md`: OBSERVED = read from a run file, DERIVED = follows from shown values, PROPOSED = not tested.
Nothing in `toolkit/`, `explorative_math/` or `emlib.py` was changed. Both interest slots in every new run are null.

Predictions were written and hashed before any run (`docs/prereg/`):

| file | sha256 |
|---|---|
| `prereg_q12_q11_q5_q9_2026-11-18.json` | `6ae97bcce1f52ae6f5c24f385f5fb16125b2f22b5ac5a4e0ca1dbfa6bffa606c` |
| `prereg_adversarial_q5slopes_2026-11-18.json` | `e31a31a1dc5964c2db8dd8bf11ae3ce66fe230b009d7dbb9594f8b56eee12bdd` |

## 1. Adversarial probe test (`experiments/adversarial_probes.py`, `runs/adversarial_probes_seed0.json`)

20 cases against `toolkit/probes.py` and `emlib.py`. Each case has a spec (a docstring, a repo rule, or a ground
truth computed by an independent method) and a prediction. Classes: HOLDS, DEFECT (contradicts docstring, repo rule
or independent truth on in-scope input), LIMIT (misbehaves outside stated scope or in a regime the bound line
names), UNDEFINED (spec silent; behaviour recorded). OBSERVED (pre-registered classes): 5 HOLDS, 7 DEFECT, 4 LIMIT, 4 UNDEFINED; 19 of 20 predicted classes were right,
the miss is A15 (below). Section 1a revises the DEFECT count to 2 after checking what each class rests on.

| id | probe | class | what was read |
|---|---|---|---|
| A01 | `three_gap_sum` | DEFECT | (1e-9, 2e-9, 7e-9) passes: the tolerance is absolute 1e-6 |
| A02 | `three_gap_sum` | DEFECT | the golden triple (1, phi, phi^2) x 1e12 fails: float spacing there exceeds 1e-6 (residual 4.9e-4) |
| A03 | `three_gap_sum` | UNDEFINED | NaN gives fail, [1,1,2] and [0,0,0] give pass; none refused |
| A04 | `angle_histogram` | DEFECT | 3 of 4 angles outside [lo, hi] vanish; no count of them is returned (conflicts with the retain-every-output rule in `CONTRIBUTING.md`) |
| A05 | `angle_histogram` | DEFECT | width 7 on [0,180] makes 26 bins reaching 182; 181 is counted |
| A06 | `polar_angles` | DEFECT | a point at the centroid is given direction 0.0; an angle of -1e-15 deg is returned as 360.0 |
| A07 | `rotational_symmetry` | DEFECT | an asymmetric 3-point cluster beside a 12-point ring of radius 1e6 reads [2,3,4,6,12]: abs_tol = 1.0 exceeds the cluster (the bound line does log abs_tol) |
| A08 | `rotational_symmetry` | HOLDS | a vertex pushed out radially flips the reading at eps = 1.0104e-6, i.e. at rel_tol |
| A09 | `rotational_symmetry` | LIMIT | greedy first-fit `_matches` vs exact bipartite matching, 20,000 near-coincident 2-fold trials: exact matching finds symmetry in 46, the probe misses 8 of them, and never claims symmetry that does not exist |
| A10 | `mirror_symmetry` vs `_from_pairs` | LIMIT | pentagon at 500 random rotations: grid probe finds 0 axes (0.5 deg step), pair probe finds 5 axes 500 of 500; 21 on-grid controls all find 5 |
| A11 | `rotational_symmetry` | LIMIT | hexagon translated by 1e9 still reads [2,3,6]; by 1e12 reads [2] only (my prediction said [] at 1e10; it was wrong in detail) |
| A12 | several | UNDEFINED | [] raises IndexError in `mirror_symmetry` and `degeneracy`; identical points read symmetric at every order 2..12; NaN and inf are accepted silently; a 3D cube about z reads [2,4]; tilted 3D cube reads []; mirror and degeneracy raise ValueError on 3D |
| A13 | `rotational_symmetry` | UNDEFINED | a hexagon with one vertex listed twice reads [] (centroid shifts); the same hexagon as a set reads [2,3,6]; no note |
| A14 | `emlib.lz76` | HOLDS | 0 mismatches with an independent Kaspar-Schuster algorithm over 3,000 random strings plus structured and known values |
| A15 | `emlib.lz76` | HOLDS | predicted too slow; n=20,000 takes 0.24 s here (growth exponent about 1.8). Prediction wrong |
| A16 | `emlib.compressibility` | UNDEFINED | symbols ["1","10"] and ["11","0"] give identical readings (tokens are joined first) |
| A17 | `smallest_period` | HOLDS | 0 mismatches with n minus longest proper border, 3,000 strings |
| A18 | `min_stretch_tag` | LIMIT | frequencies -0.5, 5.0 and NaN are accepted and assessed |
| A19 | `emlib.write_run` | DEFECT | NaN and Infinity are written as bare tokens; a strict JSON parser rejects the file |
| A20 | `write_run` / `check_run` | HOLDS | non-null human_interest, missing key, extra key all refused; same record twice gives identical gz bytes |

The run file for this experiment replaces NaN and inf by strings so that it is itself strict JSON.

## 1a. What each class is measured against (added after the first write-up)

The first table above gave one class per case without saying what it was measured against. A DEFECT can mean
four different things, and they carry different weight:

1. **Docstring**: the function's own words are contradicted.
2. **Repo rule**: a sentence in `CONTRIBUTING.md`, `CLAUDE.md` or `docs/schema.md` is contradicted. Each of those files governs a named part of the repo (`docs/schema.md` governs the Manus build, not `emlib.py`), and a rule can apply directly or by analogy.
3. **Independent truth**: a result computed outside the repo (an exact matching, an IEEE 754 spacing, RFC 8259) disagrees.
4. **Own output**: the function's returned record disagrees with itself.

Each case now has these recorded in `docs/probe_provenance_2026-11-18.json`, with the verbatim quote, its source,
and a strength (DIRECT, ANALOGY, EXTERNAL, COMPUTED). `tests/test_probe_provenance.py` checks that every repo quote
occurs in the file it cites. The run file keeps both classes: `class` is the pre-registered reading, `audited_class`
a **post-hoc revision by the same author**, with its reason in `audit_note`. I am an interested party; the quotes are
there so that you can disagree with a class without taking my word for it.

Result of the audit: 5 HOLDS, 2 DEFECT, 7 LIMIT, 6 UNDEFINED. Five of the seven first-run DEFECTs moved down.

| id | prereg class | audited class | measured against |
|---|---|---|---|
| A01 | DEFECT | LIMIT | docstring, independent truth |
| A02 | DEFECT | LIMIT | docstring, independent truth |
| A03 | UNDEFINED | UNDEFINED | docstring |
| A04 | DEFECT | UNDEFINED | no spec, repo rule (analogy) |
| A05 | DEFECT | DEFECT | own output |
| A06 | DEFECT | UNDEFINED | independent truth, no spec, repo rule (analogy) |
| A07 | DEFECT | LIMIT | docstring, repo rule |
| A08 | HOLDS | HOLDS | derived here |
| A09 | LIMIT | LIMIT | docstring, independent truth |
| A10 | LIMIT | LIMIT | docstring, repo rule |
| A11 | LIMIT | LIMIT | independent truth, repo rule |
| A12 | UNDEFINED | UNDEFINED | no spec, repo rule (analogy) |
| A13 | UNDEFINED | UNDEFINED | docstring (analogy), repo rule (analogy) |
| A14 | HOLDS | HOLDS | docstring, independent truth |
| A15 | HOLDS | HOLDS | no spec |
| A16 | UNDEFINED | UNDEFINED | docstring |
| A17 | HOLDS | HOLDS | independent truth |
| A18 | LIMIT | LIMIT | docstring |
| A19 | DEFECT | DEFECT | independent truth (external), repo rule (analogy) |
| A20 | HOLDS | HOLDS | repo rule |

What moved and why:

- **A01, A02 (three_gap_sum)**: DEFECT to LIMIT. The docstring does not say the tolerance is absolute or relative, and the output records `tol` and `residual`. The independent truth (the identity L3 = L1 + L2 is scale-free; one float step at 2.6e12 is 4.9e-4) shows the absolute 1e-6 is a unit-dependent instrument bound. No text is contradicted.
- **A04 (angle_histogram drops angles)**: DEFECT to UNDEFINED. There is no docstring. "keep every requested output" (`CONTRIBUTING.md`) reaches it only if an angle outside the range counts as requested. That is a reading, not the text. The one caller feeds angles in [0, 360).
- **A06 (polar_angles)**: DEFECT to UNDEFINED. No docstring. The direction of a zero vector is undefined (independent truth), and `-1e-15 % 360.0 == 360.0` is float behaviour that the one caller absorbs through the `a == hi` branch.
- **A07 (rotational_symmetry, far ring and small cluster)**: DEFECT to LIMIT. The literal docstring is contradicted, but `rel_tol` is in the signature, the bound line logs `abs_tol = 1.000000004553418`, and `degeneracy()` on the same set reports `coincident = True, central = True, n_at_centre = 3`. The repo already reports the limit. What is missing is that the symmetry list does not carry it.
- **A05 (angle_histogram, width 7 on [0, 180])**: stays DEFECT, on the strongest basis available: the function's own record says `hi = 180` and its 26 bins reach 182. Only for widths that do not divide `hi - lo`.
- **A19 (write_run, NaN and Infinity)**: stays DEFECT, against RFC 8259 only. Section 6: "Numeric values that cannot be represented in the grammar below (such as Infinity and NaN) are not permitted." Section 10: "The resulting text MUST strictly conform to the JSON grammar." The repo rule that matches (`docs/schema.md`: "nonfinite numbers and duplicate object keys are rejected") governs the other build. Seven committed run files in `runs/` already hold 191 bare `+Infinity` values (from `float("inf") if r is None else r` in `fold_loop.py`, `fold_shape.py`, `layered_pursuit.py`). Python reads them; JavaScript `JSON.parse` does not.

Limits of this test: 20 cases written by one author; the oracles for A14 and A17 are short and were checked against known
values; A09's trial generator targets one failure shape; A15 timing is this machine only. "DEFECT" means the probe disagrees with its basis
(section 1a says which kind), not "wrong for the repo's current use". The audited classes are a post-hoc revision by
the author of the predictions. Whether each matters for a stored run is for the maintainers to decide. Opt-in fix
proposals are in a separate PR.

## 2. Q5 across slopes (`experiments/q5_slopes.py`, `runs/q5_slopes_seed0.json.gz`)

Same pipeline as `window_extent.py` (hard window centred on the cut, gaps clustered at 1e-6, about 1,200 points).
Width w is in the internal coordinate n - s m; the canonical strip is w = 1 + |s|. Any count other than 2, and any
reading with an under-resolved tag, was re-read at 20,000 points: 257 re-reads, 0 disagreements, 0 tags.

- **P1 held (OBSERVED).** At the canonical width the count is exactly 2 for all 158 irrational slopes (8 named, 150 random in (-3, 3)).
- **P2 refuted (OBSERVED).** I predicted at most 2 for rational slopes. Slope 1 reads 1/3/2 and slope 3 reads 2/4/3 at w(1-1e-9)/w/w(1+1e-9). Slopes 0, 1e-6, 2/5, 1/2, 100 read 1 or 2. At slopes 1 and 3 the exact width puts window edges on lattice points, so the middle reading is a tie; the +-1e-9 readings are the ones to trust (1 and 2, 2 and 3).
- **P3 held where testable (OBSERVED).** Every continued-fraction convergent error with p, q <= 6 gives 2. The one exception in the table (slope phi, q=5, p=8) had n = 8, outside the scan; read directly it gives 2 (`runs/q5_slopes_followup_seed0.json`).
- **P4 refuted as written (OBSERVED).** r = 1 is phi^0 for every slope, so "not powers of phi" fails trivially. The substance: for slopes 1/phi, -1/phi and phi all 8 two-length widths are powers of phi; for sqrt2-1, sqrt3-1, 1/sqrt2, e-2 and pi-3 only r = 1 is.
- **P5 held (OBSERVED).** Every 3-length reading passes `three_gap_sum` at 1e-6.

What the scan shows beyond the predictions (EXPLORE, PROPOSED):

- Two-length widths in the scan (m <= 6, |n| <= 6): every one-sided record approximation of |s| (the m where frac(m|s|) or 1 - frac(m|s|) sets a new minimum), the width w = 1, the canonical width w = 1 + |s|, and one or two further widths (2 + |s| for slopes below 1). No record width gave 3; every other lattice-vector width tested gave 3.
- `runs/q5_slopes_followup_seed0.json`: for 5 irrational slopes outside Q(sqrt5), all 45 record widths with m <= 20 (w >= 0.02) gave 2. Moving the window off-centre changed the count in none of 55 cases (centre at -0.9, -0.4, +0.3, +0.8 of the half-width).
- The rarest gap, at two-length golden widths, is always 0.3823 = 1/phi^2 of gaps (OBSERVED in `runs/q5_special_widths_seed0.json`).
- Not tested: records beyond m = 20, windows that are not intervals, 3 or more dimensions, rational slopes other than 7, stretches below about 1,200 points.

## 3. Reproduce

```text
python3 experiments/adversarial_probes.py [--small]
python3 experiments/q5_slopes.py [--small]
python3 experiments/q5_slopes_followup.py [--small]
python3 experiments/q12_q11_checks.py
python3 experiments/q5_q9_sweeps.py q5|q9 [--small]
python3 experiments/q5_special_widths.py [--small]
python3 experiments/q9_fit_diagnostic.py
python -m unittest discover -s tests -v
```

Runs: about 4 s (adversarial), 90 s (q5_slopes, one core), 9 s (followup), 1 s (q12/q11), 60 s (q9, 100 seeds).
`--small` runs write their own files (`..._small_seed0.json`); those are not committed.
