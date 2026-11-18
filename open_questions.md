# Open questions

Unsolved prompts. **There is no answer key**, here or anywhere in the repo.

Each prompt names the probe that would say it is resolved and the toolkit
pieces that reach it. An attempt is a run: log it through `toolkit.loghook`,
whatever it shows, including a result that resolves nothing. Status for every
prompt is **OPEN** until a run file shows the probe reading it asks for.

Several prompts come from something an existing run recorded and did not
explain. Those name the run.

Co-authored. CC0.

---

## Find-the-angle challenges

**Q1. Tesseract.** Find rotation planes and angles in R⁴ where the 2D
projection of `hypercube(4)` shows the highest k that `rotational_symmetry`
reports.
- Resolved when: a logged run gives the angles and the k, and no other
  logged run beats that k.
- Tools: `primitives.hypercube`, `composers.rotate_project`,
  `probes.rotational_symmetry`.
- Recorded beside it: `runs/coxeter_plane_seed0.json` projects the
  tesseract, 16-cell, 24-cell and 5-cell onto their Coxeter planes and logs
  mirror axes and `rotational_symmetry` at `k_max` 24. It gives a projection
  basis, not rotation planes and angles, and searches no other orientation.

**Q2. 24-cell.** The same question for `cell24()`.
- Resolved when: as for Q1.
- Tools: `primitives.cell24`, `rotate_project`, `rotational_symmetry`.
- Recorded beside it: `runs/coxeter_plane_seed0.json`, as under Q1.

**Q3. 5-cell search path.** `simplex_4d` reaches the symmetric (Petrie)
view only by construction; 2000 random rotations never came near it
(`runs/simplex_4d_seed0.json`). Find a path through rotation angles from a
random start to the Petrie plane along which `pentagon_score` never rises.
Or find where every path must rise.
- Resolved when: a logged path, scored step by step.
- Tools: `experiments/simplex_4d.py`, `primitives.rotate_many`.
- Recorded beside it: `runs/mirror_sweep_seed0.json` reads mirror axes on
  360 frames of one seeded tour of the 5-cell and the 16-cell. It is a
  mirror reading along a fixed path, not a `pentagon_score` path.
  `runs/mirror_pursuit_seed0.json` logs descent paths from random starts,
  scored step by step on `mirror_residual`, not on `pentagon_score`.

**Q4. Mirror change under a φ overlay.** Overlay `icosahedron()` with its
φ-scaled copy (`composers.phi_stack`, ks (0, 1)). Find the rotation angle
at which `mirror_symmetry` changes its count of axes.
- Resolved when: two logged runs either side of the angle, with different
  axis counts.

**Q5. Three gaps.** Find a slope s for the Z² cut where
`distinct_gap_lengths` reports three lengths. Or find what about the window
keeps it at two.
- Resolved when: a logged cut run with three lengths, or a run plus a
  stated reason.
- Tools: `composers.cut_and_project`, `strip_window_1d`,
  `generator` "cut" pipeline.
- Annotation, 2026-11-15. Status: DERIVED, not run. Q5 stays OPEN.
  - OBSERVED (`runs/weighted_window_seed0.json`, slope 1/phi, windows
    centred on the cut): gap-length counts 3, 2, 3, 3, 3 at weight levels
    0.9, e^-1/2, 0.5, 0.1, 0.01. The level e^-1/2 is half-width sigma,
    which is the full width of the hard strip.
  - DERIVED: the table fits the three-distance theorem (Steinhaus): at
    most three gap lengths for any window width, two only at particular
    widths. The canonical strip width is one of them.
  - PROPOSED test: sweep the window width finely at a fixed slope and
    record the gap-length count against width. Expectation carried into
    the test: two only at the canonical width(s), three elsewhere. Which
    widths other than the canonical one give two is not stated here.
  - Not decided: whether the 3-length levels in that run count as "a
    logged cut run with three lengths". They come from centred windows of
    other widths, not from `strip_window_1d`.
  - Tools for the sweep: `composers.cut_and_project(weight=...)`,
    `gaussian_weight`; `generator` `cut_weighted` pipeline.
- Annotation, 2026-11-18 (Manus; predictions hashed first, `docs/prereg/`; details in `docs/adversarial_probes_2026-11-18.md`). Status: OBSERVED. Q5 stays OPEN.
  - OBSERVED (`runs/q5_width_sweep_seed0.json`, slope 1/phi, centred hard window, about 1,200 points, r = width / canonical width, 0.05 to 2.00 step 0.01): at most 3 lengths everywhere; 2 only at r = 1.0 on that grid.
  - OBSERVED (`runs/q5_special_widths_seed0.json`, EXPLORE): the 0.01 grid missed other two-length widths. At r = phi^k for k = -5..1 (0.0902, 0.1459, 0.2361, 0.3820, 0.6180, 1.0, 1.6180) the count is 2, and 20 other lattice-vector widths give 3. This does not match the 2026-11-15 expectation "two only at the canonical width(s)". The rarest gap is 0.3823 = 1/phi^2 of gaps at each two-length width.
  - OBSERVED (`runs/q5_slopes_seed0.json.gz`): at the canonical width the count is exactly 2 for 158 of 158 irrational slopes (8 named, 150 random in (-3, 3)); every 3-length reading passes `three_gap_sum`. Rational slopes are not limited to 2 (slope 3 reads 2, 4, 3 at w(1-1e-9), w, w(1+1e-9); slope 1 reads 1, 3, 2).
  - OBSERVED (8 slopes, m <= 6; 5 slopes to m <= 20 in `runs/q5_slopes_followup_seed0.json`, EXPLORE): every one-sided record approximation width |m|s| - n| gave 2; so did w = 1 and w = 1 + |s|; other lattice-vector widths tested gave 3. Moving the window off-centre changed the count in none of 55 cases.
  - So "three lengths" is answered (generic width, any slope tested). "What keeps it at two" is partly answered by the record-width list above; why those widths and not others is PROPOSED, not shown. Not tested: records beyond m = 20, non-interval windows, higher dimension.

## Prompts from recorded runs

**Q6. Fibonacci prefix periods.** In `runs/fibonacci_chain_seed0.json`,
the smallest period of each prefix lands on Fibonacci numbers. At another
irrational slope, which integer sequence do the prefix periods land on?
- Resolved when: logged runs at several slopes, each with its
  `period_profile` listed beside a sequence it matches, or marked "no match".

**Q7. Golden-angle asymmetry.** In `runs/golden_angle_seed0.json`,
golden+0.1° ranks below golden+0.5°; the ranking is not symmetric about the
golden angle. Does the asymmetry persist as the number of points grows, or
is it a finite-N effect?
- Resolved when: an N sweep logged, with `nn_cv` for golden ±0.1° and ±0.5°
  at each N.

**Q8. Orbit stall.** In `runs/generated_seed5.json` (alpha ≈ 0.2554), the
prefix period profile is 184 at both 256 and 512 steps, while seed 6 keeps
growing. Find which alphas stall at 512 steps and at what period.
- Resolved when: an alpha sweep logged, with the stall points listed.
- Tools: `composers.iterate`, `probes.period_profile`.
- Recorded beside it: `runs/orbit_exact_seed0.json` remakes this one orbit
  in exact arithmetic and carries it to 4096 steps. It is one alpha, not
  the sweep. See `docs/resolution_leak.md` rows 3a, 3b.

**Q9. Half-mass dimension.** For the normalised-Gaussian sampler in
`near_orthogonality`, at which dimension does the share of angles within
90 ± 5° cross 0.5? At 400 pairs, does the crossing move with the seed?
- Resolved when: a dimension sweep over several seeds, logged.
- Annotation, 2026-11-18 (Manus; predictions hashed first, `docs/prereg/`). Status: OBSERVED.
  - OBSERVED (`runs/q9_half_mass_seed0.json`): 100 seeds x dims 30..100 step 5 x 400 pairs; per-seed straight-line fit of the share against dimension, solved for 0.5. Pooled mean crossing 64.19, between-seed SD 1.82, range 59.5 to 68.9. So the crossing does move with the seed at 400 pairs, by about 2 dimensions (SD).
  - OBSERVED (`runs/q9_fit_diagnostic_seed0.json`, EXPLORE): the exact share (numerical integral) crosses 0.5 at 61.16. The same straight-line fit applied to the exact curve over 30..100 gives 64.02, so about 3 of the offset is the fit window, not sampling. Fitted over 45..80 the exact curve gives 61.96 and the empirical pooled mean 62.14 (SD 2.62).
  - Bound: the exact value is a Simpson integral, not a closed form.

**Q10. Two compressibility probes disagree.** Every run records both an
LZ76 phrase count and a zlib ratio of the same symbol string. Find two runs
that the two measures rank in opposite order.
- Resolved when: the pair, with both values, written beside the run files.
- Known limit: zlib's header makes its ratio exceed 1 on short strings, so
  compare runs of similar length.

## Across the two builds

**Q11. Where do the builds disagree, and which choice causes it?**
The repository holds two independent builds of one spec: Manus
(`explorative_math/`) and Claude Code (`toolkit/`). `compare/` logs every
field where they meet as AGREE / DIFFER / ONLY_IN_A / ONLY_IN_B
(`compare/runs/*.jsonl`, summary in `docs/two_builds.md`). For each DIFFER,
which design choice produces it? Is it the seeding rule, the normalisation,
the meaning of a composition, the measurement origin, or float rounding?
And which disagreements survive once that choice is aligned?
- Resolved when: a logged comparison run, for each DIFFER group, with the
  choice aligned and the field re-checked.
- Do not resolve it by making the builds identical.
- Annotation, 2026-11-18 (Manus; predictions hashed first, `docs/prereg/`). Status: OBSERVED for the groups below; the schema group was not re-checked.
  - OBSERVED (`runs/q11_aligned_seed0.json`): each DIFFER group re-read after aligning the named choice. Constants: both golden-angle floats within 1 ulp of a 60-digit value (A -0.19 ulp, B +0.81 ulp); aligned difference 0. Shell formula: A's `-expm1(n*log1p(-t))` is closer to the exact value in 20 of 21 differing cells; B's is closer at d = 20, t = 0.1 (errors 4.9e-17 against 6.2e-17). Hypercube and cross-polytope: identical after a common circumradius and sort. Simplex: same ambient dimension and same pairwise-dot multiset, so the same shape in a different orientation. phi composition: scaling by phi and rotating by +golden angle reproduces A's union to 3.1e-16 on 8 of 8 frames; -golden misses by 0.14 to 0.88. Orthogonality statistics: with A's stream replaced by `Random(0)` in memory, 400 of 400 pairs are identical at d = 2, 20, 50, 100, 200 and the mean-angle difference is 0.
  - So each DIFFER group is explained by one choice: stream, omitted rotation, vertex order and orientation, or float rounding.

**Q12. Same symmetry answer, different unions: real equality or probe
resolution?** Status: **PROPOSED** (not run).

In `compare/runs/phi_scale_vs_phi_scale_rotate.jsonl`, Build B's
`rotational_symmetry` probe returns the same list for both builds' unions on
all 8 frames. The unions are not the same: on cross-polytope frames 1–3,
Build A's union (`phi_scale_rotate`) has 16 distinct directions and Build B's
(`phi_scale`) has 8.

OBSERVED, from the jsonl:
- Cross-polytope frame 0: both `[2, 4]`.
- Cross-polytope frames 1–3: both `[2]`.
- Simplex frame 0: both `[3]`.
- Simplex frames 1–3: both `[]`.

Two readings are both open:
- *Real equality.* Rotating a k-fold set by a fixed angle, then adding the
  copy, keeps the same k, so the order really is the same.
- *Probe resolution.* The probe cannot see the difference. DERIVED from the
  code: `rotational_symmetry` tests only k in 2..`k_max` (default 12), at a
  relative tolerance (default 1e-6), about the centroid. So a 16-fold set
  could only show up as its divisors ≤ 12, and two sets whose orders differ
  above 12, or by less than the tolerance, read the same.
  - "16 vs 8" is a count of distinct directions, not a symmetry order.

- Test: run the probe on shapes whose symmetry orders are known and
  different (for example regular polygons with k = 5, 8, 13, 16, 24, and a
  k-gon plus its φ-scaled rotated copy). Report, for each pair, whether the
  probe separates them, and at which `k_max` and tolerance.
- Resolved when: that run is logged, with each shape's known order beside
  the probe's output.
- Tools: `toolkit.primitives.regular_polygon`, `toolkit.composers`
  (`phi_stack`; a scale-and-rotate copy), `toolkit.probes.rotational_symmetry`.
- Recorded beside it: `runs/symmetry_resolution_seed0.json` reads the probe
  on a regular 8-gon and 16-gon, and on the two unions, at `k_max` 12 and
  24 and tolerance 1e-6 and 1e-9. Two of the shapes this test lists, not
  all of them. The same run reads `mirror_symmetry` on the two unions;
  `runs/mirror_sweep_seed0.json` reads it on 360 frames per shape.
  See `docs/resolution_leak.md` rows 2a, 2b, 2c.
- Annotation, 2026-11-18 (Manus; predictions hashed first, `docs/prereg/`). Status line above is unchanged; note that `runs/symmetry_resolution_seed0.json` exists, so "not run" may be stale.
  - OBSERVED (`runs/q12_extended_seed0.json`): regular k-gons for k = 5, 8, 13, 16, 24 at k_max 12, 24, 32 read exactly the divisors of k up to k_max (13: [] at 12, [13] at 24 and 32). A k-gon united with its phi-scaled copy rotated by 7 deg or by 180/k deg reads the same list as the bare k-gon in all 30 union rows.
  - DERIVED: for a centred polygon, rotation by 360/k commutes with the rotate-and-scale that makes the copy, so the union is invariant under it; the equality is real, not probe resolution. This covers centred polygons only, not the cross-polytope frames in `compare/runs/`.

## Projection and the source lattice

**Q13. What survives the shadow?** Which measured properties survive
projection unchanged, and which need the source lattice to distinguish?
Take two different source lattices, with cuts and windows chosen so that
both give the same shadow (the same projected point set, or the same gap
word). Apply every probe to both. Sort the probes into those that read the
same on both shadows and those that differ only once the source is given.
- Status: PROPOSED. No run has been made. The pair of lattices is not
  chosen; choosing one is part of the prompt.
- Resolved when: a logged run for each source, the two shadows shown equal
  by a stated test, and each probe listed as "same on the shadow" or
  "needs the source", with its reading on both.
- Tools: `composers.cut_and_project` (hard window, or `weight=` with
  `gaussian_weight`), `primitives.lattice_zn`, `lattice_hex`, `probes`.
- Recorded beside it: `runs/weighted_window_seed0.json` logs one lattice
  under two windows. That is one source, so it is not an attempt at Q13.
- Pointers: `REFERENCES/related_studies.md` R1.

**Q14. Third-gap frequency against distance from a 2-gap width.**
Status: **PROPOSED** (not run).

At a window width that gives two gap lengths, the third length has
frequency zero. Near such a width it is rare
(`runs/window_extent_seed0.json`: 20 of 1996 gaps, where width x norm sits
8.598e-03 from the nearest m + n*slope with |m|, |n| <= 10). Is the
frequency of the third gap proportional to the residual of width x norm
from the nearest m + n*slope?
- Test: fix a slope. Sweep window widths near a 2-gap value. At each
  width log the third-gap frequency, the residual, the nearest (m, n), and
  the bounds the reading was taken at: number of points, gap cluster
  separation, the range of m and n searched, the sweep step.
- Resolved when: that sweep is logged with frequency beside residual at
  every width, and a stated check of proportionality with its reading
  (holds, fails, or not interpretable atm).
- Tools: `experiments/window_extent.py` (`band_lattice`, `nearest_combos`),
  `composers.cut_and_project`, `probes.three_gap_sum`,
  `probes.min_stretch_tag`.
- Known limit: a stretch too short for the frequency it is measuring
  reads zero. See `docs/resolution_leak.md`.

