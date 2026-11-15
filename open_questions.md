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

**Q2. 24-cell.** The same question for `cell24()`.
- Resolved when: as for Q1.
- Tools: `primitives.cell24`, `rotate_project`, `rotational_symmetry`.

**Q3. 5-cell search path.** `simplex_4d` reaches the symmetric (Petrie)
view only by construction; 2000 random rotations never came near it
(`runs/simplex_4d_seed0.json`). Find a path through rotation angles from a
random start to the Petrie plane along which `pentagon_score` never rises.
Or find where every path must rise.
- Resolved when: a logged path, scored step by step.
- Tools: `experiments/simplex_4d.py`, `primitives.rotate_many`.

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

**Q9. Half-mass dimension.** For the normalised-Gaussian sampler in
`near_orthogonality`, at which dimension does the share of angles within
90 ± 5° cross 0.5? At 400 pairs, does the crossing move with the seed?
- Resolved when: a dimension sweep over several seeds, logged.

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

## Projection and the source lattice

Q12 is void: no prompt carries that number.

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

