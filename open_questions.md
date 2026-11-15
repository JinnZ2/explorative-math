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
