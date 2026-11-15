> **Claude Code build** README (`toolkit/`, `experiments/`, `emlib.py`, `runs/*.json`). The repository also holds the Manus build, whose README is [README.md](README.md); see [docs/two_builds.md](docs/two_builds.md).

# explorative-math

Explorative geometry: map where human and AI curiosity land on the same
objects, and keep **every** output, the uninteresting ones included. CC0.
Co-authored. Python standard library only; every script runs in well under a
second.

Read `CLAUDE.md` for the rules and the run schema. `human_interest` is void in
every run until a human fills it; `DIVERGENCE_LOG.md` is empty for the same
reason.

Status tags:
- **OBSERVED**: computed in this repo, on the run file named.
- **DERIVED**: follows from computed values, with the reasoning shown.
- **PROPOSED**: not run yet.

```sh
python3 -m unittest discover -s tests     # 46 tests
python3 experiments/<name>.py [seed]      # writes runs/<name>_seed<seed>.json
python3 -m toolkit.generator 0 1 2        # writes runs/generated_seed<k>.json
```

## Experiments

### a. near_orthogonality — `runs/near_orthogonality_seed0.json`

400 pairs of random unit vectors per dimension, d = 2 to 200.
- **OBSERVED:** the share of angles within 90 ± 5° rises from 0.058 (d=2) to 0.772 (d=200).
- **OBSERVED:** sd(cos θ)·√d stays between 0.966 and 1.043.
- **DERIVED:** consistent with sd(cos θ) ≈ 1/√d.
- **Gap:** none recorded.

### b. shell_concentration — `runs/shell_concentration_seed0.json`

The share of a d-ball's volume within ε of its surface.
- **OBSERVED:** exact values come from 1 − (1 − ε)^d. A Monte Carlo cube-rejection check covers d ≤ 6, with max |z| = 3.49 over 12 rows at 4000 samples.
  - A separate 400,000-sample check at d=3, ε=0.1 gave z = −0.28, so the spread is sampling noise, not bias.
- **OBSERVED:** the shell share passes 0.5 at d=100 (ε=0.01), d=20 (ε=0.05) and d=10 (ε=0.1).
- **Void:** the Monte Carlo column is empty above d=6, because cube rejection stops accepting points there. The ball/cube volume ratio is recorded as log10 instead.

### c. fibonacci_chain — `runs/fibonacci_chain_seed0.json`

Z² cut at slope 1/φ and projected onto the line, giving a 350-gap chain.
- **OBSERVED:** there are exactly two gap lengths, with L/S = φ to 12 digits.
- **OBSERVED:** the counts are L=216 and S=134, a ratio of 1.6119.
- **OBSERVED:** the chain is a factor of the substitution word L→LS, S→L.
- **OBSERVED:** the smallest period of each prefix (length 10, 20, 40, 80, 160, 320) is 5, 13, 21, 34, 89, 233. The rational-slope control (1/2) stays at 3 for every prefix.
- **DERIVED:** a finite word always has some period, so aperiodicity can only show as the smallest period growing with length.
- **Not interpretable atm:** why the prefix periods land on Fibonacci numbers. It is recorded, not explained.

### d. golden_angle — `runs/golden_angle_seed0.json`

Vogel's model with 400 points, r = √k and θ = kα.
- **OBSERVED:** ranked by nearest-neighbour evenness (coefficient of variation, lower is more even), the golden angle is first (0.018). The two rational turns are last (2/5 turn: 0.407; 3/8 turn: 0.399), with minimum nearest-neighbour distances of about 0.13–0.20 against 1.67 for golden.
- **OBSERVED:** golden+0.1° ranks below golden+0.5°, so the ranking is not symmetric about the golden angle.
- **Not interpretable atm:** whether that asymmetry, and the size of the ±0.1° effect, is a finite-N artefact.
- **PROPOSED:** sweep N to settle it.

### e. simplex_4d — `runs/simplex_4d_seed0.json`

The regular 5-cell in R⁴, rotated in all six coordinate planes and projected onto two coordinates.
- **OBSERVED:** the Petrie projection, computed directly from the eigenplane of the cyclic vertex permutation, is a regular pentagon (score 2.2e−14).
- **OBSERVED:** over 2000 random rotations the best score was 23.45 and the median 152.8. No trial scored below the logging threshold of 5, so the logged angle sets hold only the Petrie projection.
- **PROPOSED:** search rotations toward the eigenplane, rather than sampling uniformly, and log the path.

## Toolkit layer (`toolkit/`)

| Part | Contents |
|---|---|
| primitives | vectors; rotation in any (i, j) plane in any dimension; projection onto coordinates or an orthonormal basis; angles; Gram-Schmidt; polytope vertex generators (simplex, hypercube, cross-polytope, 24-cell, icosahedron, regular polygon); lattices (Z^n, hexagonal) |
| composers | `rotate_project`, `phi_scale` / `phi_stack`, `cut_and_project` (general lattice, bases; hard window and/or a weighting function such as `gaussian_weight`), `iterate`, `overlay` (points tagged by source) |
| probes | rotational symmetry, mirror symmetry, smallest period and period profile, zlib ratio, angle histogram |
| generator | a seed draws one of three pipelines and its parameters, runs it, applies the probes, and logs it. The pipelines are: shape → φ-stack → random rotations → 2D; Z² cut at a random slope; circle-rotation orbit. A fourth, `cut_weighted` (Z² cut with a Gaussian weight), is never drawn: it runs only when asked for (`--pipeline cut_weighted`) and writes `runs/generated_cut_weighted_seed<k>.json`, so the ten drawn runs stay as they are. The drawn spec is stored in `params`, so any run can be remade. |
| loghook | `log_run()` and `@logged` write the schema for every run, with both interest slots and the divergence note left empty |

Ten generated runs are committed (`runs/generated_seed0..9.json`, seeds 0–9). None is filtered: they include shapes whose only symmetry is the trivial 2-fold one, and a 4-point simplex projection with no symmetry at all.

- **OBSERVED:** the zlib ratio exceeds 1 on short strings (zlib header overhead). It is recorded as computed; Q10 notes the limit.

## Puzzle layer

`open_questions.md` holds OPEN prompts with no answer key.
- Five are find-the-angle challenges (tesseract, 24-cell, 5-cell search path, φ-overlay mirror change, three-gap slope).
- Five come from something a recorded run showed and did not explain.
- Q11 and Q12 (PROPOSED) are across the two builds. Q13 (PROPOSED) asks which measured properties survive projection.

## Method

The method is shape → rotate → project → measure angles → compose (φ-scale or a second shape) → compare with the nature catalog.
- It stays unnamed (see `CLAUDE.md`).
- Related names: grand tour, projection pursuit, Coxeter/Petrie projections, cut-and-project, principal angles.
- "Compose" is not run yet.
- `nature_catalog.py` entries are CARRIED from general knowledge, not verified here.
