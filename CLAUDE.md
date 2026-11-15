# CLAUDE.md

Guidance for working in this repository. CC0. Co-authored.

**Two builds of one spec are kept here in parallel, on purpose** (see
`docs/two_builds.md`). The rules, layout and schema below are the
**Claude Code build**'s (`toolkit/`, `experiments/`, `emlib.py`,
`runs/*.json`). The **Manus build**'s own rules are in `CONTRIBUTING.md`
and `docs/schema.md` and govern `explorative_math/` and `examples/`. Do not
unify the builds, refactor one into the other, or delete duplicates between
them; compare them with `compare/` instead.

## Purpose

Explorative geometry. Map where human curiosity and AI curiosity land on the
same objects, and record **all** outputs, including the ones that look
uninteresting. A run nobody found interesting is still a run, and it is
kept. The log is the product; no single result is.

## Rules

```text
license          CC0, everything
python           stdlib only: math, itertools, random, fractions, numbers, statistics for
                 the maths; zlib for the compressibility probe; json, pathlib,
                 functools, gzip, unittest for I/O, the log hook and tests.
                 No numpy, no pip.
                 phone-buildable: every script runs in seconds on a phone.
                 A run that takes more than a few seconds gets a --small bound
                 option (fewer starts / evaluations, own run file); outputs are
                 never cut to make it fit.
words            no word coinage. Use an existing name or leave the slot unnamed.
gaps             leave gaps. A slot with nothing in it stays void (null),
                 never dropped and never filled to look complete.
readings         "not interpretable atm" is a legal reading. Use it rather
                 than force one.
shape            keep shape: record the object as computed before any summary.
bounds           any count / period / symmetry reading logs the instrument bound
                 it was read at (points, precision, probe step); unknown bound
                 = void, not omitted. See docs/resolution_leak.md.
response         in field/, "ideal" is the REFERENCE sheet, not a target. A cell's
                 difference from it is its RESPONSE to the local field, logged
                 beside the field value; it is not called error. The optimum is
                 relative to the environment.
                 The toy is CONSTRUCTED: never describe its numbers as findings
                 about a material.
sampling         a "none found" / "largest seen" reading logs how samples were drawn
                 and how many, and whether the target is measure zero there.
status tags      every claim is OBSERVED (computed here), DERIVED (follows from
                 computed values, reasoning shown) or PROPOSED (not yet run).
interest fields  human_interest is the human's to fill; the AI never writes it.
                 ai_interest is the AI's own reading, labelled as such, or void.
                 Automatic runs (toolkit.loghook, toolkit.generator) leave BOTH
                 interest slots void: nobody has looked yet.
co-authored      tag new files and commits co-authored.
```

## Method

The carrier's term for this method has no English name. It stays unnamed here.

```text
shape -> rotate (n-dim) -> project -> measure angles
      -> compose (phi-scale / second shape) -> compare nature catalog
```

Known related names, kept as pointers, not as a name for the method:
- grand tour
- projection pursuit
- Coxeter / Petrie projections
- cut-and-project
- principal angles

## Layout

```text
emlib.py               shared: run-file writer (+ .json.gz for large runs), load_run,
                       run_files, schema check, LZ76 + zlib
                       compressibility (vector helpers re-exported from toolkit)
exact/                 numbers with no floats: ZPhi (a + b*phi, a, b Fractions,
                       phi^2 = phi + 1; add, mul, inverse, compare, conjugate),
                       Cyc5 (plane points on fifth roots of unity), inflation
                       (Fibonacci chain, pentagon, Penrose subdivision), fold (reflect a
                       point set into its mirror wedge, with multiplicities), complement
                       (squared-distance field of an exact pattern and its exact symmetries)
field/                 hexsheet: CONSTRUCTED toy. Ideal hexagonal sheet (exact
                       reference), a field over it, a local relaxation rule, and
                       each cell's RESPONSE (never called error); saturating rule,
                       performance (strain, net pull, energy), inverse_readout
toolkit/
  primitives.py        vectors, rotation in any (i, j) plane, projection, angles,
                       Gram-Schmidt, polytopes (simplex, hypercube, cross-polytope,
                       24-cell, icosahedron, regular polygon), lattices (Z^n, hex),
                       reflect, coxeter_plane (from simple roots and Coxeter number)
  composers.py         rotate_project, phi_scale / phi_stack, cut_and_project
                       (hard window and/or weight function; gaussian_weight),
                       iterate, overlay
  probes.py            rotational / mirror symmetry (grid, and from pairs), smallest period, period
                       profile, zlib ratio, angle histogram, polar angles
                       three_gap_sum (L3 = L1 + L2: pass / fail / not_applicable)
                       min_stretch_tag (points < 3 / third-gap freq: under-resolved)
                       *_bound helpers: the bound a reading was taken at
                       mirror_residual (continuous; 0 = exact mirror)
                       mirror_axis_residuals, mirror_k_score (k distinct axes)
                       mirror_residual_at (residual on given axes)
                       degeneracy (smallest pair distance, collinearity, centre points, tags)
  complement.py        read a shape from the space around it: polar-grid distance field,
                       the turns / mirrors that leave it unchanged, voids, angular gaps;
                       on the sphere, voids by a seeded hill climb
  hull4.py             facets of the convex hull of points on a sphere in R^4, brute force
  generator.py         seed -> random pipeline (shape / cut / orbit) -> probes -> log
                       cut_weighted: asked for by name, never drawn, own run files
  loghook.py           log_run() and @logged: every run writes the schema,
                       interest slots void
nature_catalog.py      entries a run can match against (all CARRIED, not verified)
experiments/
  near_orthogonality.py   a  random unit vectors: angle distribution vs dimension
  shell_concentration.py  b  ball volume fraction near the surface vs dimension
  fibonacci_chain.py      c  Z^2 cut at slope 1/phi, projected to a 1D chain
  golden_angle.py         d  phyllotaxis points: golden angle vs perturbed angles
  simplex_4d.py           e  5-cell: rotate in 4D planes, project to 2D
  weighted_window.py      f  one Z^2 cut, two windows: hard strip and Gaussian weight
  window_extent.py        g  one logged window on a longer stretch; width beside m + n*slope
  three_gap_check.py      h  three_gap_sum probe over every logged cut run
  min_stretch_flag.py     i  under-resolved tag for logged gap readings; run files untouched
  orbit_exact.py          j  seed-5 orbit in exact arithmetic, then past 512 steps
  symmetry_resolution.py  k  rotational_symmetry at k_max 12 / 24, tolerance 1e-6 / 1e-9;
                             mirror_symmetry on the Q12 unions
  mirror_sweep.py         l  what frame 0 is; mirror axes on 360 tour frames per shape
  coxeter_plane.py        m  Coxeter-plane projections of four 4-polytopes, read for mirrors
  mirror_pursuit.py       n  coordinate descent on mirror_residual from random starts; path logged
  layered_pursuit.py      o  label layered_unanchored: descent on the sum of the k smallest
                             axis residuals, staged k = 1..h, and k = h at once; --small
  anchor_generate.py      p  one mirror by descent, then axes a + k*180/m by rule
                             (m = h, controls h-1, h+1); set beside layered_unanchored
                             and alt; --small
  degeneracy_check.py     q  coincident / collinear points in every final projection;
                             recount of starts reaching h with those excluded
  constrained_generate.py r  anchor_generate with the anchor HELD (Gauss-Newton pull-back,
                             anchor residual < 1e-6); controls; --small (full run ~1 min)
  pairing_reread.py       u  constrained_generate with the anchor's pairing read afresh at
                             every pull-back, beside the kept-pairing run; --small (~30 s)
  anchor_selection.py     w  anchors chosen by pairing type (points on the axis), re-anchored
                             otherwise; beside the anchor as found; --small (~30 s)
  axis_split.py           y  the 16-cell's 2-on-axis anchors, reach vs not, start by start:
                             on-axis points, pairing, generators, residual path, degeneracy
  fold_loop.py            x  label fold_search_space. anchor -> fold -> rerun: search candidate axes only inside the
                             wedge left by the mirrors found; wrong-fold and no-fold
                             controls; 5-cell, 16-cell, 24-cell (2 starts); --small
  fold_shape.py           z  the carrier's fold: reflect the projected points into the wedge
                             of the mirrors found, pursue on the folded set, confirm on the
                             shape, fold again; wrong-fold and no-fold controls; exact fold
                             of the pentagon; --small
  combine.py              aa fold shape, anchor by type, rule m = h and weighed multiplicities,
                             together and one at a time through one loop; --small
  collapse_gap.py         ab collapse a shape to one point, read it from the complement;
                             probes vs baseline, random controls; plane, sphere, exact; --small
  complement_fit.py       ac fit of a probe shape in a baseline's complement (residual gap),
                             pairs and triples, and the same under a field gradient; --small
  complement_fit_v2.py    ad fit at an equal point budget; points that respond to the field
                             (hexsheet relaxation) before fitting; sphere rotations searched
                             and joined sets set beside a catalogue; --small (full ~30 s)
  sphere_compounds.py     ae two 16-cells (searched, and built in exact fractions), two
                             24-cells in dual position, and the hulls of two carried
                             catalogue entries; --small
  dual_partner.py         af is a baseline's best-separated probe its dual? rank of the dual,
                             void cover, controls with no dual in the set; --small (full ~50 s)
  exact_inflation.py      s  Fibonacci / pentagon / Penrose inflation in exact forms:
                             two constructions per level must be EQUAL; --small
  field_response.py       t  hexagonal sheet under a gradient field + uniform control;
                             per-cell field value and response; --small
  field_adaptation.py     v  inverse readout of the field from offsets; ideal reference
                             sheet vs adapted sheet under one field; linear and
                             saturating rules; --small
runs/*.json            one file per run, schema below
runs/*.json.gz         the same, gzipped, for a run above emlib.GZ_THRESHOLD bytes;
                       full content; read either form with emlib.load_run
tests/                 unittest, one file per experiment plus emlib
DIVERGENCE_LOG.md      where human and AI interest differ (columns = schema)
open_questions.md      puzzle layer: unsolved prompts, no answer key
REFERENCES/            pointers to related studies, each with how far it was read
```

## Run schema (`runs/*.json`, `runs/*.json.gz`)

```text
params            inputs, including the seed
outputs           everything computed, before any summary
angles            measured angles (degrees) or null when the run has none
human_interest    null until the human fills it
ai_interest       the AI's reading, "not interpretable atm", or null (automatic runs)
compressibility   LZ76 phrase count and zlib ratio of a declared symbol string,
                  with its definition
nature_match      a nature_catalog id, or null
divergence_note   null until the two interest fields differ
```

`emlib.write_run` refuses a run that is missing a key, or that carries a
non-null `human_interest`.

Puzzles in `open_questions.md` carry no answers. Do not add one; an attempt
is a logged run, and a prompt stays OPEN until a run file shows the probe
reading it asks for.

## Commands

```sh
python3 -m unittest discover -s tests          # all tests
python3 experiments/<name>.py                  # writes runs/<name>_seed<k>.json (or .json.gz)
python3 experiments/layered_pursuit.py --small # phone bound -> runs/<name>_small_seed<k>; also anchor_generate,
                                               # constrained_generate, pairing_reread, anchor_selection,
                                               # fold_loop, fold_shape, combine, collapse_gap,
                                               # complement_fit, complement_fit_v2, sphere_compounds,
                                               # dual_partner,
                                               # exact_inflation,
                                               # field_response, field_adaptation
python3 -m toolkit.generator 0 1 2             # writes runs/generated_seed<k>.json
python3 -m toolkit.generator --pipeline cut_weighted 0 1 2   # runs/generated_cut_weighted_seed<k>.json
python3 -m compare                             # cross-build comparison -> compare/runs/*.jsonl
python3 -m explorative_math --help             # the Manus build's CLI
```
