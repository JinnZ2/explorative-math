# Resolution leak

Status: a class of reading error, defined here, with three instances, all
three now run. Each row below carries its own status; two rows split,
because the bound that was tested and the bound that was suspected were
not the same bound. The name is the
operator's. Co-authored. CC0. Claude Code build only (`toolkit/`,
`experiments/`, `runs/`).

## Class

```text
object ──► instrument (bounded) ──► reading ──► log
                 │                               ▲
                 │  points looked at             │
                 │  precision / tolerance        │  the bound does not travel
                 │  probe step / probe range     │  with the reading
                 └───────────────────────────────┘
```

A count, a period or a symmetry order is read through an instrument that
has a bound: how many points it saw, at what precision it compared, over
what step or range it searched. Change the bound and the reading can
change while the object stays the same.

The reading is a resolution leak when it is logged without that bound, so
that a property of the instrument is carried forward as a property of the
object.

It is not an error in the computation. Every value involved is computed
correctly. The loss is at the join between the reading and its log.

## Instances

| # | Reading | Where | Bound tested | Status |
|---|---|---|---|---|
| 1 | count: 2 gap lengths | `runs/generated_cut_weighted_seed1.json`, level 0.9 | points: 57 against 99..1997 | CONFIRMED (points) |
| 2a | symmetry: same list for a regular 8-gon and a regular 16-gon | `runs/symmetry_resolution_seed0.json`, `known` | `k_max` 12 against 24; tolerance 1e-6 against 1e-9 | CONFIRMED (`k_max`); tolerance: no change |
| 2b | symmetry: same list for the two unions of Q12 | `runs/symmetry_resolution_seed0.json`, `unions`; `compare/runs/phi_scale_vs_phi_scale_rotate.jsonl` | the same four bounds | not a leak at the bounds tested |
| 2c | mirror: axis lists for the two unions of Q12 | `runs/symmetry_resolution_seed0.json`, `unions_mirror` | 360 against 720 axes; tolerance 1e-6 against 1e-9; grid against no grid | not a leak at the bounds tested |
| 4 | count: 7 mirror axes on a 5-point set | `probes.mirror_symmetry_from_pairs`, first `mirror_pursuit` reading (not committed) | merge precision: accepted axes within 1e-9 degrees were one axis | CONFIRMED (precision); probe changed |
| 3a | period: 184 at 256 and 512 steps | `runs/orbit_exact_seed0.json`, `exact` | precision: binary float against exact rational | not a leak (float precision); open |
| 3b | the same reading | `runs/orbit_exact_seed0.json`, `extended` | steps: 512 against 1024, 2048, 4096 | CONFIRMED (steps) |

**1. CONFIRMED (points).** OBSERVED in `runs/window_extent_seed0.json`: the
same slope and the same window read 2 gap lengths at 57 and 99 points and
3 at 199, 501, 1001 and 1997 points. The third length is 20 of 1996 gaps on
the longest stretch. The 57-point reading was logged with its point count
in the file, but the count "2" was quoted without it.

**2a. CONFIRMED (`k_max`).** OBSERVED: at `k_max` 12 the probe returns
`[2, 4, 8]` for both the 8-gon and the 16-gon, at either tolerance. At
`k_max` 24 it returns `[2, 4, 8]` and `[2, 4, 8, 16]`. Two shapes of known,
different order read the same below the larger order and separate above it.
Moving the tolerance from 1e-6 to 1e-9 changed no output in this run.

**2b. Not a leak at the bounds tested.** OBSERVED: on all 8 frames the two
unions return the same list at all four bounds. DERIVED: a set of n points
has no rotational order above n, and each union has 10 or 16 points, so
`k_max` 24 covers every order they could have. The equal lists are
therefore not made by `k_max`. Not tested: a coarser tolerance, a centre
other than the centroid, other frames or shapes. Q12 stays as written; its
own test asks for more shapes than the two run here.

**2c. Not a leak at the bounds tested.** Prediction logged before the
reading (DERIVED, the operator's): the mirror lists separate; the rotation
lists are equal by construction. OBSERVED: the rotation lists are equal on
all 8 frames at all four bounds. The mirror lists separate on 2 of 8
frames: simplex frame 0 (Build A's union no axis; Build B's 30, 90, 150
degrees) and cross-polytope frame 0 (A none; B 0, 45, 90, 135 degrees). On
the other 6 frames both unions read no mirror axis. The readings are the
same at 360 and 720 axes, at both tolerances, and with no axis grid
(`probes.mirror_symmetry_from_pairs`, candidates from every pair of
points). The empty lists on those 6 frames are therefore not made by the
axis step. So the prediction holds on 2 frames and does not on 6; in every
frame the mirror reading separates the two unions or shows no mirror in
either.

Follow-up to 2c (`runs/mirror_sweep_seed0.json`). OBSERVED: frame 0 of
Build A's tour has every rotation angle zero, and its projected points are
the first two coordinates of the unrotated vertices. Along the same tour at
360 frames per shape, a mirror shows on 1 frame of 360 for each shape,
frame 0, for the shape alone and for Build B's union; on no frame for
Build A's union. Bound: one tour (seed 0), steps of 1/360 of each plane's
sweep, tolerance 1e-6, no axis grid. An orientation with a mirror between
two frames, or off this path, is not seen.

**3a. Not a leak (float precision); open.** OBSERVED: the 512-step word in
exact rational arithmetic equals the float word at every index; the largest
float-minus-exact difference is 0. The period profile is the same,
4, 4, 27, 39, 117, 184, 184. DERIVED: the logged alpha and x0 are both
multiples of 2^-52, so every sum below 2 is exact in a 64-bit float; this
orbit had no rounding to leak. That holds for this seed's values, not for
the orbit pipeline in general. By the rule the run was made under, the
reading persists, so it is logged as open, not as a float leak.

**3b. CONFIRMED (steps).** OBSERVED, same run: carried past 512 steps in
exact arithmetic, the prefix period reads 184 at 512, 689 at 1024, 1699 at
2048 and 3719 at 4096. The stall at 184 was the end of the orbit. Also
logged: 184 x alpha = 47.001883, 1.883e-03 from an integer. The extended
profile is itself a reading bounded at 4096 steps.

Note on row 3. DERIVED: alpha and x0 are dyadic (multiples of 2^-52), and
alpha in lowest terms has denominator 2^52. The orbit x + n*alpha (mod 1)
is therefore exactly periodic with period 2^52, and the coded word's period
divides 2^52. So 184 was never the period: it was the smallest period of a
512-letter prefix, a prefix bound. Status per 3b: CONFIRMED. The true word
period is not logged; no run here reaches 2^52 steps.

The extended reading (3b), the union reading (2b) and the no-grid mirror
check in 2c were added by the AI
co-author beside the tests that were asked for (3a, 2a, the mirror grid in 2c). The period
instance was also the AI co-author's choice for the third row; the
operator's own list was not in the relay.

**4. CONFIRMED (precision); probe changed.** Found while reading the
pursuit endpoints. A descent ends on a mirror that is exact to about 1e-10,
not to 0. Each pair of points then proposes the same axis at a slightly
different angle, all pass the matching tolerance, and the probe merged only
candidates within 1e-9 degrees: one axis was counted once per pair (3, 5,
7 axes on a 5-point set). A set of 5 distinct points cannot have 7 mirror
axes, which is how it was seen. The probe now treats accepted axes within
degrees(4 x rel_tol) as one. The three earlier runs that use it
(`symmetry_resolution`, `mirror_sweep`, `coxeter_plane`) were remade and
every value is unchanged: their mirrors were exact, so the pairs agreed.
The wrong counts were never committed.

## Lesson

The bound suspected is not always the bound that matters; test all logged
bounds, not only the hypothesized one.

```text
row   bound suspected      bound that moved the reading
2     tolerance, k_max     k_max (known shapes); none for the unions
3     float precision      steps
```

## A second entry type: sampling leak

Status: defined here with three instances, all OBSERVED in logged runs.
The name is the operator's.

```text
resolution leak    the bound is too coarse       raise the bound -> the reading moves
sampling leak      the target has measure zero   more samples of the same kind -> no change
                   in the space being sampled
```

A sampling leak is a "none found" or a "largest seen" reading, taken by
sampling orientations (at random, or along a path chosen without regard to
the target), when the thing looked for occupies a set of measure zero in
the space sampled. The sample cannot land on it. The count of samples is
logged, and it is the wrong bound: no count is large enough. Logged as a
property of the object, "none in N" reads as "there is none".

It is not cured the way a resolution leak is. What reaches the target is a
different instrument: a construction (compute the special orientation), or
a continuous reading (how far from the target), which shows approach where
a hit / miss probe shows nothing.

| # | Reading | Where | Sampled | What a different instrument showed |
|---|---|---|---|---|
| S1 | no pentagon view of the 5-cell in 2000 random rotations (lowest score 23.45, median 152.8) | `runs/simplex_4d_seed0.json` | 2000 random rotations | the Petrie projection, by construction, scores 2.2e-14 (same run) |
| S2 | a mirror on 1 frame of 360 per shape, frame 0 only. **DEGENERATE:** frame 0 has coincident points for both shapes; no frame with distinct points has a mirror. | `runs/mirror_sweep_seed0.json`, `axis_count_per_frame` | one seeded tour, 360 frames | `mirror_residual` on the same frames dips off frame 0: 5-cell 5.5e-03 at frame 137; 16-cell 9.0e-04 at frame 178 and 5.9e-04 at frame 1 (median over frames about 1e-01) |
| S3 | largest mirror count on the tour 3 (5-cell) and 4 (16-cell). **DEGENERATE:** both are read on frame 0, where points coincide (4 and 5 distinct points). Excluding degenerate frames the largest count on the tour is 0 for both. | `runs/mirror_sweep_seed0.json` (`max_axis_count_degenerate`, `max_axis_count_excluding_degenerate`) | the same tour | Coxeter-plane projections, by construction, have 5 and 8, on distinct points (`runs/coxeter_plane_seed0.json`) |

S2 and S3 do not say the tour passes near a mirror orientation at the dips;
they say the hit / miss reading could not have shown it either way.

A continuous reading is not a cure by itself. OBSERVED in
`runs/mirror_pursuit_seed0.json`: descending `mirror_residual` from 8
random starts per shape reaches a mirror on every start (one axis on the
5-cell, two on the 16-cell; one 16-cell start stops at 4.0e-06 on the
evaluation limit) and reaches the Coxeter-plane count on none of the 16.
DERIVED: the residual is 0 at any projection with one mirror, so the index
has no slope toward more axes. The many-mirror orientations are measure
zero inside the one-mirror set as well; the descent trades one sampling
leak for the next one down.

Follow-up (`runs/layered_pursuit_seed0.json.gz`): an index that asks for
more axes, the sum of the k smallest distinct-axis residuals. OBSERVED, 8
starts per shape, counts at tolerance 1e-6:

```text
                         5-cell (Coxeter 5)      16-cell (Coxeter 8)
layered_unanchored       0 of 8 reach 5          0 of 8 reach 8
  largest count seen     4 (one start)           4 (three starts, lost again later)
alt, k = h at once       0 of 8 reach 5          1 of 8 reaches 8
```

This run is labelled layered_unanchored: the intended layered pursuit
holds the found mirrors fixed and descends only the next axis, and this
one does not. Its stages do not hold what an earlier stage found: the count
falls from 1 to 0 at stage 2 on every 5-cell start, because the sum of two
residuals is lowered by giving up the exact first mirror. Nearly every
descent stops on its step limit, not its evaluation limit, at a score well
above 0: a local minimum of the index for these moves. Bound: 8 starts,
one seed, coordinate moves, axes closer than 5 degrees counted as one.

Beside it, the anchored test (`runs/anchor_generate_seed0.json.gz`): one
mirror by descent, then the axes a + k x 180/m generated by rule and
descended as a fixed set. OBSERVED, same starts, tolerance 1e-6, 700
evaluations (4000 in brackets):

```text
                         5-cell (h = 5)      16-cell (h = 8)
rule m = h               2 of 8 (2)          1 of 8 (2)
control m = h - 1        0 of 8 (0)          0 of 8 (0)
control m = h + 1        0 of 8 (0)          0 of 8 (0)
```

Bound on the word "anchored": in that run the anchor axis is a term of
the score and the rotation is free to give it up. Under rule m = h the
anchor ends exact (axis residual below 1e-6) on 6 of 8 5-cell starts and
3 of 8 16-cell starts.

The constrained run (`runs/constrained_generate_seed0.json.gz`): the same
anchors and rules, with every accepted position pulled back onto an exact
anchor mirror (Gauss-Newton on the reflected-point mismatch) and refused
if the anchor residual is not below 1e-6. OBSERVED, tolerance 1e-6:

```text
                         5-cell (h = 5)              16-cell (h = 8)
                         400 trials   2500 trials    400 trials   2500 trials
rule m = h                   0            1              1            5
control m = h - 1            0            0              0            0
control m = h + 1            0            0              0            0
anchor held at the end       8 of 8 on every rule, both limits, both shapes
```

At 400 trials nearly every descent stopped on the trial limit with its
step still large; the 2500-trial rerun is logged beside it as a bound
check, and the limit was the bound that moved the reading. At tolerance
1e-4 the 5-cell rule m = h reads 5 axes on 3 of 8 starts at 2500 trials.
Bound: 8 starts, one seed; the pairing of points across the anchor is read
at the anchor and kept; a trial here costs a pull-back, so trials and the
other runs' evaluations are not like for like.

Pairing type (`runs/pairing_reread_seed0.json.gz`). A mirror line fixes the
points on it and swaps the rest in pairs, so a pairing has a type: the
number of points on the axis. PROPOSED by the operator: the 5-cell is
weaker because the anchor's pairing is kept, and re-reading it would help.
OBSERVED, kept pairing, rule m = h, 2500 trials, by the anchor's type:

```text
5-cell    1 point on the axis    1 of 4 starts reach 5 (3 of 4 at tolerance 1e-4)
          3 points on the axis   0 of 4
16-cell   6 points on the axis   0 of 2 starts reach 8
          2 points on the axis   4 of 4
          0 points on the axis   1 of 2
```

OBSERVED, pairing re-read at every pull-back, same starts and limits:

```text
starts reaching h, rule m = h      kept pairing     re-read pairing
5-cell    400 / 2500 trials          0 / 1            0 / 2
16-cell   400 / 2500 trials          1 / 5            1 / 2
controls m = h-1, h+1                0 everywhere     0 everywhere
```

Under re-reading, one 5-cell start whose anchor had 3 points on the axis
ends with 1 on the axis and 5 mirrors. On the 16-cell, three starts that
reached 8 with the kept pairing end on 4 axes when it is re-read; 6 of
the 8 re-read 16-cell endpoints have coincident points, against 2 of 8
with the kept pairing. No 5-cell endpoint is degenerate either way. The
anchor is held (residual below 1e-6) on every start either way. Bound: 8
starts, one seed; the pairing is the nearest-point one at each trial
position; paths are not kept for the re-read descents.

Degenerate projections (`runs/degeneracy_check_seed0.json`). A projection
where two points land on each other has mirrors the polytope's symmetry
did not give it. OBSERVED over 334 final projections: 63 have coincident
points at 1e-6 and 33 more at 1e-3; none is collinear. Every start counted
above as reaching h is non-degenerate, so no count of starts reaching h
changes. What does change: frame 0 of the tour is degenerate for both
shapes (the 5-cell shows 4 distinct points, the 16-cell 5), so the "largest
mirror count on the tour, 3 and 4" of row S3 was read on coincident
points. DERIVED: collinearity cannot occur for these shapes. The vertex
sets have a scatter matrix proportional to the identity in R^4, so every
orthogonal projection to a plane has equal widths in all directions; the
collinearity reading is 1 on every row.

Rule for sampling readings: a "none found" or "largest seen" reading logs
how the samples were drawn and how many, says whether the target is
expected to be measure zero in that space, and where it is, is logged
beside a constructed or continuous reading, or marked as not having one.

## Rule

```text
any count / period / symmetry reading logs the instrument bound it was read at:

  points       how many points, gaps or steps the reading saw
  precision    the tolerance two values were compared at
  probe step   the step or range the probe searched (k_max, prefix lengths,
               sweep step)
```

A bound that is not known is logged void (null), not left out. A reading
quoted in prose (`ai_interest`, `open_questions.md`, these docs) carries
its bound in the same sentence.

Where the existing probes stand against the rule:

| Reading | Bound now logged beside it | Where |
|---|---|---|
| gap-length count | `n_points`, `gap_cluster_sep` | `generator` cut and cut_weighted; experiments c, f, g |
| `period_profile`, `smallest_period` | `period_bound`: word length `n`, `prefix_lengths` | `generator` cut, cut_weighted, orbit; experiment c |
| orbit word | `arithmetic` | `generator` orbit |
| `rotational_symmetry` | `rotational_symmetry_bound`: `n_points`, `k_tested`, `rel_tol`, `abs_tol`, `about` | `generator` shape; experiment k |
| `mirror_symmetry` | `mirror_symmetry_bound`: `n_points`, `n_axes`, `axis_step_deg`, `rel_tol`, `abs_tol`, `about` | `generator` shape |
| projected points fed to the symmetry probes | `point_rounding_decimals` | `generator` shape |

The bounds are new fields. The probes' own return values are unchanged, and
the regenerated run files carry every earlier value unchanged. A probe
called directly still returns only its reading: the caller logs the bound
with `probes.rotational_symmetry_bound`, `mirror_symmetry_bound` or
`period_bound`.

compare/ now logs them too: each record that reads `rotational_symmetry`
or `mirror_symmetry` carries `outputs.bounds`, beside `fields` and `tally`.
The bounds are not compared fields and do not enter a tally; every earlier
value in `compare/runs/*.jsonl` is unchanged.

Still without a logged bound: the Manus build, which this rule does not
govern.

## Flag for gap readings

`probes.min_stretch_tag` and `experiments/min_stretch_flag.py` apply one
bound check to every logged gap reading: a reading taken on fewer than
about 3 / f points, where f is the frequency of the third (rarest) gap
length, is tagged `under-resolved`. The factor 3 is a declared choice, not
a derived one. Rows are written to `runs/min_stretch_flag_seed0.json`;
the runs they describe are not changed. Where f is not known, the row
says so and carries no tag.
