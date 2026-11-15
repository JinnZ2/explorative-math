# Resolution leak

Status: a class of reading error, defined here, with three instances. One is
CONFIRMED by a logged run; two are suspected and not run. The name is the
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

| # | Reading | Where | Bound it was read at | Status |
|---|---|---|---|---|
| 1 | count: 2 gap lengths | `runs/generated_cut_weighted_seed1.json`, level 0.9 | 57 points | CONFIRMED |
| 2 | symmetry: same `rotational_symmetry` list for two different unions | `compare/runs/phi_scale_vs_phi_scale_rotate.jsonl`; `open_questions.md` Q12 | `k_max` 12, relative tolerance 1e-6, about the centroid | suspected |
| 3 | period: prefix period 184 at both 256 and 512 steps | `runs/generated_seed5.json`; `open_questions.md` Q8 | 512 steps; prefix lengths 8..512 in doublings | suspected |

**1. CONFIRMED.** OBSERVED in `runs/window_extent_seed0.json`: the same
slope and the same window read 2 gap lengths at 57 and 99 points and 3 at
199, 501, 1001 and 1997 points. The third length is 20 of 1996 gaps on the
longest stretch. The 57-point reading was logged with its point count in
the file, but the count "2" was quoted without it.

**2. Suspected.** Not run. DERIVED from the code (Q12): the probe tests
only k in 2..12 at a fixed tolerance, so orders above 12, or differences
below the tolerance, cannot show. Whether the equal lists are a real
equality or the probe's bound is the open reading of Q12; its test is the
one that would confirm or clear this instance.

**3. Suspected.** Not run. OBSERVED: the period profile of that run is
4, 4, 27, 39, 117, 184, 184 at prefix lengths 8..512. The profile stops at
512 because the orbit has 512 steps. Whether 184 holds past 512 is not
logged anywhere. The sweep Q8 asks for, taken to longer orbits, would
confirm or clear it.

Instance 3 was chosen by the AI co-author as the period case. The
operator's list of instances was not in the relay; if a different third
instance was meant, this row is to be replaced, not kept beside it.

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

| Probe | Bound in its output | Bound not in its output |
|---|---|---|
| gap-length count (cut pipelines, experiments c, f, g) | `n_gaps` or `n_points`; `gap_cluster_sep` in f, g | cluster separation 1e-6 in `generator` and experiment c |
| `period_profile` | prefix lengths (its keys) | none |
| `smallest_period` | none | the length of the word it was given |
| `rotational_symmetry`, `mirror_symmetry` | none | `k_max`, tolerance, step |

The right-hand column is a list of gaps, not a list of repairs made. No
existing run file was edited.

## Flag for gap readings

`probes.min_stretch_tag` and `experiments/min_stretch_flag.py` apply one
bound check to every logged gap reading: a reading taken on fewer than
about 3 / f points, where f is the frequency of the third (rarest) gap
length, is tagged `under-resolved`. The factor 3 is a declared choice, not
a derived one. Rows are written to `runs/min_stretch_flag_seed0.json`;
the runs they describe are not changed. Where f is not known, the row
says so and carries no tag.
