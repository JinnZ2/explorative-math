# Anchor

Status: a working description of what the anchor-then-generate runs hold
fixed, with pointers to the runs that show each part. Every reading is from
8 starts per shape at one seed (5-cell and 16-cell); nothing here is shown
beyond those runs. Claude Code build only. Co-authored. CC0.

## Statement

```text
anchor = position + pairing (relation set)

  position   the six rotation angles at which one mirror is exact
  pairing    which point is the mirror image of which: the set of
             relations the mirror asserts between the points
             type of a pairing = the number of points on the axis
```

Re-reading the pairing permits collapse.

Holding the position alone does not hold the anchor. The same axis can be
kept exact while the relations change, and one way for them to change is
for two points to land on each other.

## Pointers

| Part of the statement | Run | Reading (rule m = h, tolerance 1e-6) |
|---|---|---|
| position free: the anchor is given up | `runs/anchor_generate_seed0.json.gz` | anchor still exact at the end on 6 of 8 5-cell starts and 3 of 8 16-cell starts |
| position held, pairing kept | `runs/constrained_generate_seed0.json.gz` | anchor held 8 of 8; reach h on 1 of 8 (5-cell) and 5 of 8 (16-cell) at 2500 trials |
| the pairing's type sorts the starts | same run, `anchor_points_on_axis` | 5-cell: type 1 reaches 5 on 1 of 4, type 3 on 0 of 4. 16-cell: type 2 on 4 of 4, type 0 on 1 of 2, type 6 on 0 of 2 |
| position held, pairing re-read | `runs/pairing_reread_seed0.json.gz` | reach h on 2 of 8 (5-cell) and 2 of 8 (16-cell); the anchor residual stays below 1e-6 on 8 of 8 |
| re-reading permits collapse | same run, `degenerate`; `runs/degeneracy_check_seed0.json` | 16-cell endpoints with coincident points: 6 of 8 re-read against 2 of 8 kept. 5-cell: 0 of 8 either way |
| choosing the anchor by its pairing type | `runs/anchor_selection_seed0.json.gz` | 5-cell: 4 of 8 reach 5 against 1 of 8 with the anchor as found. 16-cell: 5 of 8 against 5 of 8 |
| wrong rule (m = h - 1, h + 1) | all four runs | 0 starts reach h, every variant, every limit |
| no rule: fold the axis directions by the mirrors found, search the wedge, next axis becomes the next anchor | `runs/fold_loop_seed0.json.gz` | reach h on 1 of 8 (5-cell), 3 of 8 (16-cell), 0 of 2 (24-cell). The wrong-fold control and the no-fold loop end on the same count as the fold on every start; the fold reads about half the candidate axes the no-fold loop reads |

## Bounds on these readings

- 8 starts per shape, one seed, two shapes (the fold loop adds the 24-cell
  on 2 starts).
- In `fold_loop` the fold acts on the directions of candidate axes in the
  projection plane, not on the six rotation angles; that is one reading of
  a carried method its carrier had not tested.
- The preferred types in `anchor_selection` (5-cell 1, 16-cell 2) were read
  off the kept-pairing run on the same starts. Starts that are not
  re-anchored repeat that run's descents; only the re-anchored starts are
  new draws (4 of 8 for each shape).
- A re-anchored start spends several stage-1 pursuits (up to 6 here); the
  count is logged per start.
- "Pairing" is the nearest-point pairing read at one position.
- Counts are read at tolerance 1e-6. At 1e-4 several rows are higher; both
  are in the run files.
- The 400-trial limit binds in every variant. The 2500-trial readings are
  the ones quoted; they are still a bound.

## Not shown

Whether the pairing type is what decides the outcome. On the 16-cell every
selected anchor has 2 points on the axis and 3 of 8 still do not reach 8.
