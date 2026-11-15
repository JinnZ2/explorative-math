# Complement readings

Status: notes on what the collapse / complement runs do and do not show.
Every statement points at a run. Carried names are marked; nothing here was
checked against an outside source. Claude Code build only. Co-authored. CC0.

Runs: `collapse_gap`, `complement_fit`, `complement_fit_v2`,
`sphere_compounds` (all under `runs/`, seed 0).

## Plane: what the reading gives

- The complement's distance field has the same rotational orders and the
  same number of mirrors as the shape's own points, for all 12 shapes read,
  and none for the two random controls (`collapse_gap`).
- The same shape reached by another route is the nearest probe at field
  difference 0 (5-cell on its Coxeter plane and the 5-gon; 16-cell and the
  8-gon).

## Plane: the fit score is weak at low N

Note. At an equal point budget the 2D complement fit does not separate
shapes from random point sets (`complement_fit_v2`, `equal_budget`):

```text
N = 3   baseline 5-cell:  random_8 0.2771, tesseract 0.2775, random_5 0.2776, 24-cell 0.2813
        baseline 16-cell: random_8 0.2481, random_5 0.2486, tesseract 0.2499, 24-cell 0.2524
N = 5   the first three are tesseract, 24-cell and random_8 within 0.002 of each other
N = 8   tesseract and 24-cell lead random_8 by 0.02 to 0.03
```

At N = 3 a random control is first for both baselines, and the regular
polygons rank below both controls. Without the budget the ranking follows
the number of points (`complement_fit`). So in the plane this score reads
mostly how many points are placed and where three or five of them can sit;
it is not a reading of the probe's shape at these N. Bound: the N points
are chosen greedily at the whole probe's placement; two probe sizes; one
seed for the controls.

## Sphere: pairs and what they join into

```text
pair (placement)                       joined set                          how it was reached
16-cell + tesseract (as built)         24-cell                             as built
16-cell + 16-cell                      tesseract                           built in exact fractions; also searched (60.000000 degrees)
5-cell + 5-cell                        a 5-cell and its negative, 10       searched (75.52 degrees)
24-cell + 24-cell                      a 24-cell and its dual, 48          built; also searched (44.997 degrees)
```

Two names are CARRIED from a co-author. What was computed for each set
(`sphere_compounds`, convex hull by brute force):

```text
"compound of two 5-cells"   10 vertices; its convex hull, carried name "disphenoidal
                            30-cell (dual of bitruncated 5-cell)": 40 edges, 60 faces,
                            30 cells; every cell the same tetragonal disphenoid
                            (angles 75.52 x 4, 104.48 x 2)
"disphenoidal 288-cell"     48 vertices, 336 edges, 576 faces, 288 cells; every cell the same
                            tetragonal disphenoid (angles 45 x 4, 60 x 2)
```

The 48-direction set has 288 cells, each a disphenoid, which is what the
carried name says. Whether these are the objects those names denote in the
literature was not checked: the reference pages were outside the allowed
network. A compound is not its convex hull; the 30 cells above are the hull
of the ten directions, and they are 30 disphenoids, which is what the name
carried for that hull says. That the hull is the dual of the bitruncated
5-cell was not computed here.

Bounds: the rotation search is a lower bound on the separation; the
catalogue the joined sets are compared with is short (eight entries); a
match is equality of sorted pairwise angles within 0.5 degrees.

## Pattern (PROPOSED): the best-separated partner is the dual

PROPOSED. Among a set of probe shapes, the one that can be placed farthest
from a baseline is the baseline's dual; placed there it sits in the
baseline's deepest voids and covers them, and the two together are a known
figure. Three cases, from `sphere_compounds` and `dual_partner`:

```text
baseline    dual partner (placement)     separation   joined set
16-cell     tesseract                    60           24-cell (24 directions)
5-cell      5-cell, negated              75.52        ten directions; hull carried as "disphenoidal 30-cell"
24-cell     24-cell, dual position       45           48 directions; carried as "disphenoidal 288-cell"
```

What the test logged (`runs/dual_partner_seed0.json.gz`):

```text
sphere, 4 baselines x 4 probes
  dual is best-separated                4 of 4   (alone in 2; tied in 2: 16-cell with itself at 60,
                                                  24-cell with the 16-cell and tesseract at 45)
  dual sits in and covers the voids     4 of 4, and no other probe does for any baseline
plane, 7 regular polygons
  dual (the same polygon, turned)       best-separated in 7 of 7; alone in 4, tied in 3
                                        (6-gon with the 3-gon, 8-gon with the 4-gon, 12-gon with three)
controls (dual not in the probe set)
  random 8 directions                   no probe covers its voids
  5-cell with its negative              no probe covers its voids
  rectified 5-cell                      the 5-cell covers its voids
  9-gon, 10-gon                         best-separated: 3-gon, 5-gon
```

Not shown. Separation alone does not pick the dual: it ties in 5 of 11
baselines, and a probe with fewer points can always sit farther. The void
reading is what separates the dual on the sphere. One control (the
rectified 5-cell) has a covering probe that is not its dual in this set, so
"covers the voids" is not by itself a test for duality. Bound: four
polytopes and seven polygons as probes; the rotation search is a lower
bound; duals are taken as standard, not derived here; one seed.

## Pointer, structural only

Alexander duality relates a set to its complement. No homology is computed
in these runs.
