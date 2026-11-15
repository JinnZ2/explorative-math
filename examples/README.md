# Committed seed archives

These archives are executed seed runs, not a selection of interesting outcomes. All observer slots in all records are `null`. Every `observers.jsonl`, `divergence.jsonl`, and `notes.jsonl` starts empty. Example interest reports used in regression tests exist only in disposable temporary test directories, not these seed archives.

`geometry-seed` contains all four geometry families at dimensions 2, 8, and 32, with 8 pair/point samples, automatic 2*n packed directions, 8 trajectory updates, and seed 0. Random walks and iterated maps are separate runs. Every direction pair and every path state is kept, with explicit derived summaries.

`projection-seed` contains the default 5D crosspolytope, 12 sampled frames, phi composition, 1-degree tolerance, and seed 0. Every frame contains all four catalog comparisons, whether matches or misses. Catalog coincidence is not a natural-mechanism claim.

`simple-rules-seed` exhausts the declared circular width-2 CA space, binary substitutions with nonempty right-hand sides of maximum length 2 and seeds 0 and 1, and affine modulo-8 integer maps with a,b in [-1,1]. Every run retains its initial state plus 2 updates. Repeats, fixed points, and modulo-equivalent rules are retained.

Validate any archive without installing anything:

```sh
python -m explorative_math validate --dataset examples/geometry-seed
python -m explorative_math validate --dataset examples/projection-seed
python -m explorative_math validate --dataset examples/simple-rules-seed
```

For an independent repeat, choose a new directory instead of overwriting an example:

```sh
python -m explorative_math geometry --out runs/geometry-repeat \
  --dimensions 2,8,32 --samples 8 --steps 8 --seed 0
python -m explorative_math project --out runs/projection-repeat \
  --dimension 5 --shape crosspolytope --frames 12 --composition phi --seed 0
python -m explorative_math enumerate --out runs/rules-repeat --width 2 --steps 2
```
