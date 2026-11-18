> **Two builds of one spec live in this repository, on purpose.** This README is the
> **Manus build**'s (`explorative_math/`, `examples/`, `docs/schema.md`), unchanged below.
> The **Claude Code build** (`toolkit/`, `experiments/`, `emlib.py`, `runs/*.json`,
> `open_questions.md`) has its own README: [README.claude-code-build.md](README.claude-code-build.md).
> How they differ, and the cross-build comparison: [docs/two_builds.md](docs/two_builds.md), `compare/`.
> Neither build is folded into the other; the difference between them is data.

# explorative-math

**Purpose: map whose curiosity lands where.**

Co-authored by **JinnZ2 and Manus AI**. CC0 1.0 Universal. Python 3.10+ standard library only. No pip, compiler, model API, or network connection is needed to run or build it.

This is an archive of exploration, not a feed of selected discoveries. Keep every requested output, including constants, repetitions, degeneracies, uninteresting cases, and comparison misses. Human interest and AI interest are separate observations, not a ranking built into the generator.

## Start on a computer or phone

From a downloaded or cloned repository, in any terminal with Python 3.10+:

```sh
python -m explorative_math --help
python -m explorative_math geometry --out runs/geometry \
  --dimensions 2,8,32 --samples 32 --steps 32 --seed 0
python -m explorative_math validate --dataset runs/geometry
```

A phone with a Python 3.10+ environment can run the same source. You can also build a single-file Python application using only the standard library:

```sh
python tools/build_phone.py
python dist/explorative-math.pyz --help
python dist/explorative-math.pyz geometry --out runs/phone \
  --dimensions 2,8 --samples 16 --steps 16
```

This is a Python application, not a native iOS app. It still requires a Python environment that supports standard-library modules and local files. Use smaller bounds on a phone; nothing is silently dropped to fit memory. Each dataset path must be new.

## Seed region: high-dimensional geometry

`geometry --experiment all` runs four experiment families across every requested dimension:

1. **Near-orthogonality.** Save every pair of independently sampled unit vectors, dot product, and angle. Compare the full empirical sample with the derived expectation `E[dot²] = 1/n`.
2. **Shell concentration.** Sample the unit ball uniformly by volume, keeping every point and its radius. Compare measured shell membership with the derived fraction `1 - (1 - thickness)^n`.
3. **Superposition.** Sample `k > n` directions, keep all vectors and every unordered pair angle, and report the minimum angle. This is a random packing baseline, not an optimized packing or a demonstration of a learning system's superposition. `--vectors 0` means `k = 2*n`.
4. **Trajectories.** Keep the initial point and every subsequent state for a unit-step random walk and an iterated affine map. Repeated points are retained.

Samples are finite and seeded. They are not an exhaustive enumeration of continuous geometry. Seeded runs are reproducible in the same implementation/runtime; cross-platform floating-point details can differ.

## Shape → rotate → project → measure → compose → compare

```sh
python -m explorative_math project --out runs/projections \
  --dimension 5 --shape crosspolytope --frames 12 \
  --composition phi --tolerance-degrees 1 --seed 0
```

Choose a crosspolytope, hypercube, or regular simplex. Keep the full shape, every sampled rotation, its 2D shadow, polar angles, and cyclic angular gaps. Degenerate projected points are recorded, not removed from the archive. Compare the post-composition measurements with every nature-catalog target, saving **both matches and misses**, their error, tolerance, and references.

`phi` composition uses a second, golden-ratio-scaled copy with an explicit in-plane rotation. **Uniform scaling alone does not change angles.** The union with the rotated second shape can change angular gaps. `none` and `second` are also supported. The vertex bound causes an explicit refusal, never a partial shape.

`--max-vertices` counts vertices across both shapes together. `second` uses a unit-scale copy of the same chosen shape, with the recorded in-plane rotation.

The catalog starts with the golden angle in phyllotaxis, a fivefold rotational signature relevant to many icosahedral virus capsids, a fivefold quasicrystal/Penrose reference, and hexagonal geometry. An angular coincidence is **not** evidence of a shared natural mechanism. Multiple comparisons and projection choices can produce chance matches. See [catalog and method notes](docs/methods.md) for references and qualifications.

The sampled rotation tour is related to grand-tour visualization. It does **not** sweep all projections. Projection pursuit, Coxeter-plane/Petrie projections, actual cut-and-project Penrose construction, and Grassmannian principal-angle comparison are documented directions, not falsely advertised implementations.

## Per-run schema and blank observers

Each line of `records.jsonl` contains one run:

```json
{
  "id": "stable-run-id",
  "params": {},
  "outputs": [],
  "human_interest": null,
  "ai_interest": null,
  "compressibility": null,
  "divergence_note": null,
  "status": "OBSERVED"
}
```

Real runs contain full parameters and outputs. All four observer fields remain `null` in the immutable base archive. No human or AI preference is inferred, and no compression metric is calculated automatically. `null` means unknown, not uninterested.

`OBSERVED` identifies executed outputs or an explicitly attributed observer report. `DERIVED` identifies a calculation from recorded values or a stated mathematical relation. `PROPOSED` identifies an untested suggestion. These are provenance labels, not truth certificates. Geometry summaries and comparison results carry their own `DERIVED` labels inside the observed run.

To inspect the schema and audit protocol, read [docs/schema.md](docs/schema.md). An explicit report can be added later; **these examples are commands, not pre-filled judgments**:

```sh
python -m explorative_math observe --dataset runs/geometry \
  --entry RUN_ID --slot human_interest --value true \
  --observer "your-name" --reason "your actual assessment"
python -m explorative_math show --dataset runs/geometry --entry RUN_ID
```

Interest values are `true`, `false`, or `null`. Numeric compressibility needs a named `--method`; there is no universal interpretation of its scale. `divergence_note` accepts explicit text or `null` and is never automatically filled.

`observers.jsonl` is append-only. `show` materializes the latest reports without modifying the original run. `divergence.jsonl` starts empty and logs derived events when **both** interest slots are explicitly filled and disagree. It preserves openings and later resolutions; an unknown value is never treated as disagreement. Nature matches are separate from observer divergence.

## Secondary seed region: simple rules

The original simple-rule explorers remain available:

```sh
python -m explorative_math enumerate --family all --out runs/rules \
  --width 4 --steps 4
```

They exhaust the **declared finite spaces**: all 256 elementary cellular-automaton rules on a circular binary row; every binary substitution with nonempty right-hand sides up to the configured length; and affine integer maps `(a*x+b) % modulus` over inclusive coefficient bounds and all residue seeds. Every initial state and requested update is saved, even after a fixed point or cycle is reached. Equivalent outputs and syntactically distinct rules are not deduplicated.

## Open and falsified (Claude Code build)

- H-ORDER falsified 2026-11-18 (cross-tab in runs/h_order_seed0.json.gz); total-positive embedding count: open

## Integrity, tests, and license

The manifest declares bounds and exact run/output counts. Normal completion is marked explicitly; caught failures leave an error record when possible and an incomplete manifest. Validation checks archive structure, replayed coverage for supported generators, and observer/divergence consistency. Complete runs are appended durably; an abrupt kill can lose the in-memory unfinished run and may leave a lock requiring inspection. This is not a resumable distributed database.

The [committed seed archives](examples/README.md) include geometry, the complete sampled projection pipeline, and finite simple-rule enumeration, with all observer slots blank.

```sh
python -m unittest discover -s tests -v
```

An optional GitHub Actions template is in [`tools/ci-template.yml`](tools/ci-template.yml). To enable CI, copy it to `.github/workflows/test.yml` using an account with workflow-write permission. The initial repository does not enable that workflow automatically.

See [CONTRIBUTING.md](CONTRIBUTING.md), [AUTHORS.md](AUTHORS.md), and the full [CC0 license](LICENSE). The public-domain dedication covers the original code, documentation, tests, and committed seed data.
