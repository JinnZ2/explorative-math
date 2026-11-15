# Two builds, one spec

Status:
- The comparison results below are **OBSERVED**: they are computed by `python3 -m compare` and written to `compare/runs/*.jsonl`.
- Readings of them are **DERIVED** where the reasoning is shown, and **PROPOSED** where they are next steps.
- Not established: which build's choices are better. Nothing here ranks them.

Co-authored: **JinnZ2, Manus AI, Claude**. CC0.

## Purpose

Two observers built one spec independently:
- **Build A**, by **Manus AI**: `explorative_math/`, `examples/`, `docs/schema.md`, `docs/methods.md`.
- **Build B**, by **Claude** (Claude Code): `toolkit/`, `experiments/`, `emlib.py`, `runs/*.json`, `open_questions.md`.

Both started from JinnZ2's direction. They are kept side by side **on purpose**: where they differ, the difference is data about how each observer read the spec.

Neither build is unified with, refactored into, or deduplicated against the other. Each keeps its own rules:
- Build A: `CONTRIBUTING.md` and `docs/schema.md`.
- Build B: `CLAUDE.md` and `README.claude-code-build.md`.

`compare/` feeds the same inputs to both and logs every field.

## Design deltas

| Axis | Build A (Manus) | Build B (Claude Code) |
|---|---|---|
| Centre of gravity | **archive, integrity, refusal**: manifest with exact counts, validation replay, explicit refusal of oversized shapes, append-only logs, no silent truncation | **toolkit, generator, puzzles**: composable primitives / composers / probes, a random experiment generator, an open-questions puzzle layer |
| Run storage | a dataset directory: `manifest.json`, `records.jsonl`, `observers.jsonl`, `divergence.jsonl`, `notes.jsonl`; output paths must be new | one JSON file per run in `runs/`, overwritten on re-run with the same seed |
| Run schema keys | `id, params, outputs, human_interest, ai_interest, compressibility, divergence_note, status` | `params, outputs, angles, human_interest, ai_interest, compressibility, nature_match, divergence_note` |
| `ai_interest` | always `null` in the base record; filled only by an explicit, attributed observer event; "AI-generated code is not an AI interest judgment" | `null` in automatic runs; in the five hand-made experiment runs (a–e) it holds the AI's own reading, labelled as such |
| Interest value type | boolean or `null` | free text or `null` |
| `compressibility` | `null` unless an observer reports a number with a named method | computed automatically (LZ76 phrase count + zlib ratio) with a stated definition |
| Divergence | derived events in `divergence.jsonl`, when both interest slots are explicit booleans that disagree | `divergence_note` slot plus `DIVERGENCE_LOG.md`, a hand-filled table |
| Status labels | `status` key on each run and each output | no `status` key; OBSERVED / DERIVED / PROPOSED appear in `ai_interest` text and the README |
| Nature catalog | `nature_catalog.json`: target angles, source URLs, `DERIVED`, every match and miss kept | `nature_catalog.py`: three entries marked CARRIED; a run names at most one match |
| Random streams | each run seeds `random.Random(sha256(canonical params))` | `random.Random(seed)` directly |
| Polytope normalisation | centred, **unit circumradius** | raw coordinates (simplex edge √2; hypercube ±1; cross-polytope ±1) |
| φ composition | copy scaled by φ **and** rotated in-plane by the golden angle (A states that uniform scaling alone leaves angles unchanged) | `phi_stack`: copies scaled by φ, **not** rotated |
| Implemented, not in the other | superposition packing, random walk, iterated affine map, finite simple-rule enumeration (CA, substitutions, affine), CLI, zipapp phone build | cut-and-project, Fibonacci chain, golden-angle packing, symmetry / mirror probes, 24-cell, icosahedron, lattices |
| Cut-and-project | documented as **not implemented** (`docs/methods.md`) | implemented for Z² with a strip window |

## First comparison results (OBSERVED, seed 0)

Each `compare/runs/*.jsonl` file counts fields by verdict.
- AGREE means exactly equal.
- DIFFER always carries its absolute and relative difference.
- No tolerance is applied.

| Comparison | AGREE | DIFFER | ONLY_IN_A | ONLY_IN_B | max abs diff |
|---|---|---|---|---|---|
| constants | 1 | 1 | 0 | 0 | 2.8e−14 |
| shell_formula | 15 | 21 | 0 | 0 | 3.3e−16 |
| seed_alignment | 0 | 168 | 0 | 0 | 1.8 |
| shared_stream_sampler | 340 | 0 | 0 | 0 | 0 |
| orthogonality_stats | 0 | 5 | 10 | 10 | 7.1 |
| polytopes | 30 | 15 | 0 | 3 | 1.24 |
| rotate_project | 58 | 0 | 0 | 0 | 0 |
| angle_measure | 29 | 0 | 6 | 6 | 0 |
| phi_composition | 51 | 77 | 8 | 0 | 180 |
| capabilities | 0 | 0 | 9 | 8 | — |
| schema | 6 | 5 | 2 | 2 | — |

**Readings (DERIVED):**
- **Seeds do not align.** The same seed gives different vectors in the two builds: 0 of 168 coordinates agree. Build A hashes the run's parameters into its seed and Build B uses the seed as given. So "seed 0" names two different samples, and the orthogonality statistics differ by up to 7.1° at d=2, where 400 samples spread most.
- **The samplers themselves agree.** Given the same explicit `random.Random(s)`, both unit-vector samplers return identical vectors (340 of 340 coordinates). Build A's `math.fsum` norm and Build B's `sum` norm gave the same floats here. The seed difference is entirely the seeding rule.
- **Shared explicit rotations give identical geometry.** Rotating Build A's vertices by the same plane-angle list through each build's rotator gives identical projected points (58 of 58). Build A takes degrees and Build B radians.
- **The constants are equal in exact arithmetic but not in float.** The golden angle differs by 2.8e−14° (A: 360/φ²; B: 360(1−1/φ)). The shell formula differs in 21 of 36 cells, by up to 3.3e−16 (A: `−expm1(n·log1p(−t))`; B: `1−(1−t)ⁿ`).
- **The polytopes are the same shapes at different sizes.** Vertex counts and the rotation- and scale-invariant distance sets all agree. Circumradius differs because A normalises to 1. The cross-polytope's vertex order also differs: A lists −e then +e, B lists +e then −e.
- **"φ composition" means two different operations.** On Build A's own projected frames, A's union carries about twice as many distinct directions as B's (for example 16 against 8 for the cross-polytope). That is because A rotates the copy and B does not.
  - Build B's symmetry probe gives the same k for both unions on every frame tried. A rotated copy of a k-fold set keeps the same k.
  - On simplex frame 3, A's catalog matched "capsids" and "quasicrystals" (the same 72° target, as A's docs note).
- **Angle measurement agreed on these shapes.** A measures about the origin, B about the centroid. For centred shapes they agreed exactly here; they would part for a shape not centred on the origin.

**Schema (DERIVED):** the two builds disagree on whether an AI reading can sit in `ai_interest` in a base record.
- Build A says no: interest comes only from an explicit observer event.
- Build B writes one into its five hand-made runs.

This is recorded as a disagreement and **not resolved**; neither build's rule is chosen.

## Conflict decisions made at merge (OBSERVED)

The branch histories were unrelated, so they were merged with `--allow-unrelated-histories`. Three files conflicted:

| File | Decision |
|---|---|
| `README.md` | Build A's README kept unchanged under a short label pointing to both builds. Build B's README kept unchanged as `README.claude-code-build.md`. |
| `LICENSE` | Build A's full CC0 legal code kept as `LICENSE`. Build B's short dedication kept as `LICENSE.claude-code-build.txt`. |
| `.gitignore` | Both sections kept and labelled. One rule edited: Build A's `runs/` became `runs/*` plus `!runs/*.json`. Build A ignores `runs/` (datasets are local) and Build B commits `runs/*.json`. With this rule, A's dataset **directories** under `runs/` stay ignored and B's run files stay tracked. |

Other shared files:
- **`CLAUDE.md`** (Build B's file) gained a header saying its rules are Build B's, and naming where Build A's rules live.
- **`AUTHORS.md`** (Build A's file) gained a section for Build B and the comparison. Nothing in it was removed.

The test suites of both builds pass together after the merge.

## Open (PROPOSED)

- Run Build A's own geometry and projection experiments with explicit seeds chosen so that a Build B run uses the same stream. That needs an adapter that exposes A's per-run seed; none is written.
- Compare cut-and-project and golden-angle outputs once Build A implements them. Until then they are ONLY_IN_B.
- Neither build has a principal-angle solver, so there is nothing to compare there.
- See `open_questions.md` Q11.
