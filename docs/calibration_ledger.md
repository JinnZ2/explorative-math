# Calibration ledger (retrospective, 2026-11-15)

Status: retrospective record of four episodes. **Y is a run result only**,
quoted verbatim from a file in this repository and checked by
`tests/test_calibration_ledger.py`. An episode with no run is UNRESOLVED.
Data: `docs/calibration_ledger.json`. CC0, co-authored.

**Interested party.** Claude wrote this ledger and is a party to every episode:
the A expressions are its own, and the relation column is its reading of its
own work. Recheck the quoted Y lines rather than the relation.

Fields:
- **H** — the human's mark, as the operator named the episode (CARRIED from the
  operator's message; the human's words in the chat are not in this repository).
- **A** — the AI's expression of it, and where it landed.
- **C** — the coupling. Only the operator's arrow is recorded; anything more is
  UNDECLARED.
- **Y** — an independent result from a run: a quoted line (at a git ref when
  the file is on another branch), or a verdict count recomputed by the test
  from a `compare/` jsonl.

| id | H | A | C | Y (quoted) | state | relation (Claude's reading) |
|---|---|---|---|---|---|---|
| E1 | drift (read then as a defect) | two builds kept side by side | drift → two builds | rotate_project AGREE 58 of 58; shared_stream_sampler AGREE 340 of 340; seed_alignment DIFFER 168 of 168 (compare/runs, recounted by the test); "The seed difference is entirely the seeding rule." (docs/two_builds.md) | Y_PRESENT | CONTRADICTS drift-as-defect: under shared inputs the builds agree exactly; the deltas read trace to conventions. Not conventions: φ composition (E3) and the ai_interest schema rule (recorded, unresolved) |
| E2 | OR parent | ancestry run; "void cover = dual OR parent" | "OR parent" → ancestry run | "The parent's dual covers the deepest voids in 7 of 8 forms."; "The parent, in its own position, covers no void class of any form." (docs/complement.md) | Y_PRESENT | CONTRADICTS the "OR parent" clause in the tested frame |
| E3 | φ-composition | welded-term entry; Q12 | φ-composition → weld | 16 against 8 distinct directions (open_questions.md Q12); unions equal on all frames at all bounds (resolution_leak 2b); probe separates 8- from 16-gon at k_max 24 (2a) | Y_PRESENT | MIXED: the operations differ; the symmetry probe cannot see it |
| E4 | tetra-chain | E8 link | tetra-chain → E8 link | "the (7, 13)/30 chain is E8's Coxeter element on the conjugate 4-space"; `"orbits_that_are_the_logged_chain": 2` of 8 orbits (runs/e8_h4_fold_seed0.json at a2b6ac0, branch claude/quasicrystal-projection-properties-448bx7) | Y_PRESENT | SUPPORTS, scoped: one ordering of the reflections run; 6 of 8 orbits are not the chain |

## Corrections (kept; `corrections` in the JSON holds the prior entries)

- **E1**, 2026-11-15, by the operator: first read as DOES_NOT_TEST on the
  summary table. The compare results test the defect reading directly.
- **E4**, 2026-11-15, by the operator: first entered as UNRESOLVED ("no run tests
  the E8 link"). That search covered main only; the run is on another branch.

Not recorded: anything about the episodes that is not in this repository.
