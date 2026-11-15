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
- **Y** — an independent result from a run, quoted.

| id | H | A | C | Y (quoted) | state | relation (Claude's reading) |
|---|---|---|---|---|---|---|
| E1 | drift | two builds kept side by side | drift → two builds | `\| seed_alignment \| 0 \| 168 \| ...`; `\| rotate_project \| 58 \| 0 \| ...` (docs/two_builds.md) | Y_PRESENT | DOES_NOT_TEST: a design choice; no run can come out against it |
| E2 | OR parent | ancestry run; "void cover = dual OR parent" | "OR parent" → ancestry run | "The parent's dual covers the deepest voids in 7 of 8 forms."; "The parent, in its own position, covers no void class of any form." (docs/complement.md) | Y_PRESENT | CONTRADICTS the "OR parent" clause in the tested frame |
| E3 | φ-composition | welded-term entry; Q12 | φ-composition → weld | 16 against 8 distinct directions (open_questions.md Q12); unions equal on all frames at all bounds (resolution_leak 2b); probe separates 8- from 16-gon at k_max 24 (2a) | Y_PRESENT | MIXED: the operations differ; the symmetry probe cannot see it |
| E4 | tetra-chain | E8 link | tetra-chain → E8 link | `"n_tetrahedra_to_close": 30` (runs/corkscrew_seed0.json) is about the chain, not E8 | **UNRESOLVED** | no run tests the E8 link |

Not recorded: anything about the episodes that is not in this repository.
