# Related studies

Pointers, not results. Nothing here was computed in this repository.
Co-authored. CC0.

Each entry says how far it was read. Three levels are used:

```text
ABSTRACT READ   the abstract page was fetched on 2026-11-15; the body was not read
NOT FETCHED     the page could not be reached from the build environment
CARRIED         a statement about the paper supplied by a co-author, not checked
                against the paper here
```

A line under "carried" is a reason to read the paper. It is not a citation
for the statement.

---

## R1. Topological charge vectors in plasmonic quasicrystals

```text
where     Science, 2025. doi 10.1126/science.adt2495
status    NOT FETCHED (doi.org, science.org and crossref were outside the
          allowed network). Title, authors and abstract are void here.
carried   topological charge vectors in plasmonic quasicrystals;
          the topology survives projection;
          some cases are distinguishable only through the source dimension.
```

Pointer to: `open_questions.md` Q13.

## R2. Spacetime Quasicrystals

```text
where     arXiv:2601.07769 (hep-th). Latham Boyle, Sotirios Mygdalas.
          Submitted 2026-01-12, v2 2026-02-12.
status    ABSTRACT READ
abstract  generalises self-similar quasicrystals (Penrose, Ammann-Beenker)
          from Euclidean space to Minkowski spacetime; constructs Lorentzian
          quasicrystals; names non-crystallographic orientational symmetry
          and discrete scale invariance.
carried   symmetric cut-and-project; a weighting-function window.
          The abstract uses neither "cut-and-project" nor "window".
          Whether the body does is unread.
```

Pointer to: `toolkit/composers.py` `cut_and_project(weight=...)`,
`experiments/weighted_window.py`. The Gaussian weight implemented here was
not taken from the paper's text; its form in the paper, if any, is void.

## R3. Holographic Foliations: Self-Similar Quasicrystals from Hyperbolic Honeycombs

```text
where     arXiv:2408.15316 (hep-th). Latham Boyle, Justin Kulp.
          Phys. Rev. D 111, 046001 (2025), doi 10.1103/PhysRevD.111.046001.
status    ABSTRACT READ (was: snippet only, UNVERIFIED). Body unread.
abstract  a regular tessellation of (d+1)-dimensional hyperbolic space gives
          a d-dimensional boundary with self-similar quasicrystalline
          properties: a local, invertible, self-similar substitution tiling.
          Covers the regular 3D hyperbolic tessellations; in {3,5,3} it
          reports a 2D boundary quasicrystal with 5-fold points that is not
          the Penrose tiling.
```

Pointer to: nothing in the toolkit yet. No hyperbolic primitive exists here.
