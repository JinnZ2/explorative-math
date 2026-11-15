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

## R4. The inscribed 4-polytope of maximal volume with ten vertices

```text
where     arXiv:2611.14674 (math.MG; cross-list math.CO). Csaba Medgyes.
          Submitted 2026-11-10.
status    ABSTRACT READ on 2026-11-15, from the abs page. The copy fetched ends
          mid-sentence; the end of the abstract and the body are unread.
abstract  a 4-polytope with ten vertices inscribed in the unit sphere S^3 has
          volume at most 25/384 (10 + 2 sqrt 5). Up to orthogonal
          transformations, equality holds only for the convex hull of two
          regular pentagons inscribed in S^3 in orthogonal 2-planes through
          the origin. This extends the known optimum for six to nine
          vertices (two regular polygons in orthogonal planes). The proof
          uses sum-of-squares certificates in exact arithmetic to exclude
          the other simplicial combinatorial types, with a result of
          Horvath and Langi that an optimal polytope is simplicial.
carried   seen on the arXiv math/new listing.
```

Pointer to: `experiments/ten_vertex_volume.py`, `runs/ten_vertex_volume_seed0.json`.
Computed there in exact a + b*phi forms: two regular pentagons in orthogonal
planes have volume exactly 25/384 (10 + 2 sqrt 5); the hull of a regular
5-cell and its negative has volume 25 sqrt 5 / 64, about 0.8735, below it.
The theorem itself is not checked here; only those two volumes are.


## R5. Coxeter 1985

```text
where     Canadian Mathematical Bulletin 28(4):385-393, 1985.
          doi 10.4153/CMB-1985-045-5
status    NOT FETCHED here (doi.org and the publisher were outside the allowed
          network). Title: void. Read level CARRIED from a co-author:
          ABSTRACT, the full abstract, with its characteristic equation.
carried   the equation tan(n theta) = n tan(theta);
          the mapping: theta = half the turn, n = the simplex's vertex
          count, a helix in flat R^(n-1) with a translation along its axis.
          The abstract's text was not passed on and is not quoted here.
```

Pointer to: `experiments/simplicial_helix.py`, `runs/simplicial_helix_seed0.json`.
The run takes theta as half a plane's turn and n as the simplex's vertex
count, which is the carried mapping; under it the tetrahedral helix of R^3
(turn acos(-2/3), n = 4) satisfies the equation. The 17 closed chains do
not, on any plane. Carried note: that is the expected reading, since the
chains are spherical (closed, no translation) and the equation is for flat
helices.
