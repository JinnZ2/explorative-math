# Projection methods and nature-catalog conventions

The implemented pipeline samples a finite set of Givens rotations of explicit shapes, projects onto two coordinates, records polar angles and cyclic angular gaps, composes a rotated copy, and compares every catalog target. It preserves both matches and misses. It does not treat a matching angle as a physical explanation.

## Established method names

A **grand tour**, introduced by Daniel Asimov in 1985, visualizes multivariate data through a sequence of 2D projection subspaces designed to be dense in the possible subspaces. Density is a limiting property, not a finite exhaustive sweep. This project's finite sampled rotations do not claim that convergence property. [1]

**Projection pursuit** searches for projections that score highly on an explicit structure index, such as non-Gaussianity. Choosing a visually persuasive projection from high-dimensional data can find apparently interesting structure even in Gaussian data. This project records a sampled tour; it does not optimize a pursuit index. [2]

A **Petrie polygon** of a regular polyhedron follows edges so consecutive pairs, but not consecutive triples, belong to one face. A specially chosen **Coxeter plane** makes a Coxeter element act as a rotation and can reveal regular-polytope symmetry. An arbitrary 2D projection is not automatically a Coxeter-plane construction. These are future experiment directions, not implemented modes. [3]

**Cut-and-project** selects lattice points using an acceptance rule/window before projecting them. De Bruijn's Penrose construction selects points of a 5D cubic lattice; projecting arbitrary lattice points is not sufficient to create a Penrose tiling. This repository does not implement the acceptance window, tiling topology, or a Penrose generator. [4]

The **Grassmannian** is the space of fixed-dimensional linear subspaces. The principal angles between two subspaces are obtained from the singular values of the product of their orthonormal basis matrices. They are basis-independent subspace comparisons, not the polar gaps between projected vertices computed here. A principal-angle solver is a proposed extension. [5]

## Nature catalog

Each catalog entry has a target angle, source, `DERIVED` status, and an explicit interpretation. The targets use different natural conventions. Comparing all of them against projected polar gaps is an exploratory numerical search, **not a like-for-like natural measurement or identification test**.

**Phyllotaxis:** the idealized golden angle is `360/phi² = 137.507764…` degrees. It describes the divergence between successively originating florets in a common sunflower model. Real heads can depart from the ideal. A central gap in a projected polytope is not the same growth measurement. [6]

**Virus capsids:** many capsids have icosahedral symmetry, which includes fivefold axes. A fivefold rotation is `360/5 = 72` degrees. It is not a universal angle between proteins, and not all viruses are icosahedral. A projected 72-degree gap does not establish capsid architecture or a biological mechanism. [7]

**Quasicrystals:** Shechtman's alloy discovery demonstrated fivefold symmetry and long-range order without ordinary periodicity; he received the 2011 Chemistry Nobel. The catalog uses 72 degrees only as a fivefold signature. Other quasicrystals have other symmetries. A single matching angle cannot demonstrate quasiperiodic atomic order, and phi composition alone is not a quasicrystal construction. [8] [9]

**Hexagonal patterns:** 60 degrees is the central angular spacing of a regular hexagon's vertices. Honeycomb geometry and many basalt columns show hexagonal patterns, but the formation mechanisms differ. Cooling fractures in basalt have approximately 120-degree junction angles, and column cross-sections can have three to seven sides. The catalog's 60-degree target is not that fracture-junction measurement and does not establish optimal disk packing. [10]

## Measurements, composition, and misses

The projection origin is the shape's centered origin. Zero-radius projected points have no defined polar angle; they remain in the recorded points with an explicit degeneracy. Coincident angular directions create zero gaps and are retained. The gaps sum to a full turn for a nonempty circular ordering.

Pure positive uniform scaling changes lengths, not angles. The `phi` mode therefore distinguishes the unchanged pure-scaled angles from the angular gaps of a union with a second copy that is explicitly rotated in the projection plane. This is a synthetic construction whose rotation is recorded, not a discovered natural relation.

Each frame compares every target with the measured post-composition gaps and records the nearest measurement, absolute error, threshold, and a match or miss. An empty usable angle set is a recorded miss, not an omitted comparison. Targets sharing an angle, such as the two 72-degree entries, are not independent evidence. Broad thresholds, many frames, and many targets increase opportunities for coincidental matches. There is no automatic significance claim, interest rating, or preferred projection selection.

## References

[1]: https://epubs.siam.org/doi/10.1137/0906011 "The Grand Tour: A Tool for Viewing Multidimensional Data — Daniel Asimov, 1985"
[2]: https://www.pnas.org/doi/10.1073/pnas.1801177115 "Projection pursuit in high dimensions"
[3]: https://sites.math.washington.edu/~billey/classes/reflection.groups/references/Coxeter.1934.pdf "Discrete Groups Generated by Reflections — H. S. M. Coxeter, 1934"
[4]: https://www.ams.org/publicoutreach/feature-column/fcarc-ribbons "Penrose Tilings Tied up in Ribbons — David Austin, American Mathematical Society"
[5]: https://arxiv.org/html/physics/9806030 "The Geometry of Algorithms with Orthogonality Constraints — Edelman, Arias, and Smith"
[6]: https://beatymuseum.ubc.ca/2021/11/17/fibonacci-day/ "November 23 is Fibonacci Day! — Beaty Biodiversity Museum"
[7]: https://pmc.ncbi.nlm.nih.gov/articles/PMC524849/ "Origin of icosahedral symmetry in viruses"
[8]: https://www.nist.gov/nist-and-nobel/dan-shechtman "Dan Shechtman — NIST"
[9]: https://euler.phys.cmu.edu/widom/pubs/html/EB_QC/QCarticle.html "Quasicrystals — Michael Widom"
[10]: https://www.nps.gov/subjects/volcanoes/columnar-jointing.htm "Columnar Jointing — U.S. National Park Service"
