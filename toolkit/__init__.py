"""toolkit -- the reusable layer under the experiments. CC0. Co-authored. Stdlib only.

  primitives  vectors, plane rotations in any dimension, projection, angles,
              polytope vertex generators, lattices
  composers   rotate_project, phi_scale / phi_stack, cut_and_project,
              iterate, overlay
  probes      rotational / mirror symmetry, periodicity, zlib ratio,
              angle histogram
  generator   random combination of composers -> a new experiment, logged
  loghook     every run writes the schema; interest slots void
"""
