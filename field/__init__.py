"""field -- an ideal lattice under an environmental field. CC0. Co-authored. Stdlib only.

  hexsheet   a sheet of hexagonal cells: exact ideal reference, a field over
             the sheet, a local relaxation rule, and the per-cell reading of
             how each cell answers the field (its RESPONSE).

Words, as used in this package:
  ideal      the REFERENCE sheet (regular hexagons, spacing 1). It is what
             a cell is measured against. It is not a target the cells are
             meant to reach.
  response   how a cell differs from the reference under the field where it
             sits. It is not called error.
  optimum    relative to the environment: which sheet carries less strain
             depends on the field it sits in. The reference sheet is the
             unstrained one only where the field is zero.

Everything here is CONSTRUCTED: a toy rule on a toy sheet. Its numbers are
readings of this model, not findings about any material.
"""
from . import hexsheet  # noqa: F401
