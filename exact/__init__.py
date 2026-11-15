"""exact -- numbers kept exact, no floats. CC0. Co-authored. Stdlib only.

  ZPhi   a + b*phi with a, b Fractions; phi^2 = phi + 1
  Cyc5   u + v*w with u, v ZPhi and w^2 = -(2 + phi): the plane, exactly,
         for anything built on fifth roots of unity (pentagon, Penrose)
  fold   reflect a point set into the wedge of its mirrors, with
         multiplicities; exact equality, no tolerance
"""
from .zphi import ZPhi, PHI, ONE, ZERO  # noqa: F401
from .cyc5 import Cyc5, ZETA, W, cross, dot, norm2  # noqa: F401
from . import inflation  # noqa: F401,E402
from . import fold  # noqa: F401,E402
from . import complement  # noqa: F401,E402
