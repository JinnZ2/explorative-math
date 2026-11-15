"""Tests for experiments/ring_tiling.py. CC0. Co-authored."""
import itertools
import random
import unittest

import _path  # noqa: F401
import emlib
import ring_tiling as R
from toolkit import primitives as P


class Pieces(unittest.TestCase):
    def test_inverse(self):
        inv = R.inverse([[2.0, 0.0], [0.0, 4.0]])
        self.assertEqual(inv, [[0.5, 0.0], [0.0, 0.25]])
        self.assertIsNone(R.inverse([[1.0, 2.0], [2.0, 4.0]]))

    def test_orthoplex_facets_cover_every_direction_once(self):
        orth = P.cross_polytope(5)
        facets = [[orth[2 * i + (0 if sg[i] > 0 else 1)] for i in range(5)]
                  for sg in itertools.product((1, -1), repeat=5)]
        rng = random.Random(3)
        samples = [P.unit([rng.gauss(0, 1) for _ in range(5)]) for _ in range(300)]
        cov, counts = R.coverage(samples, R.cones(facets))
        self.assertEqual(set(counts), {1})
        self.assertEqual(cov["covered_once"], 1.0)

    def test_dihedral_of_an_orthoplex_facet_is_90(self):
        orth = P.cross_polytope(5)
        self.assertEqual(R.dihedral([orth[0], orth[2], orth[4], orth[6], orth[8]]), [90.0])

    def test_logged_run(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = emlib.load_run(emlib.RUNS / "ring_tiling_seed0.json")
        self.assertIsNone(rec["ai_interest"])
        o = rec["outputs"]
        self.assertFalse(o["ring"]["ring_overlaps_itself"])
        self.assertEqual(o["regular_simplices_among_the_ring_vertices"]["n"], 24)
        self.assertTrue(o["regular_simplices_among_the_ring_vertices"]["all_are_ring_simplices"])
        self.assertEqual((o["copies"]["copies_tried"], o["copies"]["copies_placed"]), (47, 1))
        self.assertFalse(o["copies"]["tiles_the_sphere"])
        self.assertTrue(o["orthoplex_5"]["tiles_the_sphere"])
        self.assertFalse(o["dihedral"]["ring_simplex"]["whole_number"])
        self.assertTrue(o["dihedral"]["orthoplex_facet"]["whole_number"])
        small = emlib.load_run(emlib.RUNS / "ring_tiling_small_seed0.json")
        self.assertEqual(set(small["outputs"]), set(o))


if __name__ == "__main__":
    unittest.main()
