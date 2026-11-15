"""Tests for exact/volume.py and experiments/ten_vertex_volume.py. CC0. Co-authored."""
import unittest
from fractions import Fraction

import _path  # noqa: F401
import emlib
import ten_vertex_volume as T
from exact import PHI, ZPhi, ZERO
from exact import volume as XV


class ExactVolume(unittest.TestCase):
    def test_sqrt5(self):
        self.assertEqual(XV.SQRT5 * XV.SQRT5, ZPhi(5, 0))

    def test_unit_five_cell(self):
        five = XV.five_cell_unit()
        self.assertTrue(all(XV.dot(v, v) == ZPhi(1, 0) for v in five))
        self.assertEqual({XV.dot(five[i], five[j]) for i in range(5) for j in range(i + 1, 5)},
                         {ZPhi(Fraction(-1, 4), 0)})
        total = [sum((v[k] for v in five), ZERO) for k in range(4)]
        self.assertEqual(total, [ZERO] * 4)

    def test_determinant(self):
        one, zero = ZPhi(1, 0), ZERO
        eye = [[one if i == j else zero for j in range(4)] for i in range(4)]
        self.assertEqual(XV.det(eye), one)
        eye[0], eye[1] = eye[1], eye[0]
        self.assertEqual(XV.det(eye), -one)
        self.assertEqual(XV.det([[PHI, one], [one, PHI - 1]]), ZERO)     # phi (phi - 1) = 1

    def test_volume_of_the_five_cell(self):
        five = XV.five_cell_unit()
        facets = [tuple(j for j in range(5) if j != i) for i in range(5)]
        v = XV.volume_from_facets(five, facets)
        # regular 4-simplex with circumradius 1: 25 sqrt 5 / 384
        self.assertEqual(v, XV.SQRT5 * Fraction(25, 384))
        with self.assertRaises(ValueError):
            XV.volume_from_facets(five, [(0, 1, 2)])

    def test_bound_constant(self):
        # 25/384 (10 + 2 sqrt 5), written with sqrt 5 = 2 phi - 1
        self.assertEqual(T.BOUND, (ZPhi(10, 0) + XV.SQRT5 * 2) * Fraction(25, 384))


class Run(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = T.run()

    def test_a_misses_and_b_matches(self):
        a, b = self.out["A"], self.out["B"]
        self.assertEqual((a["verdict"], b["verdict"]), ("MISS", "MATCH"))
        self.assertTrue(a["below_the_bound"])
        self.assertFalse(a["above_the_bound"])
        self.assertTrue(b["equals_the_bound_exactly"])
        self.assertEqual(a["volume_exact_text"], "-25/64 + 25/32*phi")       # 25 sqrt 5 / 64
        self.assertEqual((a["n_facets"], b["n_facets"]), (30, 25))
        self.assertTrue(a["facets_close_up"] and b["facets_close_up"])

    def test_exact_and_float_volumes_agree(self):
        for k in ("A", "B"):
            e = self.out[k]
            self.assertAlmostEqual(e["volume_float_from_exact"], e["volume_float_direct"],
                                   places=10)
            self.assertLess(e["on_unit_sphere_max_deviation"], 1e-12)

    def test_a_is_the_sphere_compounds_set(self):
        self.assertTrue(self.out["A"]["same_set_as_the_sphere_compounds_entry"])
        self.assertTrue(self.out["A"]["five_cell_is_regular_and_unit_exactly"])

    def test_logged_run_and_reference(self):
        rec = emlib.load_run(emlib.RUNS / "ten_vertex_volume_seed0.json")
        self.assertIsNone(rec["ai_interest"])
        self.assertEqual(rec["outputs"]["A"]["verdict"], "MISS")
        self.assertEqual(rec["outputs"]["B"]["verdict"], "MATCH")
        text = (emlib.ROOT / "REFERENCES" / "related_studies.md").read_text()
        self.assertIn("2611.14674", text)
        self.assertIn("ABSTRACT READ", text)


if __name__ == "__main__":
    unittest.main()
