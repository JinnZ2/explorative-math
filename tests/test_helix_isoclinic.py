"""Tests for experiments/helix_isoclinic.py. CC0. Co-authored."""
import math
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
import helix_isoclinic as X  # noqa: E402

OUT = X.run(0)


class TestQ23(unittest.TestCase):
    def test_arithmetic(self):
        s2, s3, s6 = X.Q23(0, 1), X.Q23(0, 0, 1), X.Q23(0, 0, 0, 1)
        self.assertEqual(s2 * s2, 2)
        self.assertEqual(s3 * s3, 3)
        self.assertEqual(s2 * s3, s6)
        x = X.Q23(F(1, 3), -2, F(5, 7), 4)
        self.assertEqual(x * x.inverse(), 1)
        self.assertAlmostEqual(float(x), 1 / 3 - 2 * 2 ** 0.5 + 5 / 7 * 3 ** 0.5 + 4 * 6 ** 0.5)

    def test_cos15(self):
        for k in range(-5, 30):
            self.assertAlmostEqual(float(X.cos15(k)), math.cos(math.radians(15 * k)), places=12)
        self.assertEqual(X.cos15(6), 0)
        self.assertEqual(X.sin15(2), X.Q23(F(1, 2)))


class TestRun(unittest.TestCase):
    def test_tags(self):
        t = {(r["dimension_logged"], r["simplex_m"], r["n_to_close"]): r for r in OUT["tags"]}
        self.assertEqual(t[(4, 3, 4)]["tag"], "SETTING_OUTSIDE")
        self.assertEqual(t[(4, 3, 4)]["tag_2026_11_16_as_ordered"], "SETTING_IN")
        self.assertEqual(t[(6, 5, 6)]["relation"], "sphere dim < n - 1")
        for n in (6, 10, 24):
            self.assertEqual(t[(6, 4, n)]["tag"], "SETTING_IN")
            self.assertEqual(t[(6, 4, n)]["tag_2026_11_16_as_ordered"], "SETTING_OUTSIDE")
        for r in OUT["tags"]:
            self.assertEqual(r["tag"] == "SETTING_IN",
                             r["points_span_dimension"] == r["simplex_vertices_n"])
            if r["tag"] == "SETTING_IN":
                self.assertTrue(r["eq_5_1"]["holds_on_every_plane_at_1e-12"])
            else:
                self.assertIsNone(r["eq_5_1"])
        self.assertEqual(OUT["summary"]["setting_in"], 9)

    def test_exact_24(self):
        e = OUT["exact_n24"]
        self.assertTrue(e["gaps_1_to_4_exactly_equal"])
        self.assertTrue(e["cos_dihedral_equals_sqrt2_minus_1_over_2"])
        self.assertTrue(e["one_plus_sec_equals_3_plus_2_sqrt2"])
        self.assertTrue(e["eq_5_1_exact_on_both_planes"])
        self.assertAlmostEqual(sum(e["radius_squared_shares"]), 1.0)
        self.assertAlmostEqual(e["radius_squared_shares"][0], 0.7146687612194447)

    def test_7_13(self):
        c = OUT["chain_7_13"]
        self.assertTrue(c["vertex_set_is_30_vertices_of_a_600_cell"])
        self.assertEqual(c["edge_is_chord_number"], 5)
        self.assertAlmostEqual(c["edge_angle_deg"], 108.0)
        self.assertEqual(len(c["chords_of_the_600_cell_deg"]), 8)
        self.assertEqual(c["consecutive_vertices_forming_a_regular_simplex"], 4)
        self.assertIsNone(c["petrie_polygon_of_3_3_5_over_2"])

    def test_isoclinic_powers(self):
        self.assertEqual([p["s"] for p in X.isoclinic_powers(8, (1, 3))], [2, 4, 6])
        self.assertEqual(X.isoclinic_powers(5, (1, 2)), [])
        p = X.isoclinic_powers(18, (1, 5, 8))[0]
        self.assertEqual((p["s"], p["c"], p["sub_orbit_size"], p["n_sub_orbits"]), (6, 6, 3, 6))

    def test_isoclinic_rows(self):
        for r in OUT["isoclinic_logged"]:
            s = r["smallest"]
            if s and s["sub_orbit_size"] > 2:
                self.assertTrue(s["regular_polygons_on_great_circles"])
                self.assertTrue(s["sub_orbits_cover_every_vertex_once"])
                self.assertTrue(s["hopf"]["pairwise_clifford_parallel"])
        sm = OUT["summary"]
        self.assertEqual(sm["logged_with_an_isoclinic_power"], 9)
        self.assertEqual(sm["survivors_without"], [[24, [1, 7, 12]]])
        self.assertEqual(sm["r7_chains"], sm["r7_with_an_isoclinic_power"] + sm["r7_without"])
        for r in OUT["isoclinic_r7"]:
            if r["has_an_isoclinic_power"]:
                self.assertTrue(r["only_minus_identity"])

    def test_not_parallel_is_seen(self):
        # two great circles sharing a point are not Clifford parallel
        r = X.read_sub_orbits(8, (1, 3), (0.5, 0.5), X.isoclinic_powers(8, (1, 3))[0])
        self.assertTrue(r["hopf"]["pairwise_clifford_parallel"])
        a = [[1.0, 0, 0, 0], [0, 1.0, 0, 0]]
        b = [[1.0, 0, 0, 0], [0, 0, 1.0, 0]]
        mm = [[sum(x * y for x, y in zip(a[i], b[j])) for j in range(2)] for i in range(2)]
        self.assertNotAlmostEqual(mm[0][0] ** 2 + mm[1][0] ** 2, mm[0][1] ** 2 + mm[1][1] ** 2)


if __name__ == "__main__":
    unittest.main()
