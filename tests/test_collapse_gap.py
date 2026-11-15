"""Tests for toolkit/complement.py, exact/complement.py and collapse_gap. CC0. Co-authored."""
import math
import unittest
from fractions import Fraction

import _path  # noqa: F401
import collapse_gap as G
import emlib
from exact import Cyc5, ZETA
from exact import complement as XC
from exact import inflation as XI
from toolkit import complement as CP
from toolkit import primitives as P


class Plane(unittest.TestCase):
    def test_pattern_is_centred_scaled_and_turned(self):
        pat = CP.turned_to_zero(CP.centred_unit([[3.0, 3.0], [5.0, 3.0], [4.0, 6.0]]))
        self.assertAlmostEqual(max(math.hypot(*p) for p in pat), 1.0)
        far = max(pat, key=lambda p: math.hypot(*p))
        self.assertAlmostEqual(math.atan2(far[1], far[0]), 0.0, places=12)

    def test_field_reads_the_polygons_symmetry(self):
        for k in (3, 5, 8):
            r = CP.read_plane(P.regular_polygon(k, phase=0.37), 3, 840)
            self.assertEqual(r["field_rotational_orders"],
                             [d for d in range(2, 49) if k % d == 0 and 840 % d == 0])
            self.assertEqual(len(r["field_mirror_axes_deg"]), k)
            self.assertEqual(len(r["gaps_deg"]), k)
            self.assertAlmostEqual(r["centre_depth"], 1.0)
            self.assertTrue(r["centre_is_a_void"])

    def test_no_symmetry_read_from_a_lopsided_set(self):
        r = CP.read_plane([[1.0, 0.0], [-0.4, 0.9], [0.2, -0.3], [-0.8, -0.1]], 3, 840)
        self.assertEqual(r["field_rotational_orders"], [])
        self.assertEqual(r["field_mirror_axes_deg"], [])

    def test_orders_the_grid_cannot_test_are_named(self):
        r = CP.read_plane(P.regular_polygon(9), 2, 840)       # 9 does not divide 840
        self.assertNotIn(9, r["bound"]["orders_tested"])
        self.assertEqual(r["field_rotational_orders"], [3])

    def test_field_difference_is_zero_for_the_same_pattern_turned(self):
        a = CP.read_plane(P.regular_polygon(5), 3, 840)
        radii, n = CP.polar_grid(3, 840)
        fa = CP.distance_field(a["pattern"], radii, n)
        fb = CP.distance_field(CP.turned_to_zero(CP.centred_unit(
            P.regular_polygon(5, phase=1.0))), radii, n)
        self.assertLess(CP.field_difference(fa, fb, n)["rms_min_over_turns"], 1e-12)
        fc = CP.distance_field(CP.turned_to_zero(CP.centred_unit(P.regular_polygon(6))), radii, n)
        self.assertGreater(CP.field_difference(fa, fc, n)["rms_min_over_turns"], 1e-2)

    def test_gaps(self):
        self.assertEqual(len(CP.angular_gaps([[1, 0], [0, 1], [-1, 0], [0, -1], [0, 0]])), 4)
        self.assertEqual(CP.angular_gaps([[0.0, 0.0]]), [])


class Sphere(unittest.TestCase):
    def test_voids_of_the_16_cell_sit_on_tesseract_directions(self):
        v = CP.sphere_voids(P.cross_polytope(4), seed=0, starts=60, steps=300)
        self.assertAlmostEqual(v["deepest_depth_deg"], 60.0, delta=0.1)
        m = CP.match_directions([x["direction"] for x in v["voids"]], P.hypercube(4))
        self.assertLess(m["max_angle_deg"], 1.0)
        self.assertGreaterEqual(m["candidates_hit"], 12)        # found is a lower bound
        off = CP.match_directions([x["direction"] for x in v["voids"]], P.cross_polytope(4))
        self.assertGreater(off["max_angle_deg"], 30.0)

    def test_match_of_nothing_is_void(self):
        self.assertIsNone(CP.match_directions([], P.hypercube(4))["max_angle_deg"])


class Exact(unittest.TestCase):
    def test_sample_set_is_closed_under_36_degree_turns_and_the_mirror(self):
        s = set(XC.samples())
        self.assertTrue(all(z * XC.Z10 in s and z.conjugate() in s for z in s))

    def test_pentagon_decagon_and_control(self):
        pent = XI.pentagon()
        r = XC.read(pent)
        self.assertEqual(r["field_unchanged_by_turn_deg"], {"36": False, "72": True, "144": True})
        self.assertTrue(r["field_unchanged_by_mirror_in_real_axis"])
        dec = pent + [z * XC.Z10 for z in pent]
        self.assertTrue(all(XC.read(dec)["field_unchanged_by_turn_deg"].values()))
        odd = [Cyc5(1, 0), ZETA, ZETA * ZETA * 2 + Cyc5(Fraction(1, 3), 0)]
        ro = XC.read(odd)
        self.assertFalse(any(ro["field_unchanged_by_turn_deg"].values()))
        self.assertFalse(ro["field_unchanged_by_mirror_in_real_axis"])


class Run(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = G.run(seed=0, n_rings=3, n_theta=840, sphere_starts=30, sphere_steps=200,
                        baselines=("5-cell", "16-cell"))

    def test_checks_and_controls(self):
        for c in self.out["plane"]["checks"]:
            self.assertTrue(c["recovers_rotational_orders"], c["shape"])
            self.assertTrue(c["recovers_mirror_count"], c["shape"])
            self.assertEqual(c["any_symmetry_read_from_the_complement"], not c["control"])

    def test_same_shape_by_another_route_is_the_nearest_probe(self):
        want = {"5-cell": "5-gon", "16-cell": "8-gon"}
        for pb in self.out["plane"]["per_baseline"]:
            self.assertEqual(pb["probes_by_field_difference"][0], want[pb["baseline"]])
            self.assertLess(pb["deltas"][0]["field_rms_min_over_turns"], 1e-9)
            self.assertEqual(len(pb["deltas"]), len(self.out["plane"]["readings"]) - 1)

    def test_full_field_only_for_baselines_and_bounds_logged(self):
        r = self.out["plane"]["readings"]
        self.assertIn("field", r["5-cell"])
        self.assertNotIn("field", r["7-gon"])
        self.assertIn("orders_tested", r["7-gon"]["bound"])


class LoggedRun(unittest.TestCase):
    def test_observed(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = emlib.load_run(emlib.RUNS / "collapse_gap_seed0.json")
        self.assertIsNone(rec["ai_interest"])
        o = rec["outputs"]
        self.assertEqual(len(o["plane"]["checks"]), 12)
        self.assertTrue(all(c["recovers_rotational_orders"] and c["recovers_mirror_count"]
                            for c in o["plane"]["checks"]))
        cover = {r["shape"]: r["probes_fully_covered"] for r in o["sphere"]["shapes"]}
        self.assertEqual(cover["16-cell"], ["tesseract"])
        self.assertEqual(cover["tesseract"], ["16-cell"])
        self.assertEqual(cover["5-cell"], ["5-cell negated"])
        self.assertEqual(cover["24-cell"], ["16-cell + tesseract directions"])
        self.assertEqual(cover["random_8"], [])
        small = emlib.load_run(emlib.RUNS / "collapse_gap_small_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(o))


if __name__ == "__main__":
    unittest.main()
