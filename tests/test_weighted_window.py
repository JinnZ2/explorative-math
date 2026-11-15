"""Tests for experiments/weighted_window.py and the weight option. CC0. Co-authored."""
import math
import unittest

import _path  # noqa: F401
import weighted_window as E
from toolkit import composers as C
from toolkit import primitives as P


class WeightOption(unittest.TestCase):
    def setUp(self):
        self.par, self.perp = C.line_basis_2d(1 / P.PHI)
        self.lat = P.lattice_zn(2, 10)

    def test_hard_window_result_unchanged_without_weight(self):
        w = C.strip_window_1d(self.perp[0])
        self.assertEqual(C.cut_and_project(self.lat, self.par, self.perp, w),
                         C.cut_and_project(self.lat, self.par, self.perp, window=w))

    def test_weight_keeps_every_point(self):
        recs = C.cut_and_project(self.lat, self.par, self.perp,
                                 weight=C.gaussian_weight(0.5))
        self.assertEqual(len(recs), len(self.lat))
        self.assertEqual(set(recs[0]), {"point", "perp", "weight"})

    def test_gaussian_values(self):
        g = C.gaussian_weight(2.0)
        self.assertEqual(g([0.0]), 1.0)
        self.assertAlmostEqual(g([2.0]), math.exp(-0.5))
        self.assertAlmostEqual(C.gaussian_weight(1.0, center=[1.0])([1.0]), 1.0)
        with self.assertRaises(ValueError):
            C.gaussian_weight(0.0)

    def test_window_and_weight_together(self):
        w = C.strip_window_1d(self.perp[0])
        both = C.cut_and_project(self.lat, self.par, self.perp, w, C.gaussian_weight(0.5))
        hard = C.cut_and_project(self.lat, self.par, self.perp, w)
        self.assertEqual([r["point"] for r in both], hard)


class WeightedWindow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = E.run()

    def test_both_windows_logged(self):
        self.assertIn("hard", self.out)
        self.assertIn("gaussian", self.out)

    def test_hard_two_lengths_ratio_phi(self):
        ls = [c["length"] for c in self.out["hard"]["distinct_gap_lengths"]]
        self.assertEqual(len(ls), 2)
        self.assertAlmostEqual(ls[1] / ls[0], P.PHI, places=9)

    def test_weights_in_unit_interval_and_above_floor(self):
        for _, _, w in self.out["gaussian"]["points_x_q_weight"]:
            self.assertTrue(E.WEIGHT_FLOOR <= w <= 1.0)

    def test_level_is_centred_hard_window(self):
        for lv in self.out["gaussian"]["levels"]:
            n = sum(1 for _, q, _ in self.out["gaussian"]["points_x_q_weight"]
                    if abs(q) <= lv["half_width_in_q"] + 1e-12)
            self.assertEqual(n, lv["n_points"])

    def test_observed_gap_length_counts(self):
        # OBSERVED on this run, pinned so a change shows up
        got = [len(lv["distinct_gap_lengths"]) for lv in self.out["gaussian"]["levels"]]
        self.assertEqual(got, [3, 2, 3, 3, 3])


if __name__ == "__main__":
    unittest.main()
