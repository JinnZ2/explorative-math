"""Tests for experiments/corkscrew6.py and the multi-plane exact closure. CC0. Co-authored."""
import unittest
from fractions import Fraction

import _path  # noqa: F401
import corkscrew6 as C
import emlib
from exact import corkscrew as XK
from toolkit import primitives as P


class Exact(unittest.TestCase):
    def test_return_is_the_lcm_of_the_plane_periods(self):
        pts = XK.three_pentagons()
        self.assertEqual(len(set(pts)), 15)
        c = XK.closure_multi(pts, (1, 5, 10))
        self.assertEqual(c["plane_return_steps"], [10, 2, 1])
        self.assertEqual(c["point_return_step"], 10)
        self.assertTrue(c["miss_is_exactly_zero_at_return"])
        c = XK.closure_multi(pts, (2, 4, 6))
        self.assertEqual((c["set_return_step"], c["point_return_step"]), (1, 5))


class Pieces(unittest.TestCase):
    def test_period_and_lcm(self):
        self.assertEqual(C.period(36.0, Fraction(1)), 10)
        self.assertEqual(C.period(36.0, Fraction(5, 3)), 6)
        self.assertEqual(C.lcm([10, 20, 6]), 60)

    def test_follow_logs_each_plane(self):
        f = C.follow(P.cross_polytope(6), [36.0, 54.0, 60.0], 200)
        self.assertEqual(f["first_return_step_per_plane"], [10, 20, 6])
        self.assertEqual(f["point_return_step"], 60)
        self.assertLessEqual(f["set_return_step"], 60)

    def test_plane_misses_are_separate(self):
        v = P.hypercube(6)
        cur = C.turned(v, [0.0, 5.0, 0.0], 1)
        pm = C.plane_misses(v, cur, 1.0)
        self.assertEqual(pm[0], 0.0)
        self.assertGreater(pm[1], 0.0)
        self.assertEqual(pm[2], 0.0)

    def test_chains_small(self):
        c = C.chain_part(14)
        by = c["chains_using_all_three_planes"]
        self.assertEqual(by["6"]["n_values"], [7])                  # the 6-simplex itself
        self.assertEqual(by["5"]["n_values"], [6, 7, 12])
        self.assertTrue(set(by["6"]["n_values"]) <= set(by["5"]["n_values"])
                        <= set(by["4"]["n_values"]))


class LoggedRuns(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = emlib.load_run(emlib.RUNS / "corkscrew6_seed0.json")
        cls.o = cls.rec["outputs"]

    def test_observer_slots_null_and_earlier_run_untouched(self):
        self.assertIsNone(self.rec["ai_interest"])
        self.assertIsNone(self.rec["human_interest"])
        old = emlib.load_run(emlib.RUNS / "corkscrew_seed0.json")
        self.assertEqual(old["outputs"]["tetrahedra_chain"]["r4"]
                         ["n_values_with_a_closed_chain"], [4, 5, 8, 30])

    def test_observed_sweep(self):
        # OBSERVED on the committed run, pinned so a change shows up
        self.assertEqual([s["shape"] for s in self.o["sweep"]],
                         ["6-simplex", "6-orthoplex", "6-cube"])
        for s in self.o["sweep"]:
            self.assertEqual(s["by_kind"]["all rational"],
                             {"n": 6, "point_return": 6, "set_return": 6})
            self.assertEqual(s["by_kind"]["one irrational"]["point_return"], 0)
            self.assertEqual(s["by_kind"]["one irrational"]["set_return"], 0)
            self.assertEqual(s["by_kind"]["two irrational"]["point_return"], 0)
            for r in s["rows"]:
                if r["kind"] == "all rational":
                    self.assertTrue(r["point_return_matches_lcm"], r["rates"])
                    self.assertLessEqual(r["set_return_step"], r["point_return_step"])
                else:
                    self.assertIsNone(r["predicted_point_return_lcm"])

    def test_observed_miss_and_partial(self):
        m = self.o["miss"]
        self.assertTrue(m["all_single_confined"])
        self.assertTrue(m["all_double_same_as_alone"])
        for c in self.o["partial_closure"]["cases"]:
            self.assertTrue(c["rational_planes_back_at_every_lap"])
            self.assertFalse(c["third_plane_ever_back_at_a_lap"])
            self.assertIsNone(c["full_point_return_step"])

    def test_observed_chains(self):
        by = self.o["chains"]["chains_using_all_three_planes"]
        self.assertEqual(self.o["bound"]["n_max_chain"], 120)
        self.assertEqual(by["6"]["n_values"], [7])
        self.assertEqual(by["5"]["n_values"], [6, 7, 12])
        self.assertEqual(by["4"]["n_values"], [6, 7, 10, 12, 18, 24, 28])
        self.assertGreater(by["3"]["n_found"], 1000)

    def test_small_run(self):
        small = emlib.load_run(emlib.RUNS / "corkscrew6_small_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(self.o))


if __name__ == "__main__":
    unittest.main()
