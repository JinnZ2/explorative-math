"""Tests for exact/corkscrew.py and experiments/corkscrew.py. CC0. Co-authored."""
import math
import unittest
from fractions import Fraction

import _path  # noqa: F401
import corkscrew as C
import emlib
from exact import ZERO
from exact import corkscrew as XK
from toolkit import primitives as P


class Exact(unittest.TestCase):
    def test_five_cell_is_a_regular_simplex_with_no_tolerance(self):
        f = XK.five_cell()
        self.assertEqual(len(set(f)), 5)
        self.assertTrue(XK.is_regular_simplex(f))
        self.assertFalse(XK.is_regular_simplex(XK.duoprism_5_5()))

    def test_return_steps(self):
        f = XK.five_cell()
        want = {(2, 4): (1, 5), (1, 2): (2, 10), (2, 2): (5, 5), (5, 5): (2, 2)}
        for (p, q), (s, k) in want.items():
            c = XK.closure(f, p, q)
            self.assertEqual((c["set_return_step"], c["point_return_step"]), (s, k), (p, q))
            self.assertTrue(c["miss_is_exactly_zero_at_return"])
            self.assertEqual(len(c["miss_squared_of_first_point_per_step"]), k)

    def test_a_full_turn_is_the_identity_exactly(self):
        f = XK.five_cell()
        self.assertEqual(XK.turn(f, 10, 10), f)
        self.assertEqual(XK.dist2(f[0], f[0]), ZERO)


class Float(unittest.TestCase):
    def test_predicted_return(self):
        self.assertEqual(C.predicted_return(36.0, Fraction(1)), 10)
        self.assertEqual(C.predicted_return(36.0, Fraction(3, 2)), 20)
        self.assertEqual(C.predicted_return(36.0, Fraction(1618, 1000)), 5000)

    def test_follow_returns_at_the_predicted_step(self):
        f = C.follow(P.hypercube(4), 36.0, 54.0, 100)
        self.assertEqual(f["point_return_step"], 20)
        self.assertLess(f["smallest_miss"], 1e-12)

    def test_set_return_can_come_before_point_return(self):
        f = C.follow(P.cross_polytope(4), 90.0, 90.0, 20)
        self.assertEqual((f["set_return_step"], f["point_return_step"]), (1, 4))

    def test_fifth_axis_is_fixed(self):
        v = P.hypercube(5)
        cur = C.turned(v, 33.0, 71.0, 7)
        self.assertEqual([p[4] for p in cur], [p[4] for p in v])

    def test_chain_search_finds_only_whole_chains(self):
        c = C.chain_part(32)
        self.assertEqual(c["r4"]["n_values_with_a_closed_chain"], [4, 5, 8, 30])
        for f in c["r4"]["closed_chains"]:
            self.assertLess(f["largest_relative_edge_difference"], 1e-9)
        self.assertAlmostEqual(c["r3"]["twist_deg"], math.degrees(math.acos(-2 / 3)))
        for d in c["r3"]["edges_d1_d2_d3"]:
            self.assertAlmostEqual(d, 1.0, places=12)


class LoggedRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = emlib.load_run(emlib.RUNS / "corkscrew_seed0.json")
        cls.o = cls.rec["outputs"]

    def test_observer_slots_null(self):
        self.assertIsNone(self.rec["ai_interest"])
        self.assertIsNone(self.rec["human_interest"])

    def test_observed_sweep(self):
        # OBSERVED on the committed run, pinned so a change shows up
        for s in self.o["sweep"]:
            self.assertEqual((s["rational_closing"], s["rational_ratios"]), (9, 10))
            self.assertEqual(s["rational_not_closing_within_steps"], ["3/2 + 1/10000"])
            self.assertEqual((s["irrational_closing"], s["irrational_ratios"]), (0, 4))
            for r in s["rows"]:
                if r["rational"] and r["predicted_within_the_steps_run"]:
                    self.assertTrue(r["return_matches_prediction"], r["ratio"])
                if r["rational"] and not r["predicted_within_the_steps_run"]:
                    self.assertIsNone(r["return_matches_prediction"])     # void, not False

    def test_observed_error_test(self):
        for e in self.o["error_test"]["perturbations"]:
            self.assertAlmostEqual(e["log_log_slope_of_miss_against_eps"], 1.0, places=4)
            self.assertEqual(e["plane_carrying_the_miss"], [e["expected_plane"]])
            lo, hi = e["miss_over_linear_estimate_range"]
            self.assertTrue(0.999 < lo <= hi <= 1.001)
        self.assertLess(self.o["error_test"]["nominal_miss_at_return"], 1e-12)

    def test_observed_chain(self):
        r4 = self.o["tetrahedra_chain"]["r4"]
        self.assertTrue(r4["closes_at_30"])
        self.assertEqual(r4["n_values_with_a_closed_chain"], [4, 5, 8, 30])
        self.assertEqual(r4["turn_pairs_by_n"]["30"], [[1, 11], [7, 13]])
        self.assertAlmostEqual(min(r4["edge_over_radius_at_30"]), r4["one_over_phi"], places=9)
        r3 = self.o["tetrahedra_chain"]["r3"]
        self.assertFalse(r3["returns_to_its_start"])
        self.assertFalse(r3["twist_back_to_zero_within_2000_steps"])

    def test_small_run(self):
        small = emlib.load_run(emlib.RUNS / "corkscrew_small_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(self.o))


if __name__ == "__main__":
    unittest.main()
