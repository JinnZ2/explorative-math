"""Tests for bound logging, orbit_exact and symmetry_resolution. CC0. Co-authored."""
import json
import unittest

import _path  # noqa: F401
import fibonacci_chain as F
import orbit_exact as O
import symmetry_resolution as S
from toolkit import generator as G
from toolkit import primitives as P
from toolkit import probes as Q


class BoundHelpers(unittest.TestCase):
    def test_probe_return_values_unchanged(self):
        self.assertEqual(Q.rotational_symmetry(P.regular_polygon(8)), [2, 4, 8])
        self.assertEqual(Q.smallest_period("abab"), 2)

    def test_bounds(self):
        sq = P.regular_polygon(4)
        b = Q.rotational_symmetry_bound(sq, 24, 1e-9)
        self.assertEqual((b["n_points"], b["k_tested"], b["rel_tol"], b["about"]),
                         (4, [2, 24], 1e-9, "centroid"))
        self.assertAlmostEqual(b["abs_tol"], 1e-9)
        m = Q.mirror_symmetry_bound(sq)
        self.assertEqual((m["n_axes"], m["axis_step_deg"]), (360, 0.5))
        self.assertEqual(Q.period_bound("a" * 20), {"n": 20, "prefix_lengths": [8, 16]})
        self.assertIsNone(Q.rotational_symmetry_bound([])["abs_tol"])

    def test_generator_logs_bounds_for_each_pipeline(self):
        want = {"shape": {"rotational_symmetry_bound", "mirror_symmetry_bound",
                          "point_rounding_decimals"},
                "cut": {"n_points", "gap_cluster_sep", "period_bound"},
                "orbit": {"period_bound", "arithmetic"}}
        seen = set()
        for seed in range(10):
            spec, out, *_ = G.run(seed)
            self.assertTrue(want[spec["pipeline"]] <= set(out), spec["pipeline"])
            seen.add(spec["pipeline"])
        self.assertEqual(seen, set(want))
        _, out, *_ = G.run(0, "cut_weighted")
        self.assertIn("gap_cluster_sep", out)
        self.assertIn("period_bound", out["levels"][0])

    def test_fibonacci_logs_bounds(self):
        g = F.run()["golden"]
        self.assertEqual(g["gap_cluster_sep"], F.SEP)
        self.assertEqual(g["n_points"], g["n_gaps"] + 1)
        self.assertEqual(g["period_bound"]["n"], len(g["word"]))


class MirrorFromPairs(unittest.TestCase):
    def test_agrees_with_grid_on_grid_axes(self):
        sq = P.regular_polygon(4)
        self.assertEqual(Q.mirror_symmetry_from_pairs(sq), Q.mirror_symmetry(sq))

    def test_finds_axes_off_the_grid(self):
        pent = P.regular_polygon(5, phase=0.1234)
        self.assertEqual(Q.mirror_symmetry(pent), [])
        self.assertEqual(len(Q.mirror_symmetry_from_pairs(pent)), 5)

    def test_no_mirror(self):
        self.assertEqual(Q.mirror_symmetry_from_pairs([[0, 0], [1, 0], [0, 2], [3, 3]]), [])


class CompareBounds(unittest.TestCase):
    def test_compare_records_carry_bounds_outside_the_tally(self):
        from compare import cross as X
        for fn in (X.cmp_angle_measure, X.cmp_phi_scale_vs_phi_scale_rotate):
            for r in fn():
                b = r["outputs"]["bounds"]["rotational_symmetry_k"]
                for side in b.values():
                    self.assertEqual(side["k_tested"], [2, 12])
                self.assertEqual(sum(r["outputs"]["tally"][v] for v in X.VERDICTS),
                                 len(r["outputs"]["fields"]))


class OrbitExact(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = O.run(extended_steps=1024, prefixes=(512, 1024))

    def test_float_word_is_the_logged_word(self):
        self.assertTrue(self.out["float"]["word_equals_logged"])

    def test_observed_exact_equals_float(self):
        # OBSERVED on this run, pinned so a change shows up
        e = self.out["exact"]
        self.assertTrue(e["word_equals_float"])
        self.assertIsNone(e["first_index_where_words_differ"])
        self.assertEqual(e["period_profile"]["512"], 184)
        self.assertTrue(self.out["float_precision_check"]["persists"])

    def test_observed_extended_profile(self):
        self.assertEqual(self.out["extended"]["period_profile"], {"512": 184, "1024": 689})
        self.assertEqual(self.out["extended"]["period_bound"]["n"], 1025)

    def test_orbit_period_is_alpha_denominator(self):
        p = self.out["exact_orbit_period"]
        self.assertTrue(p["alpha_denominator_is_power_of_two"])
        self.assertEqual(p["orbit_period"], 2 ** p["alpha_denominator_log2"])
        self.assertEqual(p["alpha_denominator_log2"], 52)

    def test_first_difference(self):
        self.assertEqual(O.first_difference("0101", "0111"), 2)
        self.assertIsNone(O.first_difference("01", "01"))


class SymmetryResolution(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = S.run()

    def test_observed_known_pair(self):
        # OBSERVED on this run, pinned so a change shows up
        self.assertEqual(self.out["known_separate_by_bound"],
                         [[12, 1e-6, False], [12, 1e-9, False],
                          [24, 1e-6, True], [24, 1e-9, True]])

    def test_every_row_logs_its_bound(self):
        for r in self.out["known"]["rows"]:
            self.assertEqual(r["a_bound"]["k_tested"], [2, r["k_max"]])
            self.assertEqual(r["a_bound"]["rel_tol"], r["rel_tol"])

    def test_observed_unions_do_not_separate(self):
        self.assertEqual(len(self.out["unions"]), 8)
        self.assertEqual([n for _, _, n, _ in self.out["union_frames_separating_by_bound"]],
                         [0, 0, 0, 0])
        for f in self.out["unions"]:
            self.assertLessEqual(max(f["n_points"]), 24)

    def test_observed_mirror_on_unions(self):
        # OBSERVED on this run, pinned so a change shows up
        self.assertEqual([c for _, _, c, _ in self.out["mirror_frames_separating_by_bound"]],
                         [2, 2, 2, 2])
        self.assertEqual([c for _, c, _ in self.out["mirror_frames_separating_from_pairs"]],
                         [2, 2])
        self.assertTrue(self.out["prediction"]["rotation_lists_equal_on_all_frames_and_bounds"])
        for f in self.out["unions_mirror"]:
            for r in f["rows"]:
                self.assertEqual(r["a_bound"]["n_axes"], r["n_axes"])
            # grid and no-grid readings agree on these frames
            self.assertEqual(len(f["rows"][0]["a"]), len(f["from_pairs"][0]["a"]))
            self.assertEqual(len(f["rows"][0]["b"]), len(f["from_pairs"][0]["b"]))

    def test_json_serialisable(self):
        json.dumps(self.out)


if __name__ == "__main__":
    unittest.main()
