"""Tests for field/hexsheet.py and field_response. CC0. Co-authored."""
import math
import unittest
from fractions import Fraction

import _path  # noqa: F401
import field_response as F
from field import hexsheet as H


class Sheet(unittest.TestCase):
    def test_cell_count_and_exact_reference(self):
        for r in (1, 2, 5):
            self.assertEqual(len(H.cells(r)), 3 * r * (r + 1) + 1)
        x, yu = H.ideal_exact((3, -1))
        self.assertEqual((x, yu), (Fraction(5, 2), Fraction(-1)))
        for o in H.ideal_offsets():
            self.assertAlmostEqual(math.hypot(*o), 1.0, places=12)

    def test_gradient_field_has_zero_mean(self):
        sheet = H.cells(4)
        f = H.gradient_field(sheet, 0.5)
        self.assertAlmostEqual(sum(f.values()) / len(f), 0.0, places=12)

    def test_pearson_and_slope(self):
        self.assertAlmostEqual(H.pearson([1, 2, 3], [2, 4, 6]), 1.0)
        self.assertAlmostEqual(H.slope([1, 2, 3], [2, 4, 6]), 2.0)
        self.assertIsNone(H.pearson([1, 1, 1], [2, 4, 6]))
        self.assertIsNone(H.slope([1, 1, 1], [2, 4, 6]))


class Response(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = F.run(radius=3, sweeps=800)

    def test_every_cell_logs_field_value_and_response(self):
        for key in ("gradient", "control_uniform"):
            cells = self.out[key]["cells"]
            self.assertEqual(len(cells), self.out["bound"]["n_cells"])
            for c in cells:
                for k in ("field_value", "size_response", "shape_response",
                          "offset_response", "interior", "ideal_exact"):
                    self.assertIn(k, c)
                self.assertEqual(len(c["offset_response"]), 6)
                self.assertEqual(c["interior"], all(o is not None for o in c["offset_response"]))

    def test_deviation_is_named_response_not_error(self):
        import json
        text = json.dumps(self.out)
        self.assertNotIn("error", text.lower().replace("not error", ""))
        self.assertIn("size_response", text)

    def test_control_uniform_field_gives_no_response(self):
        p = self.out["control_uniform"]["probe"]
        self.assertLess(p["max_abs_size_response_all_cells"], 1e-12)
        self.assertLess(p["max_shape_response"], 1e-12)
        self.assertIsNone(p["correlation_size_response_vs_field"])      # void, not 0

    def test_gradient_response_follows_the_field(self):
        p = self.out["gradient"]["probe"]
        self.assertGreater(p["correlation_size_response_vs_field"], 0.99)
        self.assertAlmostEqual(p["slope_size_response_vs_field"], F.ALPHA, delta=F.ALPHA * 0.05)

    def test_mean_cell_is_nearer_the_ideal_than_single_cells(self):
        m = self.out["gradient"]["probe"]["mean_cell"]
        self.assertLess(m["max_abs_mean_offset_response"],
                        m["largest_single_cell_offset_response"])
        for x in m["mean_opposite_pair_lengths"]:
            self.assertAlmostEqual(x, 1.0, delta=1e-3)

    def test_relaxation_bound_is_logged(self):
        b = self.out["gradient"]["relaxation_bound"]
        self.assertEqual(b["sweeps"], 800)
        self.assertIn("max_net_pull_left", b)

    def test_logged_runs(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = F.emlib.load_run(F.emlib.RUNS / "field_response_seed0.json")
        o = rec["outputs"]
        self.assertTrue(o["constructed"])
        self.assertEqual(o["bound"]["n_cells"], 127)
        self.assertGreater(o["gradient"]["probe"]["correlation_size_response_vs_field"], 0.9999)
        self.assertLess(o["control_uniform"]["probe"]["max_abs_size_response_all_cells"], 1e-12)
        self.assertIsNone(rec["human_interest"])
        small = F.emlib.load_run(F.emlib.RUNS / "field_response_small_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(o))


if __name__ == "__main__":
    unittest.main()
