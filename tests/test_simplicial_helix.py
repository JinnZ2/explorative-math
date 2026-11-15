"""Tests for experiments/simplicial_helix.py. CC0. Co-authored."""
import math
import unittest

import _path  # noqa: F401
import emlib
import simplicial_helix as S


class Equation(unittest.TestCase):
    def test_r3_tetrahedral_helix_holds_at_half_turn_with_n_4(self):
        twist = math.degrees(math.acos(-2.0 / 3.0))
        self.assertTrue(S.residuals(twist / 2.0, 4)["holds"])
        self.assertFalse(S.residuals(twist, 4)["holds"])
        self.assertFalse(S.residuals(twist / 2.0, 3)["holds"])

    def test_bounded_form_agrees_with_the_tan_form(self):
        r = S.residuals(20.0, 4)
        t = math.radians(20.0)
        self.assertAlmostEqual(r["residual_tan_form"] * math.cos(t) * math.cos(4 * t),
                               r["residual_bounded_form"])
        self.assertIsNone(S.residuals(90.0, 5)["residual_tan_form"])       # void, not a number

    def test_logged_run(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = emlib.load_run(emlib.RUNS / "simplicial_helix_seed0.json")
        self.assertIsNone(rec["ai_interest"])
        s = rec["outputs"]["summary"]
        self.assertEqual((s["chains"], s["holding_on_every_plane"], s["survivors"],
                          s["survivors_holding_on_every_plane"]), (17, 0, 4, 0))
        self.assertEqual(s["planes_holding"], s["planes_holding_that_are_theta_90"])
        self.assertTrue(rec["outputs"]["reference_r3"]["half_turn"]["holds"])
        self.assertEqual(sum(1 for c in rec["outputs"]["chains"] if c["survivor_unmatched"]), 4)


if __name__ == "__main__":
    unittest.main()
