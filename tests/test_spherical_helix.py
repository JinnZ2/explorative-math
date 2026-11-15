"""Tests for experiments/spherical_helix.py. CC0. Co-authored."""
import math
import unittest

import _path  # noqa: F401
import emlib
import spherical_helix as S


class Pieces(unittest.TestCase):
    def test_flat_turns(self):
        self.assertAlmostEqual(S.flat_helix_turns(4)["turns_deg"][0],
                               math.degrees(math.acos(-2 / 3)), places=6)
        t5 = S.flat_helix_turns(5)["turns_deg"]
        self.assertAlmostEqual(t5[0], math.degrees(math.acos(-1 / 4)), places=6)
        self.assertAlmostEqual(t5[1], 180.0, places=6)
        self.assertEqual(S.flat_helix_turns(6)["coefficients"], [5, 8, 9, 8, 5])

    def test_five_cell_itself_is_one_window(self):
        r = S.read("5-cell", 5, (1, 2), (0.5, 0.5), 4)
        self.assertTrue(r["every_window_is_a_regular_simplex"])
        self.assertTrue(r["step_is_one_isometry"])
        self.assertEqual(r["points_span_dimension"], 4)

    def test_logged_run(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = emlib.load_run(emlib.RUNS / "spherical_helix_seed0.json")
        self.assertIsNone(rec["ai_interest"])
        o = rec["outputs"]
        sm = o["summary"]
        self.assertEqual((sm["survivors"], sm["every_window_regular"],
                          sm["step_is_one_isometry"], sm["ring_every_window_regular"],
                          sm["ring_step_is_one_isometry"]), (4, 4, 4, True, True))
        ring = o["ring_reference"]
        self.assertAlmostEqual(ring["hinge_angle_deg"]["min"], 164.478, places=2)
        self.assertEqual(ring["sphere_dimension_minus_simplex_dimension"], 0)
        got = {x["n_to_close"]: [] for x in o["survivors"]}
        for x in o["survivors"]:
            got[x["n_to_close"]].append((x["points_span_dimension"],
                                         x["sphere_dimension_minus_simplex_dimension"]))
            self.assertEqual(x["sphere_dimension"], x["points_span_dimension"] - 1)
            self.assertFalse(x["next_vertex_is_at_the_edge_distance"])
        # the N = 24 chain has one plane turned by 180 degrees: its points span R^5
        self.assertEqual(got, {18: [(6, 1)], 24: [(5, 0)], 28: [(6, 1), (6, 1)]})
        for x in o["survivors"]:
            self.assertTrue(x["hinge_angle_deg"]["same_at_every_facet"])


if __name__ == "__main__":
    unittest.main()
