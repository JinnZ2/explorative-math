"""Tests for experiments/e8_orbits.py. CC0. Co-authored."""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
import e8_orbits as X  # noqa: E402

OUT = X.run(0)


class TestE8Orbits(unittest.TestCase):
    def test_eight_orbits_of_thirty(self):
        self.assertEqual(len(OUT["orbits"]), 8)
        for r in OUT["orbits"]:
            self.assertEqual(r["n_points"], 30)
            self.assertTrue(r["step_cosine_same_all_round"])
            self.assertTrue(r["fold_space"]["opposite_point_is_15_steps_on"])

    def test_step_pairs(self):
        got = sorted(tuple(r["steps_deg"]) for r in OUT["orbits"])
        want = sorted([(36.0, 108.0)] * 2 + [(72.0, 144.0)] * 2
                      + [(90.0, 90.0)] * 2 + [(120.0, 120.0)] * 2)
        self.assertEqual(got, want)
        for r in OUT["orbits"]:
            self.assertTrue(r["fold_space"]["step_angle_is_a_listed_exact_angle"])

    def test_petrie_fold_space_only(self):
        s = OUT["summary"]
        self.assertEqual(len(s["petrie_in_fold_space"]), 2)
        self.assertEqual(s["petrie_in_conjugate_space"], [])
        for r in OUT["orbits"]:
            self.assertEqual(r["fold_space"]["is_a_petrie_polygon_of_a_600_cell"],
                             r["steps_deg"][0] == 36.0)

    def test_windows(self):
        for r in OUT["orbits"]:
            for sp in ("fold_space", "conjugate_space"):
                w = r[sp]["largest_regular_simplex_window"]
                self.assertEqual(w, 4 if r["steps_deg"][0] == 36.0 else None)
                self.assertFalse(r[sp]["four_consecutive_are_a_plane_square"])
                self.assertFalse(r[sp]["eight_consecutive_are_a_16_cell"])

    def test_pairing(self):
        p = OUT["pairing"]
        self.assertEqual(len(p["phi_times_pairs"]), 4)
        self.assertEqual(len(p["pairs_summing_to_180"]), 2)
        for a, b in p["phi_times_pairs"]:
            self.assertEqual(OUT["orbits"][a]["steps_deg"], OUT["orbits"][b]["steps_deg"])

    def test_read_on_known_cosines(self):
        from fractions import Fraction
        from exact import ZPhi, ZERO
        q = ZPhi(Fraction(-1, 4))
        r = X.read([q] * 4 + [ZERO] * 11)
        self.assertEqual(r["largest_regular_simplex_window"], 5)
        self.assertTrue(r["five_consecutive_are_a_regular_5_cell"])
        self.assertFalse(r["is_a_petrie_polygon_of_a_600_cell"])


if __name__ == "__main__":
    unittest.main()
