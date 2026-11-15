"""Tests for experiments/coxeter_exponents.py. CC0. Co-authored."""
import unittest

import _path  # noqa: F401
import coxeter_exponents as X
import emlib


class Computed(unittest.TestCase):
    def test_exponents_have_the_right_count_and_sum(self):
        # for every group: as many exponents as the rank, summing to rank * h / 2,
        # and symmetric under m -> h - m
        for g in X.groups():
            self.assertTrue(g["n_exponents_equals_rank"], g["group"])
            self.assertTrue(g["exponents_sum_to_rank_h_over_2"], g["group"])
            self.assertEqual(sorted(g["h"] - m for m in g["exponents"]), g["exponents"])
            self.assertEqual(g["exponents"][0], 1)

    def test_a_few_computed_values(self):
        got = {g["group"]: (g["h"], g["exponents"]) for g in X.groups()}
        self.assertEqual(got["A4"], (5, [1, 2, 3, 4]))
        self.assertEqual(got["B4"], (8, [1, 3, 5, 7]))
        self.assertEqual(got["D4"], (6, [1, 3, 3, 5]))
        self.assertEqual(got["H4"], (30, [1, 11, 19, 29]))
        self.assertEqual(got["E8"][0], 30)

    def test_logged_run(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = emlib.load_run(emlib.RUNS / "coxeter_exponents_seed0.json")
        self.assertIsNone(rec["ai_interest"])
        s = rec["outputs"]["summary"]
        self.assertEqual((s["chains"], s["match"], len(s["no_match"])), (17, 13, 4))
        rows = {(r["dimension"], r["simplex_m"], r["n"], tuple(r["turns"])): r
                for r in rec["outputs"]["chains"]}
        self.assertEqual(rows[(4, 3, 30, (1, 11))]["matches_exactly"], ["H4"])
        self.assertEqual(rows[(4, 3, 30, (7, 13))]["matches"], ["E8"])
        self.assertEqual(rows[(6, 4, 18, (1, 5, 8))]["verdict"], "no match")


if __name__ == "__main__":
    unittest.main()
