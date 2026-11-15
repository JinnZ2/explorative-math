"""Tests for the pairing re-read option, pairing_reread and the tour degeneracy tags.
CC0. Co-authored."""
import unittest

import _path  # noqa: F401
import constrained_generate as K
import emlib
import pairing_reread as R


class PairingType(unittest.TestCase):
    def test_points_on_axis(self):
        self.assertEqual(K.points_on_axis([0, 3, 4, 1, 2]), 1)
        self.assertEqual(K.points_on_axis([0, 1, 2, 4, 3]), 3)
        self.assertEqual(K.points_on_axis([1, 0]), 0)

    def test_reread_option_logs_the_pairing_and_holds_the_anchor(self):
        out = R.run(seed=0, starts=1, limits=(40,), anchor_max_evals=800)
        for s in out["shapes"]:
            for label in R.RULES:
                g = s["starts"][0]["rules"][label]["40"]
                self.assertEqual(g["pairing"], "re-read at every pull-back")
                self.assertIsNotNone(g["pull_backs_with_a_different_pairing"])
                self.assertTrue(g["anchor_held"])
                self.assertEqual(g["anchor_points_on_axis"], K.points_on_axis(g["anchor_pairing"]))

    def test_kept_option_is_unchanged_and_says_so(self):
        out = K.run(seed=0, starts=1, max_trials=30, anchor_max_evals=800,
                    bound_check_trials=None)
        g = out["shapes"][0]["starts"][0]["rules"]["h"]
        self.assertEqual(g["pairing"], "kept from the anchor")
        self.assertIsNone(g["pull_backs_with_a_different_pairing"])


class LoggedRuns(unittest.TestCase):
    def test_observed_reread_counts(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = emlib.load_run(emlib.RUNS / "pairing_reread_seed0.json")
        self.assertIsNone(rec["ai_interest"])
        self.assertIsNone(rec["human_interest"])
        by = {s["shape"]: s["compare"] for s in rec["outputs"]["shapes"]}
        self.assertEqual(by["5-cell"]["anchor_points_on_axis"], [1, 1, 1, 1, 3, 3, 3, 3])
        self.assertEqual(by["16-cell"]["anchor_points_on_axis"], [6, 6, 2, 2, 0, 2, 0, 2])
        got = {n: {lim: c["reread"][lim]["h"]["starts_reaching_h"] for lim in ("400", "2500")}
               for n, c in by.items()}
        self.assertEqual(got, {"5-cell": {"400": 0, "2500": 2}, "16-cell": {"400": 1, "2500": 2}})
        for c in by.values():
            self.assertTrue(c["starts_equal_constrained_generate_starts"])
            self.assertEqual(c["kept"]["2500"]["h"]["starts_reaching_h"],
                             {5: 1, 8: 5}[c["target_h"]])
            for lim in ("400", "2500"):
                for label in ("h_minus_1", "h_plus_1"):
                    self.assertEqual(c["reread"][lim][label]["starts_reaching_h"], 0)
                self.assertEqual(c["reread"][lim]["h"]["starts_holding_anchor"], 8)

    def test_kept_run_values_and_new_pairing_fields(self):
        rec = emlib.load_run(emlib.RUNS / "constrained_generate_seed0.json")
        for s in rec["outputs"]["shapes"]:
            self.assertEqual(s["compare"]["constrained_bound_check"]["h"]["starts_reaching_h"],
                             {5: 1, 8: 5}[s["h"]])
            self.assertEqual(len(s["compare"]["constrained"]["h"]["anchor_points_on_axis"]), 8)

    def test_tour_maximum_carries_the_degeneracy_tag(self):
        rec = emlib.load_run(emlib.RUNS / "mirror_sweep_seed0.json")
        for s in rec["outputs"]["sweep"]:
            self.assertTrue(s["max_axis_count_degenerate"]["shape"])
            self.assertEqual(s["max_axis_count_excluding_degenerate"]["shape"], 0)
            self.assertIn(s["max_axis_count"]["shape"], (3, 4))        # value kept as read
        for f in rec["outputs"]["frame_zero"]:
            self.assertTrue(f["degeneracy"]["tags"]["1e-06"]["coincident"])

    def test_small_run(self):
        small = emlib.load_run(emlib.RUNS / "pairing_reread_small_seed0.json")
        full = emlib.load_run(emlib.RUNS / "pairing_reread_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(full["outputs"]))


if __name__ == "__main__":
    unittest.main()
