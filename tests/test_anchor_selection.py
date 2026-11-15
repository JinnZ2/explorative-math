"""Tests for experiments/anchor_selection.py and docs/anchor.md pointers. CC0. Co-authored."""
import re
import unittest
from pathlib import Path

import _path  # noqa: F401
import anchor_selection as S
import emlib


class Selection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = S.run(seed=0, starts=2, limits=(40,), anchor_max_evals=1000, max_attempts=6)

    def test_attempts_stop_at_the_preferred_type_or_the_limit(self):
        for s in self.out["shapes"]:
            want = S.PREFERRED[s["shape"]]
            for r in s["starts"]:
                self.assertLessEqual(r["attempts"], 6)
                self.assertEqual(r["attempts"], len(r["attempt_types"]))
                self.assertEqual(r["re_anchored"], r["attempts"] > 1)
                self.assertTrue(all(t != want for t in r["attempt_types"][:-1]))
                self.assertEqual(r["selected_anchor_is_preferred_type"],
                                 r["attempt_types"][-1] == want)
                self.assertEqual(r["anchor"]["points_on_axis"], r["attempt_types"][-1])
                self.assertEqual(r["first_anchor_points_on_axis"], r["attempt_types"][0])

    def test_pairing_is_kept_and_anchor_held(self):
        for s in self.out["shapes"]:
            for r in s["starts"]:
                for label in S.RULES:
                    g = r["rules"][label]["40"]
                    self.assertEqual(g["pairing"], "kept from the anchor")
                    self.assertTrue(g["anchor_held"])

    def test_deterministic(self):
        again = S.run(seed=0, starts=2, limits=(40,), anchor_max_evals=1000, max_attempts=6)
        a = [r["attempt_types"] for s in self.out["shapes"] for r in s["starts"]]
        b = [r["attempt_types"] for s in again["shapes"] for r in s["starts"]]
        self.assertEqual(a, b)


class LoggedRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = emlib.load_run(emlib.RUNS / "anchor_selection_seed0.json")
        cls.by = {s["shape"]: s["compare"] for s in cls.rec["outputs"]["shapes"]}

    def test_observer_slots_null_and_same_starts(self):
        self.assertIsNone(self.rec["ai_interest"])
        self.assertIsNone(self.rec["human_interest"])
        for c in self.by.values():
            self.assertTrue(c["starts_equal_constrained_generate_starts"])
            self.assertEqual(c["starts_left_on_another_type"], [])

    def test_observed_reach(self):
        # OBSERVED on the committed run, pinned so a change shows up
        got = {n: (c["selected"]["2500"]["h"]["starts_reaching_h"],
                   c["random"]["2500"]["h"]["starts_reaching_h"]) for n, c in self.by.items()}
        self.assertEqual(got, {"5-cell": (4, 1), "16-cell": (5, 5)})
        for c in self.by.values():
            for lim in ("400", "2500"):
                for label in ("h_minus_1", "h_plus_1"):
                    self.assertEqual(c["selected"][lim][label]["starts_reaching_h"], 0)
        self.assertEqual(self.by["5-cell"]["re_anchored_starts"], [4, 5, 6, 7])
        self.assertEqual(self.by["16-cell"]["re_anchored_starts"], [0, 1, 4, 6])

    def test_unselected_starts_repeat_the_random_run(self):
        for c in self.by.values():
            sel = c["selected"]["2500"]["h"]["final_counts"]
            rnd = c["random"]["2500"]["h"]["final_counts"]
            for i, (a, b) in enumerate(zip(sel, rnd)):
                if i not in c["re_anchored_starts"]:
                    self.assertEqual(a, b, i)

    def test_small_run(self):
        small = emlib.load_run(emlib.RUNS / "anchor_selection_small_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(self.rec["outputs"]))


class AnchorDoc(unittest.TestCase):
    def test_statement_and_run_pointers_exist(self):
        text = (emlib.ROOT / "docs" / "anchor.md").read_text()
        self.assertIn("anchor = position + pairing (relation set)", text)
        self.assertIn("Re-reading the pairing permits collapse.", text)
        for name in re.findall(r"`(runs/[^`]+)`", text):
            self.assertTrue((emlib.ROOT / name).exists() or
                            (emlib.ROOT / (name + ".gz")).exists(), name)


if __name__ == "__main__":
    unittest.main()
