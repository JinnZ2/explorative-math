"""Tests for the seed-1 holdout and experiments/axis_split.py. CC0. Co-authored."""
import unittest

import _path  # noqa: F401
import anchor_selection as S
import axis_split as X
import emlib


class Holdout(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = emlib.load_run(emlib.RUNS / "anchor_selection_seed1.json")
        cls.by = {s["shape"]: s["compare"] for s in cls.rec["outputs"]["shapes"]}

    def test_nothing_retuned(self):
        self.assertEqual(S.PREFERRED, {"5-cell": 1, "16-cell": 2})
        self.assertEqual(S.PREFERRED_FROM_SEED, 0)
        self.assertTrue(self.rec["params"]["holdout"])
        self.assertEqual(self.rec["params"]["preferred_points_on_axis"], S.PREFERRED)
        self.assertEqual(self.rec["params"]["max_attempts"], S.MAX_ATTEMPTS)
        seed0 = emlib.load_run(emlib.RUNS / "anchor_selection_seed0.json")
        self.assertFalse(seed0["params"]["holdout"])
        for k in ("trial_limits", "max_attempts", "preferred_points_on_axis", "rules"):
            self.assertEqual(self.rec["params"][k], seed0["params"][k], k)

    def test_reference_is_the_same_seed_and_starts_differ_from_seed_0(self):
        self.assertEqual(self.rec["params"]["reference"], "runs/constrained_generate_seed1.json")
        seed0 = emlib.load_run(emlib.RUNS / "anchor_selection_seed0.json")
        for s1, s0 in zip(self.rec["outputs"]["shapes"], seed0["outputs"]["shapes"]):
            self.assertTrue(s1["compare"]["starts_equal_constrained_generate_starts"])
            self.assertNotEqual(s1["starts"][0]["start_angles_deg"],
                                s0["starts"][0]["start_angles_deg"])

    def test_observed_reach_on_the_holdout(self):
        # OBSERVED on the committed run, pinned so a change shows up
        got = {n: (c["selected"]["2500"]["h"]["starts_reaching_h"],
                   c["random"]["2500"]["h"]["starts_reaching_h"]) for n, c in self.by.items()}
        self.assertEqual(got, {"5-cell": (7, 4), "16-cell": (1, 2)})
        self.assertIsNone(self.rec["ai_interest"])


class AxisSplit(unittest.TestCase):
    def test_labels_and_pairing_description(self):
        self.assertEqual([X.label(i) for i in range(4)], ["+e0", "-e0", "+e1", "-e1"])
        d = X.describe_pairing([0, 1, 3, 2, 6, 7, 4, 5])
        self.assertEqual(d["points_on_axis"], ["+e0", "-e0"])
        self.assertTrue(d["on_axis_is_antipodal_pair"])
        self.assertEqual(d["swaps_with_own_antipode"], 1)        # +e1 <-> -e1
        self.assertEqual(d["swaps_between_generators"], 2)       # e2 <-> e3, twice
        self.assertEqual(d["generator_permutation"], [(2, 3)])

    def test_coxeter_reference(self):
        c = X.coxeter_reference()
        self.assertEqual(c["n_mirror_axes"], 8)
        self.assertAlmostEqual(c["generators"]["radius_ratio_min_over_max"], 1.0, places=9)
        self.assertEqual(sorted(len(m["on_axis_indices"]) for m in c["mirrors"]),
                         [0, 0, 0, 0, 2, 2, 2, 2])

    def test_logged_runs_remake_the_selection_counts_and_keep_the_path(self):
        for seed, n_reach in ((0, 5), (1, 1)):
            rec = emlib.load_run(emlib.RUNS / f"axis_split_seed{seed}.json")
            o = rec["outputs"]
            self.assertIsNone(rec["ai_interest"])
            self.assertTrue(o["counts_equal_logged_selection_run"])
            self.assertTrue(o["all_anchors_have_preferred_type"])
            self.assertEqual(len(o["reached"]["starts"]), n_reach)
            self.assertEqual(len(o["reached"]["starts"]) + len(o["not_reached"]["starts"]), 8)
            for r in o["starts"]:
                self.assertEqual(len(r["residual_path"]), len(r["angles_path"]))
                self.assertGreater(len(r["residual_path"]), 1)
                for k in ("anchor_pairing", "generators_at_anchor", "generators_at_end",
                          "degeneracy_at_anchor", "degeneracy_at_end"):
                    self.assertIn(k, r)
                self.assertTrue(r["anchor_pairing"]["on_axis_is_antipodal_pair"])

    def test_small_run(self):
        out = X.run(seed=0, starts=1, trials=40, anchor_max_evals=1000, max_attempts=6)
        self.assertEqual(len(out["starts"]), 1)
        self.assertEqual(set(out), set(emlib.load_run(
            emlib.RUNS / "axis_split_seed0.json")["outputs"]))


if __name__ == "__main__":
    unittest.main()
