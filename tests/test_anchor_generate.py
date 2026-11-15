"""Tests for probes.mirror_residual_at and experiments/anchor_generate.py. CC0. Co-authored."""
import unittest

import _path  # noqa: F401
import anchor_generate as A
import emlib
from toolkit import primitives as P
from toolkit import probes as Q


class ResidualAt(unittest.TestCase):
    def test_zero_on_mirror_axes_and_positive_off_them(self):
        sq = P.regular_polygon(4)                  # mirrors at 0, 45, 90, 135
        r = Q.mirror_residual_at(sq, [0.0, 45.0, 20.0])
        self.assertAlmostEqual(r[0], 0.0, places=12)
        self.assertAlmostEqual(r[1], 0.0, places=12)
        self.assertGreater(r[2], 1e-2)

    def test_order_and_length_follow_the_axes_given(self):
        pts = [[1.0, 0.0], [-1.0, 0.2], [0.3, 1.0], [0.9, -1.4]]
        a = Q.mirror_residual_at(pts, [10.0, 70.0, 130.0])
        b = Q.mirror_residual_at(pts, [130.0, 10.0])
        self.assertEqual(len(a), 3)
        self.assertEqual([a[2], a[0]], b)
        self.assertEqual(Q.mirror_residual_at(pts, []), [])


class GeneratedAxes(unittest.TestCase):
    def test_dihedral_rule(self):
        ax = A.generated_axes(10.0, 5)
        self.assertEqual(len(ax), 5)
        for got, want in zip(ax, [10.0, 46.0, 82.0, 118.0, 154.0]):
            self.assertAlmostEqual(got, want)
        self.assertTrue(all(0.0 <= x < 180.0 for x in A.generated_axes(170.0, 8)))


class AnchorGenerate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = A.run(seed=0, starts=1, max_evals=80, anchor_max_evals=300,
                        bound_check_evals=None)

    def test_three_rules_from_the_anchor(self):
        for s in self.out["shapes"]:
            r = s["starts"][0]
            self.assertEqual([r["rules"][k]["m"] for k in ("h_minus_1", "h", "h_plus_1")],
                             [s["h"] - 1, s["h"], s["h"] + 1])
            for rule in r["rules"].values():
                self.assertEqual(rule["path"][0]["angles_deg"], r["anchor"]["final_angles_deg"])
                self.assertAlmostEqual(rule["generated_axes_deg"][0],
                                       r["anchor"]["axis_deg"] % 180.0)
                sc = [p["score"] for p in rule["path"]]
                self.assertEqual(sc, sorted(sc, reverse=True))
                self.assertEqual(len(rule["final_residual_per_axis"]), rule["m"])
                self.assertAlmostEqual(sum(rule["final_residual_per_axis"]), rule["final_score"])

    def test_bound_check_off_is_void(self):
        for s in self.out["shapes"]:
            self.assertIsNone(s["compare"]["anchor_generate"]["h"]["bound_check"])

    def test_missing_layered_reference_is_void(self):
        out = A.run(seed=0, starts=1, max_evals=20, anchor_max_evals=60,
                    layered_path="/nonexistent.json", bound_check_evals=None)
        for s in out["shapes"]:
            self.assertIsNone(s["compare"]["layered_unanchored"])
            self.assertIsNone(s["compare"]["alt"])
            self.assertIsNone(s["compare"]["starts_equal_layered_pursuit_starts"])


class LoggedRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = emlib.load_run(emlib.RUNS / "anchor_generate_seed0.json")
        cls.by = {s["shape"]: s["compare"] for s in cls.rec["outputs"]["shapes"]}

    def test_gz_and_observer_slots_null(self):
        self.assertTrue((emlib.RUNS / "anchor_generate_seed0.json.gz").exists())
        self.assertFalse((emlib.RUNS / "anchor_generate_seed0.json").exists())
        self.assertIsNone(self.rec["human_interest"])
        self.assertIsNone(self.rec["ai_interest"])
        self.assertIsNone(self.rec["divergence_note"])

    def test_same_starts_and_bounds_logged(self):
        for c in self.by.values():
            self.assertTrue(c["starts_equal_layered_pursuit_starts"])
        b = self.rec["outputs"]["bound"]
        self.assertEqual(b["generate"]["max_evals"], A.MAX_EVALS)
        self.assertEqual(b["generate"]["bound_check_max_evals"], A.BOUND_CHECK_EVALS)
        self.assertIn("evals_not_like_for_like", b)

    def test_observed_counts(self):
        # OBSERVED on the committed run, pinned so a change shows up
        got = {name: {k: c["anchor_generate"][k]["starts_reaching_h"]
                      for k in ("h_minus_1", "h", "h_plus_1")} for name, c in self.by.items()}
        self.assertEqual(got, {"5-cell": {"h_minus_1": 0, "h": 2, "h_plus_1": 0},
                               "16-cell": {"h_minus_1": 0, "h": 1, "h_plus_1": 0}})
        self.assertEqual(self.by["16-cell"]["anchor_generate"]["h"]["bound_check"]
                         ["starts_reaching_h"], 2)
        self.assertEqual(self.by["5-cell"]["layered_unanchored"]["starts_reaching_h"], 0)
        self.assertEqual(self.by["16-cell"]["alt"]["starts_reaching_h"], 1)


class LabelsAndSmall(unittest.TestCase):
    def test_layered_run_is_labelled_unanchored_and_values_kept(self):
        rec = emlib.load_run(emlib.RUNS / "layered_pursuit_seed0.json")
        self.assertEqual(rec["params"]["label"], "layered_unanchored")
        self.assertFalse(rec["params"]["anchored"])
        self.assertEqual(rec["outputs"]["label"], "layered_unanchored")
        by = {s["shape"]: s for s in rec["outputs"]["shapes"]}
        self.assertEqual(by["16-cell"]["alt_final_counts"], [0, 0, 0, 0, 0, 8, 0, 0])
        self.assertEqual(by["5-cell"]["layered_final_counts"], [0, 4, 0, 0, 0, 0, 0, 0])

    def test_anchor_kept_is_logged_and_no_constraint_is_claimed(self):
        rec = emlib.load_run(emlib.RUNS / "anchor_generate_seed0.json")
        self.assertIn("none", rec["outputs"]["bound"]["generate"]["anchor_constraint"])
        got = {s["shape"]: s["compare"]["anchor_generate"]["h"]["starts_keeping_anchor"]
               for s in rec["outputs"]["shapes"]}
        self.assertEqual(got, {"5-cell": 6, "16-cell": 3})   # OBSERVED, pinned
        for s in rec["outputs"]["shapes"]:
            self.assertIn("layered_unanchored", s["compare"])
            self.assertNotIn("layered", s["compare"])

    def test_small_runs_are_separate_files_with_every_output_key(self):
        for name in ("layered_pursuit", "anchor_generate"):
            full = emlib.load_run(emlib.RUNS / f"{name}_seed0.json")
            small = emlib.load_run(emlib.RUNS / f"{name}_small_seed0.json")
            self.assertTrue(small["params"]["small"])
            self.assertFalse(full["params"]["small"])
            self.assertLess(small["params"]["starts"], full["params"]["starts"])
            self.assertEqual(set(small["outputs"]), set(full["outputs"]))
            fs, ss = full["outputs"]["shapes"][0], small["outputs"]["shapes"][0]
            self.assertEqual(set(ss), set(fs))
            self.assertEqual(set(ss["starts"][0]), set(fs["starts"][0]))
            self.assertEqual(len(ss["starts"]), small["params"]["starts"])


if __name__ == "__main__":
    unittest.main()
