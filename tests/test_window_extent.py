"""Tests for window_extent, three_gap_check and probes.three_gap_sum. CC0. Co-authored."""
import json
import tempfile
import unittest
from pathlib import Path

import _path  # noqa: F401
import three_gap_check as T
import window_extent as E
from toolkit import probes as Q


class ThreeGapSum(unittest.TestCase):
    def test_pass_fail_not_applicable(self):
        self.assertEqual(Q.three_gap_sum([1.0, 2.0, 3.0])["verdict"], "pass")
        self.assertEqual(Q.three_gap_sum([3.0, 1.0, 2.0])["verdict"], "pass")
        f = Q.three_gap_sum([1.0, 2.0, 3.5])
        self.assertEqual(f["verdict"], "fail")
        self.assertAlmostEqual(f["residual"], 0.5)
        for ls in ([], [1.0], [1.0, 2.0], [1.0, 2.0, 3.0, 4.0]):
            r = Q.three_gap_sum(ls)
            self.assertEqual(r["verdict"], "not_applicable")
            self.assertIsNone(r["residual"])


class WindowExtent(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = E.run(targets=(57, 200, 500))

    def test_window_matches_source_run(self):
        src = json.loads((E.emlib.ROOT / E.SOURCE).read_text())
        lv = src["outputs"]["levels"][E.LEVEL_INDEX]
        self.assertAlmostEqual(self.out["source_window"]["half_width_in_q"],
                               lv["half_width_in_q"], places=12)
        # the 57-point stretch is the source run's own stretch
        self.assertEqual(self.out["stretches"][0]["n_points"], lv["n_points"])

    def test_stretches_nested_and_at_most_three_lengths(self):
        ns = [s["n_points"] for s in self.out["stretches"]]
        self.assertEqual(ns, sorted(ns))
        for s in self.out["stretches"]:
            self.assertLessEqual(s["n_gap_lengths"], 3)
            self.assertNotEqual(s["three_gap_sum"]["verdict"], "fail")

    def test_observed_counts(self):
        # OBSERVED on this run, pinned so a change shows up
        self.assertEqual([s["n_gap_lengths"] for s in self.out["stretches"]], [2, 3, 3])

    def test_control_canonical_width_is_one_plus_slope(self):
        r = self.out["combinations"]["control_canonical_width_times_norm"]["nearest"][0]
        self.assertEqual((r["m"], r["n"]), (1, 1))
        self.assertAlmostEqual(r["residual"], 0.0, places=6)

    def test_nearest_combos_sorted(self):
        rows = E.nearest_combos(0.5, 0.25)
        self.assertEqual(abs(rows[0]["residual"]), 0.0)
        res = [abs(r["residual"]) for r in rows]
        self.assertEqual(res, sorted(res))


class ThreeGapCheck(unittest.TestCase):
    def test_rows_kept_for_every_entry(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "a_seed0.json").write_text(json.dumps({"outputs": {
                "x": {"distinct_gap_lengths": [1.0, 2.0, 3.0]},
                "y": [{"distinct_gap_lengths": [{"length": 1.0, "count": 2},
                                                {"length": 2.0, "count": 1}]}],
                "z": {"distinct_gap_lengths": [1.0, 2.0, 4.0]},
                "v": {"distinct_gap_lengths": None}}}))
            (d / "three_gap_check_seed0.json").write_text("{}")   # own output: skipped
            out = T.run(runs_dir=d)
        self.assertEqual([r["verdict"] for r in out["rows"]],
                         ["pass", "not_applicable", "fail", "not_applicable"])
        self.assertEqual(out["rows"][1]["where"], "y[0].distinct_gap_lengths")

    def test_no_fail_in_logged_runs(self):
        # OBSERVED on the committed runs, pinned so a change shows up
        out = T.run()
        self.assertEqual(out["verdict_counts"].get("fail", 0), 0)
        self.assertGreater(out["verdict_counts"]["pass"], 0)


if __name__ == "__main__":
    unittest.main()
