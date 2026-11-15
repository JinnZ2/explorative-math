"""Tests for probes.min_stretch_tag and experiments/min_stretch_flag.py. CC0. Co-authored."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import _path  # noqa: F401
import min_stretch_flag as M
from toolkit import probes as Q


class MinStretchTag(unittest.TestCase):
    def test_tag_and_threshold(self):
        r = Q.min_stretch_tag(57, 0.01)
        self.assertEqual(r["tag"], "under-resolved")
        self.assertAlmostEqual(r["min_points"], 300.0)
        self.assertIsNone(Q.min_stretch_tag(300, 0.01)["tag"])
        self.assertTrue(Q.min_stretch_tag(300, 0.01)["assessed"])

    def test_unknown_frequency_is_not_assessed(self):
        for f in (None, 0, 0.0):
            r = Q.min_stretch_tag(57, f)
            self.assertFalse(r["assessed"])
            self.assertIsNone(r["tag"])
            self.assertIsNone(r["min_points"])


class MinStretchFlag(unittest.TestCase):
    def test_sources_of_frequency(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "short_seed0.json").write_text(json.dumps({"outputs": {
                "two": {"distinct_gap_lengths": [1.0, 2.0], "word": "AB" * 10},
                "three": {"distinct_gap_lengths": [1.0, 2.0, 3.0], "word": "ABC" + "B" * 97},
                "void": {"distinct_gap_lengths": None}}}))
            (d / "long_seed0.json").write_text(json.dumps({"outputs": {"s": [
                {"distinct_gap_lengths": [{"length": 1.0, "count": 10},
                                          {"length": 2.0, "count": 890},
                                          {"length": 3.0, "count": 100}]}]}}))
            out = M.run(runs_dir=d, groups={},
                        references={("short_seed0.json", "two"): ("long_seed0.json", "s[0]")})
        rows = {(r["run"], r["where"]): r for r in out["rows"]}
        two = rows[("short_seed0.json", "two")]
        self.assertEqual((two["freq_source"], two["tag"]), ("reference", "under-resolved"))
        self.assertAlmostEqual(two["third_gap_freq"], 0.01)
        three = rows[("short_seed0.json", "three")]
        self.assertEqual((three["freq_source"], three["tag"]), ("self", "under-resolved"))
        void = rows[("short_seed0.json", "void")]
        self.assertEqual((void["freq_source"], void["assessed"], void["tag"]),
                         ("void", False, None))
        self.assertIsNone(rows[("long_seed0.json", "s[0]")]["tag"])
        self.assertEqual(out["n_rows"], 4)

    def test_logged_runs_flag_and_are_not_edited(self):
        files = M.emlib.run_files()
        before = {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
        out = M.run()
        after = {f.name: hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
        self.assertEqual(before, after)
        # OBSERVED on the committed runs, pinned so a change shows up
        tagged = {(r["run"], r["where"]) for r in out["under_resolved"]}
        self.assertIn(("generated_cut_weighted_seed1.json", "levels[0]"), tagged)
        self.assertEqual(out["n_under_resolved"], 4)
        self.assertEqual(out["n_rows"], out["n_assessed"] + out["n_not_assessed"])


if __name__ == "__main__":
    unittest.main()
