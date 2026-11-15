"""Tests for gz run files, the k-axis score and layered_pursuit. CC0. Co-authored."""
import gzip
import json
import tempfile
import unittest
from pathlib import Path

import _path  # noqa: F401
import emlib
import layered_pursuit as L
from toolkit import primitives as P
from toolkit import probes as Q


class GzRuns(unittest.TestCase):
    def _write(self, d, big, gz=None):
        outputs = {"x": list(range(60000 if big else 10))}
        return emlib.write_run("t", 0, {}, outputs, None, None, None, None, out_dir=d, gz=gz)

    def test_small_is_plain_and_large_is_gz_with_full_content(self):
        with tempfile.TemporaryDirectory() as d:
            small = self._write(d, False)
            self.assertEqual(small.name, "t_seed0.json")
            big = self._write(d, True)
            self.assertEqual(big.name, "t_seed0.json.gz")
            self.assertFalse(small.exists())            # one file per run
            rec = emlib.load_run(big)
            self.assertEqual(rec["outputs"]["x"], list(range(60000)))
            self.assertEqual(tuple(rec), emlib.SCHEMA)
            with gzip.open(big, "rt", encoding="utf-8") as fh:
                self.assertEqual(json.load(fh), rec)

    def test_loader_reads_either_name(self):
        with tempfile.TemporaryDirectory() as d:
            big = self._write(d, True)
            self.assertEqual(emlib.load_run(Path(d) / "t_seed0.json"), emlib.load_run(big))
            plain = self._write(d, True, gz=False)
            self.assertEqual(plain.name, "t_seed0.json")
            self.assertFalse(big.exists())
            self.assertEqual(emlib.load_run(Path(d) / "t_seed0.json.gz"), emlib.load_run(plain))
            with self.assertRaises(FileNotFoundError):
                emlib.load_run(Path(d) / "absent_seed0.json")

    def test_gz_bytes_are_reproducible(self):
        with tempfile.TemporaryDirectory() as d:
            a = self._write(d, True).read_bytes()
            b = self._write(d, True).read_bytes()
            self.assertEqual(a, b)

    def test_run_files_lists_both_forms(self):
        with tempfile.TemporaryDirectory() as d:
            self._write(d, True)
            emlib.write_run("u", 0, {}, {}, None, None, None, None, out_dir=d)
            self.assertEqual([p.name for p in emlib.run_files(d)],
                             ["t_seed0.json.gz", "u_seed0.json"])

    def test_committed_runs_all_load(self):
        files = emlib.run_files()
        self.assertTrue(any(f.name.endswith(".json.gz") for f in files))
        names = [f.name[:-3] if f.name.endswith(".gz") else f.name for f in files]
        self.assertEqual(len(names), len(set(names)))    # no run in both forms
        for f in files:
            self.assertEqual(tuple(emlib.load_run(f)), emlib.SCHEMA)


class KScore(unittest.TestCase):
    def test_k_mirrors_score_zero_and_one_more_does_not(self):
        sq = P.regular_polygon(4, phase=0.2)
        self.assertAlmostEqual(Q.mirror_k_score(sq, 4), 0.0, places=12)
        self.assertIsNone(Q.mirror_k_score(sq, 5))       # a square offers 4 axes only
        hexa = P.regular_polygon(6, phase=0.2)
        self.assertAlmostEqual(Q.mirror_k_score(hexa, 6), 0.0, places=12)
        lopsided = [[1.0, 0.0], [-1.0, 0.2], [0.3, 1.0], [0.9, -1.4]]
        self.assertGreater(Q.mirror_k_score(lopsided, 2), 1e-3)

    def test_axes_are_distinct(self):
        rows = Q.mirror_axis_residuals(P.regular_polygon(5, phase=0.1), min_sep_deg=5.0)
        axes = sorted(r["axis_deg"] for r in rows)
        self.assertTrue(all(b - a > 5.0 for a, b in zip(axes, axes[1:])))
        self.assertEqual(sum(1 for r in rows if r["residual"] < 1e-12), 5)
        self.assertEqual([r["residual"] for r in rows], sorted(r["residual"] for r in rows))

    def test_too_few_candidate_axes_is_void(self):
        self.assertIsNone(Q.mirror_k_score([[1.0, 0.0], [-1.0, 0.0]], 5))

    def test_first_row_is_the_mirror_residual(self):
        pts = [[1.0, 0.0], [-1.0, 0.2], [0.3, 1.0], [0.9, -1.4]]
        self.assertAlmostEqual(Q.mirror_axis_residuals(pts)[0]["residual"],
                               Q.mirror_residual(pts)["residual"])


class LayeredPursuit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = L.run(seed=0, starts=1, max_evals=120)

    def test_stage_structure(self):
        for s in self.out["shapes"]:
            r = s["starts"][0]
            self.assertEqual([st["k"] for st in r["layered_stages"]], list(range(1, s["h"] + 1)))
            self.assertEqual(len(r["layered_axis_count_per_stage"]), s["h"])
            self.assertEqual(r["alt"]["k"], s["h"])
            for a, b in zip(r["layered_stages"], r["layered_stages"][1:]):
                self.assertEqual(b["path"][0]["angles_deg"], a["final_angles_deg"])
            self.assertEqual(r["layered_stages"][0]["path"][0]["angles_deg"],
                             r["start_angles_deg"])
            self.assertEqual(r["alt"]["path"][0]["angles_deg"], r["start_angles_deg"])

    def test_scores_never_rise_within_a_descent(self):
        for s in self.out["shapes"]:
            r = s["starts"][0]
            for st in r["layered_stages"] + [r["alt"]]:
                sc = [p["score"] for p in st["path"]]
                self.assertEqual(sc, sorted(sc, reverse=True))

    def test_same_starts_as_mirror_pursuit(self):
        # with one start per shape only the first shape's draw lines up
        self.assertTrue(self.out["shapes"][0]["starts_equal_mirror_pursuit_starts"])
        rec = emlib.load_run(emlib.RUNS / "layered_pursuit_seed0.json")
        for s in rec["outputs"]["shapes"]:
            self.assertTrue(s["starts_equal_mirror_pursuit_starts"])

    def test_observed_counts_in_the_logged_run(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = emlib.load_run(emlib.RUNS / "layered_pursuit_seed0.json")
        by = {s["shape"]: s for s in rec["outputs"]["shapes"]}
        self.assertEqual(by["5-cell"]["layered_starts_reaching_coxeter_count"], 0)
        self.assertEqual(by["5-cell"]["alt_starts_reaching_coxeter_count"], 0)
        self.assertEqual(by["16-cell"]["layered_starts_reaching_coxeter_count"], 0)
        self.assertEqual(by["16-cell"]["alt_final_counts"], [0, 0, 0, 0, 0, 8, 0, 0])
        for s in by.values():
            self.assertEqual(len(s["starts"]), 8)

    def test_missing_references_are_void(self):
        out = L.run(seed=0, starts=1, max_evals=30, pursuit_path="/nonexistent.json",
                    coxeter_path="/nonexistent.json")
        for s in out["shapes"]:
            self.assertIsNone(s["coxeter_plane_mirror_count"])
            self.assertIsNone(s["starts_equal_mirror_pursuit_starts"])
            self.assertIsNone(s["alt_starts_reaching_coxeter_count"])


if __name__ == "__main__":
    unittest.main()
