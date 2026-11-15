"""Tests for experiments/dual_partner.py and the carried hull name. CC0. Co-authored."""
import unittest

import _path  # noqa: F401
import dual_partner as D
import emlib


class Pieces(unittest.TestCase):
    def test_ranks_share_a_place_on_a_tie(self):
        rows = [{"n": "a", "v": 60.0}, {"n": "b", "v": 60.01}, {"n": "c", "v": 45.0}]
        D.rank_with_ties(rows, "v")
        self.assertEqual({r["n"]: r["v_rank"] for r in rows}, {"a": 1, "b": 1, "c": 3})
        D.rank_with_ties(rows, "v", reverse=False)
        self.assertEqual({r["n"]: r["v_rank"] for r in rows}, {"c": 1, "a": 2, "b": 2})

    def test_polygon_separation(self):
        sep, turn = D.polygon_separation(5, 5)
        self.assertAlmostEqual(sep, 36.0, places=6)
        self.assertAlmostEqual(turn, 36.0, places=6)
        self.assertAlmostEqual(D.polygon_separation(6, 3)[0], 30.0, places=6)
        self.assertAlmostEqual(D.polygon_separation(8, 4)[0], 22.5, places=6)

    def test_plane_part(self):
        pl = D.plane_part((3, 4, 5), (9,))
        by = {b["baseline"]: b for b in pl["baselines"]}
        self.assertTrue(by["5-gon"]["dual_is_the_only_best"])
        self.assertIsNone(by["9-gon"]["dual_is_best_separated"])       # control: void, not False
        self.assertFalse(by["9-gon"]["dual_in_probe_set"])


class LoggedRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = emlib.load_run(emlib.RUNS / "dual_partner_seed0.json")
        cls.o = cls.rec["outputs"]

    def test_observer_slots_null(self):
        self.assertIsNone(self.rec["ai_interest"])
        self.assertIsNone(self.rec["human_interest"])

    def test_observed_summary(self):
        # OBSERVED on the committed run, pinned so a change shows up
        s = self.o["summary"]
        self.assertEqual((s["sphere_baselines"], s["sphere_dual_best_separated"],
                          s["sphere_dual_the_only_best"],
                          s["sphere_dual_sits_in_and_covers_voids"]), (4, 4, 2, 4))
        self.assertEqual(s["sphere_non_dual_probes_that_also_cover"],
                         {"5-cell": [], "16-cell": [], "tesseract": [], "24-cell": []})
        self.assertEqual(s["sphere_controls_with_a_covering_probe"],
                         {"random_8": [], "rectified 5-cell": ["5-cell"],
                          "5-cell + its negative (as built)": []})
        self.assertEqual((s["plane_baselines"], s["plane_dual_best_separated"],
                          s["plane_dual_the_only_best"]), (7, 7, 4))

    def test_every_baseline_logs_ranks_and_search_effort(self):
        for b in self.o["sphere"]["baselines"]:
            self.assertEqual(len(b["probes"]), 4)
            for r in b["probes"]:
                for k in ("separation_deg", "separation_deg_rank", "residual_gap_deg",
                          "residual_gap_deg_rank", "sits_in_and_covers_the_voids", "search"):
                    self.assertIn(k, r)
            self.assertEqual(b["control"], b["dual"] is None)
        for b in self.o["plane"]["baselines"]:
            self.assertEqual(b["control"], not b["dual_in_probe_set"])

    def test_small_run(self):
        small = emlib.load_run(emlib.RUNS / "dual_partner_small_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(self.o))

    def test_carried_hull_name_is_marked_unchecked(self):
        rec = emlib.load_run(emlib.RUNS / "sphere_compounds_seed0.json")
        c = rec["outputs"]["catalogue_additions"]["compound of two 5-cells"]
        self.assertEqual(c["hull_name_carried"],
                         "disphenoidal 30-cell (dual of bitruncated 5-cell)")
        self.assertFalse(c["hull_name_checked_against_a_source"])
        self.assertEqual((c["computed"]["vertices"], c["computed"]["cells"]), (10, 30))
        self.assertTrue(c["computed"]["cells_are_tetragonal_disphenoids"])
        self.assertIsNone(rec["outputs"]["catalogue_additions"]["disphenoidal 288-cell"]
                          ["hull_name_carried"])

    def test_doc_carries_the_pattern_as_proposed(self):
        text = (emlib.ROOT / "docs" / "complement.md").read_text()
        self.assertIn("Pattern (PROPOSED)", text)
        self.assertIn("disphenoidal 30-cell", text)
        self.assertTrue((emlib.RUNS / "dual_partner_seed0.json.gz").exists())


if __name__ == "__main__":
    unittest.main()
