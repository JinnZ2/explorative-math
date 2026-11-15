"""Tests for experiments/fold_loop.py. CC0. Co-authored."""
import math
import unittest

import _path  # noqa: F401
import emlib
import fold_loop as F
from toolkit import primitives as P


class Wedge(unittest.TestCase):
    def test_width(self):
        self.assertEqual([F.wedge_width(n) for n in (1, 2, 3, 4, 5, 8)],
                         [90.0, 90.0, 30.0, 45.0, 18.0, 22.5])

    def test_in_wedge_wraps_at_180(self):
        self.assertTrue(F.in_wedge(5.0, 170.0, 30.0))
        self.assertFalse(F.in_wedge(25.0, 170.0, 30.0))
        self.assertTrue(F.in_wedge(170.0, 170.0, 30.0))
        self.assertAlmostEqual(F.circ(1.0, 179.0), 2.0)

    def test_score_reads_only_candidates_in_the_wedge_and_off_found_axes(self):
        sq = P.regular_polygon(4)                       # mirrors at 0, 45, 90, 135
        r, deg, k = F.wedge_score(sq, 0.0, 90.0, found=[0.0])
        self.assertAlmostEqual(r, 0.0, places=12)
        self.assertTrue(F.SEP_DEG < deg <= 90.0)
        _, _, k_all = F.wedge_score(sq, 0.0, 180.0, found=[])
        self.assertLess(k, k_all)
        none, deg, k = F.wedge_score(sq, 10.0, 20.0, found=[])   # no candidate in 10..30
        self.assertIsNone(none)
        self.assertEqual(k, 0)

    def test_pull_back_holds_two_mirrors_at_once(self):
        verts = P.cross_polytope(4)
        planes = P.all_planes(4)
        ang = [0.0] * 6                                  # axis-aligned: mirrors at 0 and 90
        pts = F.K.project(verts, planes, ang)
        held = [(0.0, F.K.pairing(pts, 0.0)), (90.0, F.K.pairing(pts, 90.0))]
        moved = [0.3, -0.2, 0.25, 0.1, -0.3, 0.2]
        back, _ = F.pull_back(verts, planes, moved, held)
        self.assertTrue(F.held_ok(verts, planes, back, held))


class Loop(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = F.run(seed=0, starts=2, level_trials=60, anchor_max_evals=1000, extra_starts=0)

    def test_every_endpoint_and_level_has_a_degeneracy_reading_and_held_mirrors(self):
        for s in self.out["shapes"]:
            for r in s["starts"]:
                for v in F.VARIANTS:
                    g = r["variants"][v]
                    self.assertIn("min_pair_distance", g["degeneracy"])
                    for lv in g["levels"]:
                        self.assertIn("degeneracy", lv)
                        self.assertLess(lv["max_held_residual"], F.ANCHOR_TOL)
                        for k in ("domain_size_deg", "axes_found_after", "trials",
                                  "projection_evals", "score_evals"):
                            self.assertIn(k, lv)

    def test_domain_sizes_by_variant(self):
        for s in self.out["shapes"]:
            for r in s["starts"]:
                for lv in r["variants"]["no_fold"]["levels"]:
                    self.assertEqual(lv["domain_size_deg"], 180.0)
                for v in ("fold", "wrong_fold"):
                    for lv in r["variants"][v]["levels"]:
                        self.assertEqual(lv["domain_size_deg"],
                                         F.wedge_width(max(lv["axes_found_before"], 1)))
                self.assertIsNotNone(r["variants"]["wrong_fold"]["wrong_fold_offset_deg"])
                self.assertIsNone(r["variants"]["fold"]["wrong_fold_offset_deg"])

    def test_extra_shape_not_run_is_logged(self):
        self.assertEqual(self.out["bound"]["extra_shapes_not_run"], ["24-cell"])
        self.assertEqual([s["shape"] for s in self.out["shapes"]], ["5-cell", "16-cell"])


class LoggedRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = emlib.load_run(emlib.RUNS / "fold_loop_seed0.json")
        cls.by = {s["shape"]: s["compare"] for s in cls.rec["outputs"]["shapes"]}

    def test_observer_slots_null_and_same_starts(self):
        self.assertIsNone(self.rec["ai_interest"])
        self.assertIsNone(self.rec["human_interest"])
        for name in ("5-cell", "16-cell"):
            self.assertEqual(self.by[name]["starts_equal_reference_starts"],
                             {"constrained_generate": True, "anchor_selection": True})

    def test_observed_counts(self):
        # OBSERVED on the committed run, pinned so a change shows up
        self.assertEqual(set(self.by), {"5-cell", "16-cell", "24-cell"})
        reach = {n: [c[k]["starts_reaching_h"] for k in ("b_fold", "c_wrong_fold", "d_no_fold")]
                 for n, c in self.by.items()}
        self.assertEqual(reach, {"5-cell": [1, 1, 1], "16-cell": [3, 3, 3], "24-cell": [0, 0, 0]})
        for c in self.by.values():
            self.assertEqual(c["b_fold"]["final_counts"], c["c_wrong_fold"]["final_counts"])
            self.assertEqual(c["b_fold"]["final_counts"], c["d_no_fold"]["final_counts"])
            self.assertLess(sum(c["b_fold"]["candidate_axes_total"]),
                            sum(c["d_no_fold"]["candidate_axes_total"]))
            self.assertEqual(c["b_fold"]["starts_reaching_h"],
                             c["b_fold"]["starts_reaching_h_excluding_degenerate"])
        a = self.by["16-cell"]["a_full_space_anchored_generate"]
        self.assertEqual(a["constrained_generate"]["starts_reaching_h"], 5)
        self.assertEqual(a["anchor_selection"]["starts_reaching_h"], 5)
        self.assertIsNone(self.by["24-cell"]["a_full_space_anchored_generate"]
                          ["constrained_generate"])

    def test_small_run(self):
        small = emlib.load_run(emlib.RUNS / "fold_loop_small_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(self.rec["outputs"]))
        self.assertTrue(math.isfinite(small["outputs"]["shapes"][0]["starts"][0]["anchor"]
                                      ["final_residual"]))


if __name__ == "__main__":
    unittest.main()
