"""Tests for experiments/combine.py, the weighed mirror test and the centre reading.
CC0. Co-authored."""
import unittest

import _path  # noqa: F401
import combine as C
import emlib
from exact import Cyc5, ZETA
from exact import fold as XF
from exact import inflation as XI
from toolkit import primitives as P
from toolkit import probes as Q


class Weighing(unittest.TestCase):
    def test_exact_pentagon_folded_mirror_does_not_keep_multiplicities(self):
        f = XF.half_plane_fold(XI.pentagon(), Cyc5(1, 0))
        self.assertTrue(XF.is_mirror([z for z, _ in f], ZETA))
        self.assertFalse(XF.is_mirror_weighed(f, ZETA))
        # a set whose mirror does keep them
        sq = [(Cyc5(1, 0), 2), (Cyc5(-1, 0), 2)]
        self.assertTrue(XF.is_mirror_weighed(sq, ZETA * ZETA.conjugate() * Cyc5(1, 0)))

    def test_weighed_residual_matches_only_equal_multiplicity(self):
        reps = [[1.0, 0.0], [0.0, 1.0]]            # mirror at 45 degrees swaps them
        self.assertAlmostEqual(C.residual_on_axis(reps, [1, 1], 1.0, 45.0, True), 0.0)
        self.assertAlmostEqual(C.residual_on_axis(reps, [2, 1], 1.0, 45.0, False), 0.0)
        self.assertGreater(C.residual_on_axis(reps, [2, 1], 1.0, 45.0, True), 0.5)

    def test_rule_candidates_stay_inside_the_wedge_and_off_found_axes(self):
        score = C.make_score(True, "rule", False, 0.0, 8)
        pts = P.regular_polygon(8)
        r, deg, k, n = score(pts, [0.0, 90.0])        # wedge 0..90, found axes excluded
        self.assertEqual(k, 3)                         # 22.5, 45, 67.5
        self.assertAlmostEqual(r, 0.0, places=9)
        score_nf = C.make_score(False, "rule", False, 0.0, 8)
        self.assertEqual(score_nf(pts, [0.0, 90.0])[2], 6)


class Centre(unittest.TestCase):
    def test_square_with_a_centre_point(self):
        pts = P.regular_polygon(4) + [[0.0, 0.0]]
        d = Q.degeneracy(pts)
        self.assertAlmostEqual(d["min_radius"], 0.0, places=12)
        self.assertTrue(d["tags"]["1e-06"]["central"])
        self.assertEqual(d["tags"]["1e-06"]["n_at_centre"], 1)
        self.assertFalse(d["tags"]["1e-06"]["degenerate"])      # no two points coincide
        self.assertEqual(len(Q.mirror_symmetry_from_pairs(pts)), 4)

    def test_no_centre_point(self):
        d = Q.degeneracy(P.regular_polygon(5))
        self.assertFalse(d["tags"]["1e-06"]["central"])
        self.assertAlmostEqual(d["min_radius"], 1.0)

    def test_logged_five_cell_start_5(self):
        # OBSERVED on the committed runs, pinned so a change shows up
        rec = emlib.load_run(emlib.RUNS / "degeneracy_check_seed0.json")
        row = next(r for r in rec["outputs"]["rows"]
                   if r["run"] == "fold_shape_seed0.json" and r["shape"] == "5-cell"
                   and r["where"] == "starts[5].variants.fold_shape")
        self.assertEqual((row["n_mirror_axes"], row["degenerate"], row["central"],
                          row["n_at_centre"]), (4, False, True, 1))
        self.assertEqual(sum(g["reaching_h_central"] for g in rec["outputs"]["recount"]), 0)


class Combine(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = C.run(seed=0, starts=1, level_trials=40, anchor_max_evals=1000, max_attempts=4)

    def test_every_variant_runs_and_logs_cost_and_degeneracy(self):
        names = [v[0] for v in C.VARIANTS]
        for s in self.out["shapes"]:
            self.assertEqual(list(s["compare"]["variants"]), names)
            for g in s["starts"][0]["variants"].values():
                self.assertIn("min_radius", g["degeneracy"])
                for k in ("trials_total", "candidate_axes_total", "anchor_evals", "final_count"):
                    self.assertIn(k, g)
        self.assertEqual(self.out["bound"]["shapes_not_run"], ["24-cell"])
        self.assertFalse(self.out["exact_weighing_check"]["that_mirror_keeps_multiplicities"])


class LoggedRuns(unittest.TestCase):
    def test_observed_reach(self):
        # OBSERVED on the committed runs, pinned so a change shows up
        want = {0: {"5-cell": [1, 0, 3, 1, 5, 4, 0], "16-cell": [3, 4, 3, 4, 5, 4, 5]},
                1: {"5-cell": [1, 0, 1, 1, 6, 4, 0], "16-cell": [2, 2, 2, 2, 1, 2, 2]}}
        for seed, per in want.items():
            rec = emlib.load_run(emlib.RUNS / f"combine_seed{seed}.json")
            self.assertIsNone(rec["ai_interest"])
            self.assertEqual(rec["params"]["holdout"], seed != 0)
            for s in rec["outputs"]["shapes"]:
                got = [s["compare"]["variants"][v[0]]["starts_reaching_h"] for v in C.VARIANTS]
                self.assertEqual(got, per[s["shape"]], (seed, s["shape"]))
                for v in C.VARIANTS:
                    g = s["compare"]["variants"][v[0]]
                    self.assertEqual(g["starts_reaching_h"],
                                     g["starts_reaching_h_excluding_degenerate"])

    def test_fold_only_repeats_the_fold_shape_run(self):
        rec = emlib.load_run(emlib.RUNS / "combine_seed0.json")
        for s in rec["outputs"]["shapes"]:
            self.assertTrue(s["compare"]["fold_only_repeats_the_fold_shape_run"])

    def test_small_run(self):
        small = emlib.load_run(emlib.RUNS / "combine_small_seed0.json")
        full = emlib.load_run(emlib.RUNS / "combine_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(full["outputs"]))


if __name__ == "__main__":
    unittest.main()
