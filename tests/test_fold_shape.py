"""Tests for exact/fold.py and experiments/fold_shape.py. CC0. Co-authored."""
import math
import unittest

import _path  # noqa: F401
import emlib
import fold_shape as F
from exact import Cyc5, PHI, ZETA
from exact import fold as XF
from exact import inflation as XI
from toolkit import primitives as P


class ExactFold(unittest.TestCase):
    def setUp(self):
        self.pts = XI.pentagon()
        self.one = Cyc5(1, 0)
        self.z10 = -(ZETA * ZETA * ZETA)

    def test_reflect_is_an_involution_and_fixes_its_axis(self):
        for z in self.pts:
            self.assertEqual(XF.reflect(XF.reflect(z, ZETA), ZETA), z)
        self.assertEqual(XF.reflect(ZETA * PHI, ZETA), ZETA * PHI)
        self.assertEqual(XF.reflect(ZETA, self.one), ZETA.conjugate())

    def test_mirrors_of_the_pentagon(self):
        self.assertTrue(XF.is_mirror(self.pts, self.one))
        self.assertTrue(XF.is_mirror(self.pts, self.z10))
        self.assertFalse(XF.is_mirror(self.pts, self.one * 2 + ZETA))

    def test_one_mirror_fold_and_its_own_mirror(self):
        f = XF.half_plane_fold(self.pts, self.one)
        self.assertEqual(sorted(m for _, m in f), [1, 2, 2])
        self.assertEqual(sum(m for _, m in f), 5)
        fp = [z for z, _ in f]
        self.assertTrue(XF.is_mirror(fp, ZETA))             # 72 degrees
        self.assertFalse(XF.is_mirror(fp, self.z10))        # 36 degrees
        self.assertTrue(XF.is_mirror(self.pts, ZETA))       # and the shape has it

    def test_full_fold_is_a_single_point(self):
        f = XF.fold(self.pts, self.one, self.z10)
        self.assertEqual([(z, m) for z, m in f], [(self.one, 5)])

    def test_every_folded_point_lies_in_the_wedge(self):
        tris = XI.penrose_seed()
        for _ in range(3):
            tris = XI.penrose_subdivide(tris)
        pts = list({v for _, a, b, c in tris for v in (a, b, c)})
        for z, _ in XF.fold(pts, self.one, self.z10):
            self.assertGreaterEqual(XF.side(z, self.one), 0)
            self.assertLessEqual(XF.side(z, self.z10), 0)
        self.assertEqual(sum(m for _, m in XF.fold(pts, self.one, self.z10)), len(pts))

    def test_fold_does_not_check_that_the_wedge_is_a_mirror_wedge(self):
        # a point can land between two lines that are not mirrors of anything;
        # fold_one only reflects until the point is inside. The caller checks
        # the mirrors (is_mirror) before folding by them.
        odd = self.one * 3 + ZETA
        z, steps = XF.fold_one(ZETA * ZETA, self.one, odd)
        self.assertGreaterEqual(XF.side(z, self.one), 0)
        self.assertLessEqual(XF.side(z, odd), 0)
        self.assertFalse(XF.is_mirror(self.pts, odd))


class FloatFold(unittest.TestCase):
    def test_fold_by_mirrors_lands_in_the_wedge_and_keeps_radii(self):
        pts = P.regular_polygon(8, phase=0.3)
        deg0 = math.degrees(0.3)
        axes = [deg0 + 22.5 * k for k in range(8)]
        f = F.fold_by_mirrors(pts, axes)
        for (x, y), (u, v) in zip(pts, f):
            self.assertAlmostEqual(math.hypot(x, y), math.hypot(u, v))
            t = (math.degrees(math.atan2(v, u)) - deg0) % 360.0
            self.assertTrue(t <= 22.5 + 1e-9 or t >= 360.0 - 1e-9)
        reps, mult = F.merge(f, 1.0)
        self.assertEqual((len(reps), mult), (1, [8]))

    def test_multiplicities_sum_to_the_point_count(self):
        pts = [[1.0, 0.2], [1.0, -0.2], [-0.5, 0.7], [-0.5, -0.7], [0.3, 0.0]]
        reps, mult, _ = F.folded(pts, [0.0], "fold_shape", 40.0)
        self.assertEqual(sorted(mult), [1, 2, 2])
        reps, mult, _ = F.folded(pts, [0.0], "no_fold", 40.0)
        self.assertEqual(len(reps), 5)
        reps, mult, _ = F.folded(pts, [0.0], "wrong_fold", 40.0)
        self.assertEqual(sum(mult), 5)

    def test_spacing_defect(self):
        self.assertEqual(F.spacing_defect([10.0]), 0.0)
        self.assertAlmostEqual(F.spacing_defect([0.0, 45.0, 90.0, 135.0]), 0.0)
        self.assertAlmostEqual(F.spacing_defect([0.0, 50.0]), 40.0)

    def test_folded_score_finds_the_folded_sets_mirror(self):
        pent = P.regular_polygon(5)
        reps, _, scale = F.folded(pent, [0.0], "fold_shape", 0.0)
        r, deg, _ = F.folded_score(reps, scale, [0.0])
        self.assertAlmostEqual(r, 0.0, places=9)
        self.assertAlmostEqual(deg, 72.0, places=6)


class Loop(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = F.run(seed=0, starts=2, level_trials=60, anchor_max_evals=1000, extra_starts=0)

    def test_levels_log_points_multiplicities_axes_and_wedge(self):
        for s in self.out["shapes"]:
            for r in s["starts"]:
                for v in F.VARIANTS:
                    g = r["variants"][v]
                    self.assertIn("min_pair_distance", g["degeneracy"])
                    self.assertIn(g["stop"], ("no new structure", "single point",
                                              "max_levels reached"))
                    for lv in g["levels"]:
                        for k in ("folded_point_count", "multiplicities", "wedge_angle_deg",
                                  "axes_found_before", "axes_found_after", "degeneracy"):
                            self.assertIn(k, lv)
                        self.assertEqual(sum(lv["multiplicities"]), lv["n_points"])
                        self.assertEqual(len(lv["multiplicities"]), lv["folded_point_count"])

    def test_exact_fold_block(self):
        e = self.out["exact_fold"]
        self.assertEqual(e["level_1"]["multiplicities"], [2, 2, 1])
        self.assertEqual(e["level_1"]["mirrors_of_the_folded_set_deg"], ["72"])
        self.assertEqual(e["level_1"]["each_is_a_mirror_of_the_shape"], [True])
        self.assertEqual(e["level_2"]["multiplicities"], [5])
        self.assertEqual(e["level_2"]["stop"], "single point")
        self.assertIn("16-cell", e["not_exact_here"])


class LoggedRuns(unittest.TestCase):
    def test_observed_counts(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = emlib.load_run(emlib.RUNS / "fold_shape_seed0.json")
        self.assertIsNone(rec["ai_interest"])
        self.assertEqual(rec["outputs"]["label"], "fold_shape")
        by = {s["shape"]: s["compare"] for s in rec["outputs"]["shapes"]}
        reach = {n: [c[k]["starts_reaching_h"]
                     for k in ("b_fold_shape", "c_wrong_fold", "d_no_fold")]
                 for n, c in by.items()}
        self.assertEqual(reach, {"5-cell": [1, 0, 1], "16-cell": [3, 0, 3], "24-cell": [0, 0, 0]})
        for c in by.values():
            g = c["b_fold_shape"]
            self.assertEqual(g["starts_reaching_h"], g["starts_reaching_h_excluding_degenerate"])
            for stop, final in zip(g["stops"], g["final_folded_point_counts"]):
                self.assertEqual(stop == "single point", final == 1)

    def test_search_space_run_is_labelled_and_kept(self):
        rec = emlib.load_run(emlib.RUNS / "fold_loop_seed0.json")
        self.assertEqual(rec["outputs"]["label"], "fold_search_space")
        self.assertFalse(rec["outputs"]["tests_the_carriers_method"])
        by = {s["shape"]: s["compare"] for s in rec["outputs"]["shapes"]}
        self.assertEqual(by["16-cell"]["b_fold"]["final_counts"], [2, 4, 4, 8, 4, 8, 4, 8])

    def test_small_run(self):
        small = emlib.load_run(emlib.RUNS / "fold_shape_small_seed0.json")
        full = emlib.load_run(emlib.RUNS / "fold_shape_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(full["outputs"]))
        self.assertEqual(small["outputs"]["bound"]["extra_shapes_not_run"], ["24-cell"])


if __name__ == "__main__":
    unittest.main()
