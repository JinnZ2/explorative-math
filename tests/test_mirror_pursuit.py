"""Tests for probes.mirror_residual, the axis merge, the residual sweep and
mirror_pursuit. CC0. Co-authored."""
import json
import math
import unittest

import _path  # noqa: F401
import mirror_pursuit as U
import mirror_sweep as M
from toolkit import primitives as P
from toolkit import probes as Q


class MirrorResidual(unittest.TestCase):
    def test_zero_on_an_exact_mirror(self):
        r = Q.mirror_residual(P.regular_polygon(5, phase=0.1234))
        self.assertAlmostEqual(r["residual"], 0.0, places=12)
        self.assertEqual(r["n_points"], 5)
        self.assertIn("matching", r)

    def test_positive_without_a_mirror_and_shrinks_on_approach(self):
        base = [[1.0, 0.0], [-1.0, 0.0], [0.3, 1.0], [0.3, -1.0]]     # mirror: the x axis
        far = [[1.0, 0.0], [-1.0, 0.0], [0.3, 1.0], [0.9, -1.4]]
        near = [[1.0, 0.0], [-1.0, 0.0], [0.3, 1.0], [0.31, -1.01]]
        self.assertAlmostEqual(Q.mirror_residual(base)["residual"], 0.0, places=12)
        rf, rn = Q.mirror_residual(far)["residual"], Q.mirror_residual(near)["residual"]
        self.assertGreater(rn, 0.0)
        self.assertGreater(rf, rn)

    def test_scale_free(self):
        pts = [[1.0, 0.0], [-1.0, 0.2], [0.3, 1.0], [0.9, -1.4]]
        big = [[100 * x, 100 * y] for x, y in pts]
        self.assertAlmostEqual(Q.mirror_residual(pts)["residual"],
                               Q.mirror_residual(big)["residual"])


class AxisMerge(unittest.TestCase):
    def test_one_near_exact_mirror_is_one_axis(self):
        # mirror about the x axis, each partner off by about 1e-10
        pts = [[1.0, 0.5], [1.0 + 1e-10, -0.5], [-0.7, 0.2], [-0.7, -0.2 - 1e-10],
               [0.2, 1.3], [0.2 - 1e-10, -1.3]]
        axes = Q.mirror_symmetry_from_pairs(pts)
        self.assertEqual(len(axes), 1)

    def test_exact_counts_unchanged(self):
        for k in (3, 4, 5, 8, 12):
            self.assertEqual(len(Q.mirror_symmetry_from_pairs(P.regular_polygon(k, phase=0.3))), k)

    def test_axis_near_zero_and_180_is_one_axis(self):
        pts = [[1.0, 0.5], [1.0, -0.5 - 1e-10], [-0.7, 0.2], [-0.7 + 1e-10, -0.2]]
        self.assertEqual(len(Q.mirror_symmetry_from_pairs(pts)), 1)


class ResidualSweep(unittest.TestCase):
    def test_residual_logged_per_frame_and_zero_at_frame_zero(self):
        out = M.run(frames=36)
        for s in out["sweep"]:
            for k in M.SETS:
                self.assertEqual(len(s["residual_per_frame"][k]), 36)
            self.assertAlmostEqual(s["residual_per_frame"]["shape"][0], 0.0, places=12)
            self.assertNotEqual(s["residual_min_excluding_frame_0"]["shape"][0], 0)
            self.assertGreater(s["residual_min_excluding_frame_0"]["shape"][1], 0.0)

    def test_local_minima(self):
        self.assertEqual(M.local_minima([3.0, 1.0, 2.0, 0.5, 4.0]), [[3, 0.5], [1, 1.0]])
        self.assertEqual(M.local_minima([0.0, 1.0, 2.0]), [[0, 0.0]])

    def test_observed_dips_in_the_logged_run(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = json.loads((M.emlib.RUNS / "mirror_sweep_seed0.json").read_text())
        got = {s["shape"]: s["residual_min_excluding_frame_0"]["shape"][0]
               for s in rec["outputs"]["sweep"]}
        self.assertEqual(got, {"simplex": 137, "crosspolytope": 1})


class MirrorPursuit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = U.run(seed=0, starts=2)

    def test_path_never_rises_and_is_logged_from_the_start(self):
        for s in self.out["shapes"]:
            for r in s["starts"]:
                res = [p["residual"] for p in r["path"]]
                self.assertEqual(res, sorted(res, reverse=True))
                self.assertEqual(r["path"][0]["angles_deg"], r["start_angles_deg"])
                self.assertEqual(r["path"][-1]["residual"], r["final_residual"])
                self.assertLessEqual(r["final_residual"], r["start_residual"])
                self.assertIn(r["stop"], ("step below step_min", "max_evals reached",
                                          "residual below residual_stop"))

    def test_deterministic(self):
        again = U.run(seed=0, starts=2)
        self.assertEqual(json.dumps(again), json.dumps(self.out))

    def test_observed_counts_in_the_logged_run(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = json.loads((U.emlib.RUNS / "mirror_pursuit_seed0.json").read_text())
        by = {s["shape"]: s for s in rec["outputs"]["shapes"]}
        self.assertEqual(by["5-cell"]["final_mirror_counts"], [1] * 8)
        self.assertEqual(by["16-cell"]["final_mirror_counts"], [2, 2, 2, 2, 0, 2, 2, 2])
        self.assertEqual(by["5-cell"]["starts_reaching_coxeter_count"], 0)
        self.assertEqual(by["16-cell"]["starts_reaching_coxeter_count"], 0)
        # no point set is given more mirror axes than it has points
        for s in by.values():
            for r in s["starts"]:
                self.assertLessEqual(r["final_counts"][0]["n_mirror_axes"], s["n_vertices"])

    def test_missing_coxeter_reference_is_void(self):
        out = U.run(seed=0, starts=1, max_evals=50, coxeter_path="/nonexistent.json")
        for s in out["shapes"]:
            self.assertIsNone(s["coxeter_plane_mirror_count"])
            self.assertIsNone(s["starts_reaching_coxeter_count"])
        self.assertTrue(math.isfinite(out["shapes"][0]["starts"][0]["final_residual"]))


if __name__ == "__main__":
    unittest.main()
