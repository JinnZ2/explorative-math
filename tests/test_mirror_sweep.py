"""Tests for mirror_sweep, coxeter_plane and primitives.coxeter_plane. CC0. Co-authored."""
import json
import math
import unittest

import _path  # noqa: F401
import coxeter_plane as K
import mirror_sweep as M
from toolkit import primitives as P


class CoxeterPlanePrimitive(unittest.TestCase):
    def test_reflect(self):
        self.assertEqual(P.reflect([1.0, 2.0], [1.0, 0.0]), [-1.0, 2.0])

    def test_basis_is_orthonormal_plane(self):
        for roots, h in ((K.A4, 5), (K.B4, 8), (K.F4, 12)):
            u, v = P.coxeter_plane(roots, h)
            self.assertAlmostEqual(P.dot(u, u), 1.0)
            self.assertAlmostEqual(P.dot(v, v), 1.0)
            self.assertAlmostEqual(P.dot(u, v), 0.0)

    def test_wrong_coxeter_number_is_refused(self):
        with self.assertRaises(ValueError):
            P.coxeter_plane(K.B4, 5)

    def test_coxeter_element_rotates_the_plane_by_360_over_h(self):
        u, v = P.coxeter_plane(K.B4, 8)
        w = u
        for r in K.B4:
            w = P.reflect(w, r)
        self.assertAlmostEqual(P.dot(w, u), math.cos(2 * math.pi / 8))
        self.assertAlmostEqual(abs(P.dot(w, v)), math.sin(2 * math.pi / 8))


class CoxeterPlaneRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = K.run()

    def test_observed_mirror_count_equals_h(self):
        # OBSERVED on this run, pinned so a change shows up
        got = [(c["name"], c["readings"][0]["n_mirror_axes"]) for c in self.out["cases"]]
        self.assertEqual(got, [("5-cell", 5), ("16-cell", 8), ("tesseract", 8), ("24-cell", 12)])
        for c in self.out["cases"]:
            self.assertTrue(c["n_mirror_axes_equals_h"])
            self.assertIn(c["coxeter_number_h"], c["readings"][0]["rotational_symmetry"])

    def test_every_reading_logs_its_bounds(self):
        for c in self.out["cases"]:
            for r in c["readings"]:
                self.assertEqual(r["rotational_symmetry_bound"]["k_tested"], [2, K.K_MAX])
                self.assertEqual(r["mirror_bound"]["rel_tol"], r["rel_tol"])
                self.assertEqual(r["mirror_grid_bound"]["n_axes"], K.N_AXES)

    def test_missing_sweep_reference_is_void(self):
        out = K.run(sweep_path="/nonexistent/sweep.json")
        for c in out["cases"]:
            self.assertIsNone(c["sweep_reference"])
            self.assertIsNone(c["more_axes_than_sweep_max"])
        json.dumps(out)


class MirrorSweep(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = M.run(frames=36)

    def test_observed_frame_zero_is_identity(self):
        for f in self.out["frame_zero"]:
            self.assertTrue(f["all_rotation_angles_zero"])
            self.assertTrue(f["projected_equals_first_two_coordinates"])

    def test_counts_cover_every_frame_and_log_the_bound(self):
        for s in self.out["sweep"]:
            self.assertEqual(s["frames"], 36)
            for k in M.SETS:
                self.assertEqual(len(s["axis_count_per_frame"][k]), 36)
            self.assertEqual(s["bound"]["rel_tol"], M.REL_TOL)

    def test_observed_only_frame_zero_has_a_mirror(self):
        # OBSERVED at 36 frames, pinned; the logged run reads 360
        for s in self.out["sweep"]:
            self.assertEqual([m["frame_index"] for m in s["mirror_frames"]], [0])
            self.assertEqual(s["frames_with_any_mirror"]["a_union"], 0)

    def test_logged_run_reads_360_frames(self):
        rec = json.loads((M.emlib.RUNS / "mirror_sweep_seed0.json").read_text())
        for s in rec["outputs"]["sweep"]:
            self.assertEqual(s["frames"], 360)
            self.assertEqual(s["frames_with_any_mirror"]["b_union"], 1)


if __name__ == "__main__":
    unittest.main()
