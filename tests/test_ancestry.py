"""Tests for experiments/ancestry.py. CC0. Co-authored."""
import math
import unittest

import _path  # noqa: F401
import ancestry as A
import complement_fit_v2 as V
import emlib
from toolkit import complement as CP
from toolkit import primitives as P


class Forms(unittest.TestCase):
    def test_vertex_counts_and_one_sphere(self):
        want = {"5-cell": (10, 20), "16-cell": (24, 48), "tesseract": (32, 64),
                "24-cell": (96, 192)}
        for name, par in A.parents().items():
            rect = A.rectified(par["vertices"])
            trunc, t = A.truncated(par["vertices"], par["face"])
            self.assertEqual((len(rect), len(trunc)), want[name], name)
            self.assertLess(A.on_one_sphere(rect), 1e-9)
            self.assertLess(A.on_one_sphere(trunc), 1e-9)
            self.assertTrue(0.0 < t < 0.5)

    def test_truncation_fraction(self):
        self.assertAlmostEqual(A.truncated(P.simplex(4), 3)[1], 1 / 3)
        self.assertAlmostEqual(A.truncated(P.hypercube(4), 4)[1], 1 / (2 + math.sqrt(2)))

    def test_rectified_16_cell_is_the_24_cell(self):
        rect = CP.unit_dirs(A.rectified(P.cross_polytope(4)))
        self.assertTrue(V.same_signature(V.signature(rect),
                                         V.signature(CP.unit_dirs(P.cell24())), tol=1e-6))

    def test_rectified_5_cell_matches_the_catalogue_entry(self):
        rect = CP.unit_dirs(A.rectified(P.simplex(4)))
        self.assertTrue(V.same_signature(V.signature(rect),
                                         V.signature(V.catalogue()["rectified 5-cell"]), tol=1e-6))

    def test_cover(self):
        cube = CP.unit_dirs(P.hypercube(4))
        c = A.cover(cube, P.hypercube(4))
        self.assertTrue(c["sits_in_and_covers"])
        half = A.cover(cube[:8], P.hypercube(4))
        self.assertTrue(half["every_void_has_a_probe_direction"])
        self.assertFalse(half["every_probe_direction_is_on_a_void"])
        self.assertFalse(half["sits_in_and_covers"])
        self.assertFalse(A.cover(cube, P.cross_polytope(4))["sits_in_and_covers"])

    def test_depth_classes(self):
        cl = A.depth_classes([{"depth_deg": 60.0, "direction": [1]},
                              {"depth_deg": 59.6, "direction": [2]},
                              {"depth_deg": 45.0, "direction": [3]}])
        self.assertEqual([len(c["directions"]) for c in cl], [2, 1])


class LoggedRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = emlib.load_run(emlib.RUNS / "ancestry_seed0.json")
        cls.o = cls.rec["outputs"]

    def test_observer_slots_null(self):
        self.assertIsNone(self.rec["ai_interest"])

    def test_observed_summary(self):
        # OBSERVED on the committed run, pinned so a change shows up
        s = self.o["summary"]
        self.assertEqual(s["forms"], 8)
        self.assertEqual(s["deepest_covered_by_parent"], [])
        self.assertEqual(len(s["deepest_covered_by_dual"]), 7)
        self.assertEqual(s["deepest_covered_by_parent_and_dual_together"], ["rectified 16-cell"])
        self.assertEqual(s["deepest_covered_by_none_of_the_three"], [])
        self.assertEqual(s["some_class_covered_by_parent"], [])
        self.assertEqual(set(map(tuple, s["plain_deepest_covered_by"].values())), {("dual",)})

    def test_every_class_logs_the_three_probes_and_the_bound(self):
        self.assertEqual(len(self.o["baselines"]), 12)
        for r in self.o["baselines"]:
            self.assertIn("found_is_a_lower_bound", r["bound"])
            for c in r["void_classes"]:
                for k in ("parent", "dual", "parent_and_dual_together"):
                    self.assertIn("sits_in_and_covers", c[k])

    def test_small_run_and_doc(self):
        small = emlib.load_run(emlib.RUNS / "ancestry_small_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(self.o))
        text = (emlib.ROOT / "docs" / "complement.md").read_text()
        self.assertIn("void cover = dual OR parent", text)
        self.assertIn("stays PROPOSED", text)


if __name__ == "__main__":
    unittest.main()
