"""Tests for toolkit/hull4.py and experiments/sphere_compounds.py. CC0. Co-authored."""
import re
import unittest

import _path  # noqa: F401
import complement_fit_v2 as V
import emlib
import sphere_compounds as S
from toolkit import complement as CP
from toolkit import hull4 as H4
from toolkit import primitives as P


class Hull(unittest.TestCase):
    def test_known_cell_counts(self):
        want = {"5-cell": (P.simplex(4), 5, {"4": 5}), "16-cell": (P.cross_polytope(4), 16, {"4": 16}),
                "tesseract": (P.hypercube(4), 8, {"8": 8}), "24-cell": (P.cell24(), 24, {"6": 24})}
        for name, (verts, cells, sizes) in want.items():
            s = H4.summary(CP.unit_dirs(verts))
            self.assertEqual((s["n_facets"], s["facets_by_vertex_count"]), (cells, sizes), name)

    def test_closed_detects_a_missing_cell(self):
        fs = H4.facets(CP.unit_dirs(P.cross_polytope(4)))
        self.assertTrue(H4.closed(fs))
        self.assertFalse(H4.closed(fs[:-1]))
        self.assertIsNone(H4.closed(H4.facets(CP.unit_dirs(P.hypercube(4)))))

    def test_edge_limit_can_miss_facets_and_closed_says_so(self):
        d = CP.unit_dirs(P.cross_polytope(4))            # edges are at 90 degrees
        self.assertEqual(H4.facets(d, max_edge_deg=60.0), [])


class Compounds(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = S.run(seed=0, starts=3, jiggles=500, polish_moves=3000)

    def test_two_16_cells_built_exactly(self):
        b = self.out["two_16_cells"]["built_outright"]
        self.assertTrue(b["each_set_is_a_16_cell"])
        self.assertEqual(b["dots_between_the_sets"], ["-1/2", "1/2"])
        self.assertTrue(b["union_is_the_tesseract_vertex_set"])
        self.assertEqual(b["hull"]["facets_by_vertex_count"], {"8": 8})
        self.assertTrue(b["against_catalogue"][0]["match_at_0.5_deg"])

    def test_search_never_exceeds_the_built_separation(self):
        s = self.out["two_16_cells"]["search"]
        self.assertLessEqual(s["separation_after_polish_deg"], 60.0 + 1e-6)
        self.assertGreaterEqual(s["separation_after_polish_deg"],
                                s["separation_after_search_deg"] - 1e-9)

    def test_catalogue_additions_are_marked_carried_and_computed(self):
        c = self.out["catalogue_additions"]
        for entry in c.values():
            self.assertFalse(entry["name_checked_against_a_source"])
            self.assertTrue(entry["computed"]["cells_close_up"])
        d = c["disphenoidal 288-cell"]["computed"]
        self.assertEqual((d["vertices"], d["edges"], d["faces_by_euler"], d["cells"]),
                         (48, 336, 576, 288))
        self.assertTrue(d["cells_are_tetragonal_disphenoids"])
        self.assertTrue(d["every_vertex_sees_the_same_angles"])
        t = c["compound of two 5-cells"]["computed"]
        self.assertEqual((t["vertices"], t["edges"], t["cells"]), (10, 40, 30))

    def test_two_24_cells(self):
        w = self.out["two_24_cells"]
        self.assertTrue(w["built_outright"]["dual_is_a_24_cell"])
        self.assertAlmostEqual(w["built_outright"]["separation_deg"], 45.0, places=6)
        x = w["searched_in_complement_fit_v2"]
        self.assertTrue(x["same_as_the_built_placement_at_0.5_deg"])
        self.assertEqual(x["hull"]["n_facets"], 288)


class LoggedRun(unittest.TestCase):
    def test_observed(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = emlib.load_run(emlib.RUNS / "sphere_compounds_seed0.json")
        self.assertIsNone(rec["ai_interest"])
        s = rec["outputs"]["two_16_cells"]["search"]
        self.assertTrue(s["matches_tesseract_at_0.5_deg"])
        self.assertLess(s["short_of_target_deg"], 1e-4)
        fit = emlib.load_run(emlib.RUNS / "complement_fit_v2_seed0.json")
        row = next(r for r in fit["outputs"]["sphere"]["pairs"]
                   if r["baseline"] == "24-cell" and r["probe"] == "24-cell")
        self.assertEqual(row["matches_in_catalogue"], ["24-cell + its dual (as built)"])
        small = emlib.load_run(emlib.RUNS / "sphere_compounds_small_seed0.json")
        self.assertEqual(set(small["outputs"]), set(rec["outputs"]))

    def test_doc_run_pointers(self):
        text = (emlib.ROOT / "docs" / "complement.md").read_text()
        self.assertIn("weak at low N", text)
        for name in re.findall(r"`([a-z_0-9]+)`", text):
            if name in ("collapse_gap", "complement_fit", "complement_fit_v2", "sphere_compounds"):
                self.assertTrue((emlib.RUNS / f"{name}_seed0.json.gz").exists(), name)
        self.assertEqual(len(V.catalogue()), 8)


if __name__ == "__main__":
    unittest.main()
