"""Tests for probes.degeneracy, degeneracy_check and constrained_generate. CC0. Co-authored."""
import json
import tempfile
import unittest
from pathlib import Path

import _path  # noqa: F401
import constrained_generate as K
import degeneracy_check as D
import emlib
from toolkit import primitives as P
from toolkit import probes as Q


class Degeneracy(unittest.TestCase):
    def test_clean_set(self):
        d = Q.degeneracy(P.regular_polygon(5))
        self.assertFalse(d["tags"]["1e-06"]["degenerate"])
        self.assertAlmostEqual(d["collinearity"], 1.0, places=9)
        self.assertEqual(d["tags"]["1e-06"]["n_distinct"], 5)

    def test_coincident_points(self):
        d = Q.degeneracy([[1.0, 0.0], [1.0, 0.0], [-1.0, 0.5], [0.0, -1.0]])
        self.assertEqual(d["min_pair_distance"], 0.0)
        self.assertEqual(d["closest_pair"], [0, 1])
        self.assertTrue(d["tags"]["1e-06"]["coincident"])
        self.assertEqual(d["tags"]["1e-06"]["n_distinct"], 3)

    def test_collinear_points(self):
        d = Q.degeneracy([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [-3.0, -3.0]])
        self.assertAlmostEqual(d["collinearity"], 0.0, places=7)
        self.assertTrue(d["tags"]["0.001"]["collinear"])
        self.assertFalse(d["tags"]["0.001"]["coincident"])

    def test_near_coincident_is_tagged_only_at_the_looser_tolerance(self):
        d = Q.degeneracy([[1.0, 0.0], [1.0, 1e-4], [-1.0, 0.5], [0.0, -1.0]])
        self.assertFalse(d["tags"]["1e-06"]["degenerate"])
        self.assertTrue(d["tags"]["0.001"]["degenerate"])


class DegeneracyCheck(unittest.TestCase):
    def test_absent_sources_are_listed_and_recount_excludes_degenerate(self):
        with tempfile.TemporaryDirectory() as d:
            sq = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]]
            dbl = sq + sq                                   # 8 points, 4 distinct
            rec = {"outputs": {"shapes": [{"shape": "16-cell", "coxeter_plane_mirror_count": 4,
                                           "starts": [
                {"start_index": 0, "final_projected": dbl,
                 "final_counts": [{"rel_tol": 1e-6, "n_mirror_axes": 4}]},
                {"start_index": 1, "final_projected": P.regular_polygon(8),
                 "final_counts": [{"rel_tol": 1e-6, "n_mirror_axes": 4}]}]}]}}
            (Path(d) / "mirror_pursuit_seed0.json").write_text(json.dumps(rec))
            out = D.run(runs_dir=d)
        self.assertIn("coxeter_plane_seed0.json", out["sources_absent"])
        self.assertEqual(out["n_rows"], 2)
        g = out["recount"][0]
        self.assertEqual((g["reaching_h"], g["reaching_h_degenerate"],
                          g["reaching_h_excluding_degenerate"]), (2, 1, 1))
        self.assertEqual(g["final_counts_excluding_degenerate"], [None, 4])

    def test_observed_on_the_logged_runs(self):
        # OBSERVED on the committed runs, pinned so a change shows up
        out = D.run()
        self.assertEqual(out["sources_absent"], [])
        self.assertTrue(all(abs(r["collinearity"] - 1.0) < 1e-6 for r in out["rows"]))
        for g in out["recount"]:
            self.assertEqual(g["reaching_h"], g["reaching_h_excluding_degenerate"], g["method"])
        frame0 = [r for r in out["rows"] if r["method"] == "tour_frame_0"]
        self.assertEqual(len(frame0), 2)
        self.assertTrue(all(r["degenerate"] for r in frame0))
        cox = [r for r in out["rows"] if r["method"] == "coxeter_plane_by_construction"]
        self.assertFalse(any(r["degenerate"] for r in cox))


class PullBack(unittest.TestCase):
    def test_solve(self):
        x = K.solve([[2.0, 1.0], [1.0, 3.0]], [3.0, 5.0])
        self.assertAlmostEqual(x[0], 0.8)
        self.assertAlmostEqual(x[1], 1.4)
        self.assertIsNone(K.solve([[0.0, 0.0], [0.0, 0.0]], [1.0, 1.0]))

    def test_constraint_is_zero_on_an_exact_mirror(self):
        pts = [[1.0, 0.5], [1.0, -0.5], [-0.7, 0.2], [-0.7, -0.2]]
        pair = K.pairing(pts, 0.0)
        self.assertEqual(pair, [1, 0, 3, 2])
        self.assertAlmostEqual(max(abs(v) for v in K.constraint(pts, 0.0, pair)), 0.0)

    def test_pull_back_restores_a_disturbed_mirror(self):
        import mirror_pursuit as U
        verts = P.cross_polytope(4)
        planes = P.all_planes(4)
        anc = U.descend(verts, [10.0, 50.0, 80.0, 200.0, 310.0, 33.0], max_evals=1500)
        a = anc["final_reading"]["axis_deg"]
        pair = K.pairing(K.project(verts, planes, anc["final_angles_deg"]), a)
        moved = [x + 0.5 for x in anc["final_angles_deg"]]
        before = Q.mirror_residual_at(K.project(verts, planes, moved), [a])[0]
        back, _ = K.pull_back(verts, planes, moved, a, pair)
        after = Q.mirror_residual_at(K.project(verts, planes, back), [a])[0]
        self.assertGreater(before, 1e-4)
        self.assertLess(after, K.ANCHOR_TOL)


class ConstrainedGenerate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = K.run(seed=0, starts=1, max_trials=40, anchor_max_evals=800,
                        bound_check_trials=None)

    def test_anchor_is_held_on_every_logged_position(self):
        for s in self.out["shapes"]:
            for g in s["starts"][0]["rules"].values():
                self.assertTrue(g["anchor_held"])
                self.assertLess(g["max_anchor_residual_on_path"], K.ANCHOR_TOL)
                sc = [p["score"] for p in g["path"]]
                self.assertEqual(sc, sorted(sc, reverse=True))
                self.assertEqual(len(g["final_residual_per_axis"]), g["m"])
                self.assertIn("degeneracy", g)

    def test_bound_check_off_is_void(self):
        for s in self.out["shapes"]:
            self.assertIsNone(s["compare"]["constrained_bound_check"])


class LoggedConstrainedRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = emlib.load_run(emlib.RUNS / "constrained_generate_seed0.json")
        cls.by = {s["shape"]: s["compare"] for s in cls.rec["outputs"]["shapes"]}

    def test_observer_slots_null_and_same_starts(self):
        self.assertIsNone(self.rec["human_interest"])
        self.assertIsNone(self.rec["ai_interest"])
        for c in self.by.values():
            self.assertTrue(c["starts_equal_anchor_generate_starts"])

    def test_observed_counts(self):
        # OBSERVED on the committed run, pinned so a change shows up
        def reach(key):
            return {n: {k: c[key][k]["starts_reaching_h"]
                        for k in ("h_minus_1", "h", "h_plus_1")} for n, c in self.by.items()}
        self.assertEqual(reach("constrained"),
                         {"5-cell": {"h_minus_1": 0, "h": 0, "h_plus_1": 0},
                          "16-cell": {"h_minus_1": 0, "h": 1, "h_plus_1": 0}})
        self.assertEqual(reach("constrained_bound_check"),
                         {"5-cell": {"h_minus_1": 0, "h": 1, "h_plus_1": 0},
                          "16-cell": {"h_minus_1": 0, "h": 5, "h_plus_1": 0}})
        for c in self.by.values():
            for key in ("constrained", "constrained_bound_check"):
                for g in c[key].values():
                    self.assertEqual(g["starts_holding_anchor"], 8)

    def test_small_run_is_a_separate_file_with_the_same_output_keys(self):
        small = emlib.load_run(emlib.RUNS / "constrained_generate_small_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(self.rec["outputs"]))
        self.assertEqual(set(small["outputs"]["shapes"][0]), set(self.rec["outputs"]["shapes"][0]))


if __name__ == "__main__":
    unittest.main()
