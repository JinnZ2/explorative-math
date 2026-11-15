"""Tests for experiments/complement_fit_v2.py. CC0. Co-authored."""
import math
import random
import unittest

import _path  # noqa: F401
import complement_fit as F
import complement_fit_v2 as V
import emlib
from toolkit import complement as CP
from toolkit import primitives as P


class Budget(unittest.TestCase):
    def test_greedy_picks_n_points_and_never_widens_the_gap(self):
        radii, _ = CP.polar_grid(F.N_RINGS, 840)
        samples = V.sample_points(radii, 840, 40)
        w = F.weights(radii, 840, 40, "none")
        base = V.field_at(CP.turned_to_zero(CP.centred_unit(P.regular_polygon(5))), samples)
        cand = V.placed(P.regular_polygon(8), 0.5, 15.0)
        g3, idx3, _ = V.greedy_points(base, cand, 3, samples, w)
        g5, idx5, _ = V.greedy_points(base, cand, 5, samples, w)
        self.assertEqual((len(idx3), len(idx5)), (3, 5))
        self.assertEqual(idx5[:3], idx3)                       # greedy: a prefix
        self.assertEqual(len(set(idx5)), 5)
        self.assertLess(g5, g3)
        self.assertLess(g3, F.gap(base, w))

    def test_placed_matches_the_shifted_field(self):
        # a probe placed by V.placed has the field F.place read at that turn
        n = 840
        radii, _ = CP.polar_grid(F.N_RINGS, n)
        pat = CP.turned_to_zero(CP.centred_unit(P.regular_polygon(3)))
        fine = CP.distance_field(pat, radii, n)
        shift = 70                                              # 30 degrees
        a = F.sample(fine, n, 40, shift)
        b = V.field_at(V.placed(pat, 1.0, 360.0 * shift / n), V.sample_points(radii, n, 40))
        self.assertLess(max(abs(x - y) for ra, rb in zip(a, b) for x, y in zip(ra, rb)), 1e-9)


class Response(unittest.TestCase):
    def test_points_move_with_the_sheet_and_the_centre_barely(self):
        move, bound = V.response_map(4, V.SHEET_ALPHA, V.SHEET_SLOPE, 600)
        self.assertGreater(bound["largest_displacement_plane_units"], 1e-3)
        right, left = move([0.9, 0.0]), move([-0.9, 0.0])
        self.assertNotAlmostEqual(right[0] - 0.9, left[0] + 0.9, places=4)
        self.assertEqual(len(move([0.2, 0.3])), 2)

    def test_zero_field_leaves_points_where_they_are(self):
        move, _ = V.response_map(3, 0.0, 1.0, 200)
        p = move([0.4, -0.3])
        self.assertAlmostEqual(p[0], 0.4, places=9)
        self.assertAlmostEqual(p[1], -0.3, places=9)


class Sphere(unittest.TestCase):
    def test_signature_is_unchanged_by_rotation(self):
        d = CP.unit_dirs(P.hypercube(4))
        r = V.rotate4(d, [13.0, 47.0, 211.0, 5.0, 300.0, 88.0])
        self.assertTrue(V.same_signature(V.signature(d), V.signature(r), tol=1e-6))
        self.assertFalse(V.same_signature(V.signature(d),
                                          V.signature(CP.unit_dirs(P.cross_polytope(4)))))

    def test_catalogue(self):
        cat = V.catalogue()
        self.assertEqual({k: len(v) for k, v in cat.items()},
                         {"5-cell": 5, "16-cell": 8, "tesseract": 16, "24-cell": 24,
                          "rectified 5-cell": 10, "5-cell + its negative (as built)": 10,
                          "16-cell + tesseract (as built)": 24,
                          "24-cell + its dual (as built)": 48})
        s = {k: V.signature(v) for k, v in cat.items()}
        self.assertTrue(V.same_signature(s["24-cell"], s["16-cell + tesseract (as built)"]))
        self.assertFalse(V.same_signature(s["rectified 5-cell"],
                                          s["5-cell + its negative (as built)"]))

    def test_search_is_a_lower_bound_on_the_negated_simplex(self):
        # the largest separation of a 5-cell from a copy of itself is acos(1/4)
        # = 75.52 degrees (the negated copy). A short search gets near it and
        # never above it; how near depends on the effort, which is why the
        # effort is logged with every pair.
        sim = CP.unit_dirs(P.simplex(4))
        top = math.degrees(math.acos(0.25))
        s = V.search_rotation(sim, sim, random.Random(1), starts=4, jiggles=800)
        self.assertGreater(s["separation_deg"], 70.0)
        self.assertLessEqual(s["separation_deg"], top + 1e-6)
        neg = [[-x for x in v] for v in sim]
        self.assertAlmostEqual(V.separation(sim, neg), top, places=6)

    def test_search_effort_shrinks_for_large_pairs(self):
        self.assertEqual(V.search_effort(5, 5, 12), (12, V.ROT_JIGGLES))
        self.assertEqual(V.search_effort(24, 24, 12), (4, V.ROT_JIGGLES // 4))


class LoggedRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = emlib.load_run(emlib.RUNS / "complement_fit_v2_seed0.json")
        cls.o = cls.rec["outputs"]

    def test_observer_slots_null(self):
        self.assertIsNone(self.rec["ai_interest"])
        self.assertIsNone(self.rec["human_interest"])

    def test_budget(self):
        for x in self.o["equal_budget"]:
            for bu in x["budgets"]:
                for r in bu["singles"]:
                    self.assertEqual(r["n_points_used"], bu["budget"])
                    self.assertEqual(len(r["points_chosen"]), bu["budget"])
                for r in bu["pairs"]:
                    self.assertEqual(sum(r["points_from_each"]), bu["budget"])
            # OBSERVED, pinned: at N = 8 only shapes with 8 or more points are ranked
            n8 = next(bu for bu in x["budgets"] if bu["budget"] == 8)
            self.assertEqual(len(n8["singles"]), 6)

    def test_response(self):
        # OBSERVED, pinned so a change shows up
        got = {x["baseline"]: (x["shift_singles"]["n_changed_rank"],
                               x["shift_singles"]["largest_move"],
                               x["shift_singles"]["best_after"])
               for x in self.o["field_response"]["baselines"]}
        self.assertEqual(got, {"5-cell": (4, 1, "24-cell"), "16-cell": (4, 1, "24-cell")})

    def test_sphere(self):
        rows = {(r["baseline"], r["probe"]): r for r in self.o["sphere"]["pairs"]}
        self.assertEqual(rows[("5-cell", "5-cell")]["matches_in_catalogue"],
                         ["5-cell + its negative (as built)"])
        self.assertIn("24-cell", rows[("16-cell", "tesseract")]["matches_in_catalogue"])
        self.assertEqual(self.o["sphere"]["catalogue_entries_that_match_each_other"],
                         [["24-cell", "16-cell + tesseract (as built)"]])
        for r in rows.values():
            self.assertGreaterEqual(r["separation_after_search_deg"],
                                    r["separation_as_built_deg"] - 1e-9)
            self.assertIn("rotation_starts", r)

    def test_small_run(self):
        small = emlib.load_run(emlib.RUNS / "complement_fit_v2_small_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(self.o))


if __name__ == "__main__":
    unittest.main()
