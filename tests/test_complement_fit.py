"""Tests for experiments/complement_fit.py. CC0. Co-authored."""
import unittest

import _path  # noqa: F401
import complement_fit as F
import emlib
from toolkit import complement as CP
from toolkit import primitives as P


class Pieces(unittest.TestCase):
    def setUp(self):
        self.n = 840
        self.radii, _ = CP.polar_grid(F.N_RINGS, self.n)

    def field(self, pts, scale=1.0):
        pat = CP.turned_to_zero(CP.centred_unit(pts))
        return CP.distance_field([[scale * x, scale * y] for x, y in pat], self.radii, self.n)

    def test_weights(self):
        none = F.weights(self.radii, self.n, 20, "none")
        zero = F.weights(self.radii, self.n, 20, "uniform_zero")
        grad = F.weights(self.radii, self.n, 20, "gradient")
        self.assertEqual(none, zero)
        self.assertTrue(all(x == 1.0 for row in none for x in row))
        self.assertLess(grad[-1][0], 1.0)                    # x > 0: larger local unit
        self.assertGreater(grad[-1][len(grad[-1]) // 2], 1.0)   # x < 0: smaller local unit

    def test_placing_a_probe_never_widens_the_gap_and_finds_the_offset_turn(self):
        base = F.sample(self.field(P.regular_polygon(5)), self.n, 20)
        w = F.weights(self.radii, self.n, 20, "none")
        g0 = F.gap(base, w)
        g, scale, turn, new = F.place(base, {1.0: self.field(P.regular_polygon(5))},
                                      self.n, 20, 7, w)
        self.assertLess(g, g0)
        self.assertAlmostEqual(turn % 72.0, 36.0)            # a pentagon between the points
        self.assertTrue(all(a <= b + 1e-15 for ra, rb in zip(new, base) for a, b in zip(ra, rb)))
        self.assertAlmostEqual(F.gap(new, w), g)

    def test_ranking(self):
        rows = F.ranking([{"name": "b", "residual_gap": 0.2}, {"name": "a", "residual_gap": 0.1}])
        self.assertEqual([(r["name"], r["rank"]) for r in rows], [("a", 1), ("b", 2)])


class Run(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = F.run(seed=0, n_theta=840, stride=40, turn_step=14, baselines=("5-cell",),
                        combo_probes=("3-gon", "5-gon", "8-gon"), sphere_starts=20,
                        sphere_steps=150)

    def test_singles_cover_every_probe_and_never_exceed_the_baseline_gap(self):
        b = self.out["plane"][0]
        for env in F.ENVIRONMENTS:
            e = b["environments"][env]
            self.assertEqual(len(e["singles"]), 12)
            self.assertIn("self", e["ranking"])
            for r in e["singles"]:
                self.assertLessEqual(r["residual_gap"], e["baseline_gap"] + 1e-12)
                self.assertAlmostEqual(r["gap_closed"], e["baseline_gap"] - r["residual_gap"])
            self.assertEqual([r["rank"] for r in e["singles"]], list(range(1, 13)))

    def test_control_environment_keeps_the_order(self):
        b = self.out["plane"][0]
        self.assertTrue(b["shift_none_to_uniform_zero"]["same_order"])
        self.assertEqual(b["environments"]["none"]["ranking"],
                         b["environments"]["uniform_zero"]["ranking"])
        self.assertNotIn("pairs", b["environments"]["uniform_zero"])

    def test_combinations(self):
        e = self.out["plane"][0]["environments"]["none"]
        self.assertEqual(len(e["pairs"]), 3)
        self.assertEqual(len(e["triples"]), 1)
        singles = {r["shape"]: r["residual_gap"] for r in e["singles"]}
        for r in e["pairs"] + e["triples"]:
            self.assertEqual(len(r["placement"]), len(r["probes"]))
            self.assertLessEqual(r["residual_gap"], min(singles[p] for p in r["probes"]) + 1e-12)
            gaps = [s["gap_after"] for s in r["placement"]]
            self.assertEqual(gaps, sorted(gaps, reverse=True))

    def test_shift_log(self):
        s = self.out["plane"][0]["shift_none_to_gradient"]
        for k in ("singles", "singles_per_point", "pairs", "triples"):
            self.assertEqual(s[k]["n_changed_rank"],
                             sum(1 for m in s[k]["moves"] if m["moved_by"]))

    def test_sphere(self):
        b = self.out["sphere"]["baselines"][0]
        self.assertEqual(len(b["probes"]), 7)
        for r in b["probes"]:
            self.assertLessEqual(r["residual_gap"], b["baseline_deepest_void_deg"] + 1.0)


class LoggedRun(unittest.TestCase):
    def test_observed(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = emlib.load_run(emlib.RUNS / "complement_fit_seed0.json")
        self.assertIsNone(rec["ai_interest"])
        by = {b["baseline"]: b for b in rec["outputs"]["plane"]}
        self.assertEqual(set(by), {"5-cell", "16-cell"})
        for b in by.values():
            self.assertTrue(b["shift_none_to_uniform_zero"]["same_order"])
            self.assertEqual(b["environments"]["none"]["ranking"][0], "24-cell")
            self.assertEqual(b["environments"]["none"]["ranking_by_gap_closed_per_point"][0],
                             "3-gon")
        self.assertEqual(by["5-cell"]["shift_none_to_gradient"]["singles"]["n_changed_rank"], 0)
        self.assertEqual(by["16-cell"]["shift_none_to_gradient"]["singles"]["n_changed_rank"], 3)
        self.assertEqual(by["16-cell"]["shift_none_to_gradient"]["singles"]["largest_move"], 2)
        sph = {b["baseline"]: b["ranking"] for b in rec["outputs"]["sphere"]["baselines"]}
        self.assertEqual(sph["5-cell"][-1], "self")
        small = emlib.load_run(emlib.RUNS / "complement_fit_small_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(rec["outputs"]))


if __name__ == "__main__":
    unittest.main()
