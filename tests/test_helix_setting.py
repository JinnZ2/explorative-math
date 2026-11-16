"""Tests for experiments/helix_setting.py. CC0. Co-authored."""
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
import helix_setting as H  # noqa: E402
import spherical_helix as SH  # noqa: E402
from toolkit import primitives as P  # noqa: E402

OUT = H.run(0)


def dihedral_from_normals(vs):
    """Angle between two facets of the spherical simplex on vs, from normals."""
    def normal(skip):
        rest = [v for i, v in enumerate(vs) if i != skip]
        basis = P.gram_schmidt([list(v) for v in rest], tol=1e-12)
        w = list(vs[skip])
        for b in basis:
            w = P.sub(w, P.scale(b, P.dot(w, b)))
        return P.scale(w, 1.0 / P.norm(w))      # inward normal of the facet opposite `skip`
    a, b = normal(0), normal(1)
    return math.degrees(math.acos(max(-1.0, min(1.0, -P.dot(a, b)))))


class TestHelixSetting(unittest.TestCase):
    def test_dihedral_formula_against_normals(self):
        pts = SH.orbit(30, (1, 11), (OUT["chains"][3]["eq_5_1"] and
                                     self._share(30, 1), 1 - self._share(30, 1)))
        c = sum(x * y for x, y in zip(pts[0], pts[1]))
        self.assertAlmostEqual(H.dihedral(c, 4)["deg"], dihedral_from_normals(pts[:4]), places=6)
        self.assertAlmostEqual(H.dihedral(c, 4)["deg"], 72.0, places=6)

    def _share(self, n, p):
        import emlib
        r4 = emlib.load_run(emlib.ROOT / "runs/corkscrew_seed0.json")["outputs"][
            "tetrahedra_chain"]["r4"]["closed_chains"]
        return next(c["r1_squared_share"] for c in r4
                    if c["n_tetrahedra_to_close"] == n and c["p"] == p)

    def test_flat_limit(self):
        # c -> 1: cos d -> 1 / (n - 1), and eq 5.1 becomes tan n th = n tan th
        self.assertAlmostEqual(H.dihedral(1.0, 4)["cos"], 1.0 / 3.0)
        th = math.degrees(math.acos(-2.0 / 3.0)) / 2
        self.assertLess(abs(H.residual(th, 4, 1.0 / 3.0)), 1e-12)
        self.assertGreater(abs(H.residual(th + 1.0, 4, 1.0 / 3.0)), 1e-3)

    def test_tags(self):
        tags = [r["tag"] for r in OUT["chains"]]
        self.assertEqual(tags.count("SETTING_IN"), 8)
        self.assertEqual(tags.count("SETTING_OUTSIDE"), 8)
        self.assertEqual(tags.count(None), 1)
        for r in OUT["chains"]:
            if r["tag"] == "SETTING_OUTSIDE":
                self.assertIsNone(r["eq_5_1"])
            if r["survivor_unmatched"]:
                self.assertEqual(r["tag"], "SETTING_OUTSIDE")

    def test_eq_5_1_in_setting(self):
        w = {}
        for r in OUT["chains"]:
            if r["tag"] == "SETTING_IN":
                self.assertTrue(r["eq_5_1"]["holds_on_every_plane"])
                w[(r["n_to_close"], tuple(r["turns"]))] = r["eq_5_1"]["w"]
        self.assertAlmostEqual(w[(30, (1, 11))], 5.0, places=6)
        self.assertAlmostEqual(w[(30, (7, 13))], 2.5, places=6)
        self.assertAlmostEqual(w[(5, (1, 2))], 3.0, places=6)
        self.assertAlmostEqual(w[(8, (1, 3))], 4.0, places=6)

    def test_slack(self):
        s = OUT["slack"]
        self.assertEqual([(c["n"], c["turns"]) for c in s["5"]["listed"]],
                         [(6, [1, 2, 3]), (10, [1, 3, 5]), (24, [1, 7, 12])])
        self.assertEqual([(c["n"], c["turns"]) for c in s["6"]["listed"]],
                         [(7, [1, 2, 3]), (12, [1, 3, 5]), (18, [1, 5, 8]),
                          (28, [1, 8, 13]), (28, [5, 9, 12])])
        self.assertTrue(OUT["persistence"]["five_and_six_d_together_equal_the_logged_6d_list"])
        for c in s["5"]["listed"]:
            self.assertTrue(c["eq_5_1"]["holds_on_every_plane"])
        for c in s["7"]["listed"][:40]:          # every listed chain is a chain
            pts = SH.orbit(c["n"], c["turns"], c["radius_squared_shares"])
            d = [math.dist(pts[0], pts[j]) for j in range(1, 5)]
            self.assertLess(max(d) - min(d), 1e-8)
            self.assertEqual(SH.rank(pts), 7)

    def test_projection_and_gram(self):
        v = {(p["n"], tuple(p["turns"])): p["verdict"] for p in OUT["projection"]}
        self.assertEqual(v[(24, (1, 7, 12))], ["D13"])
        self.assertEqual(sum(1 for x in v.values() if x == "none"), 3)
        self.assertTrue(OUT["prediction"]["holds"])
        for g in OUT["groups"]:
            self.assertTrue(g["n_exponents_equals_rank"])
        for s in OUT["vertex_sets"]:
            self.assertEqual(s["n_distinct"], s["n_vertices"])
            self.assertEqual(s["antipodal_pairs"], 0)
            self.assertTrue(s["row_0_holds_every_value"])

    def test_solve(self):
        self.assertIsNone(H.solve([[1.0, 2.0], [2.0, 4.0]], [1.0, 2.0]))
        x = H.solve([[2.0, 1.0], [1.0, 3.0]], [3.0, 5.0])
        self.assertAlmostEqual(x[0], 0.8)
        self.assertAlmostEqual(x[1], 1.4)


if __name__ == "__main__":
    unittest.main()
