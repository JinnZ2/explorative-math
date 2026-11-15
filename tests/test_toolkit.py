import json
import math
import tempfile
import unittest

import _path  # noqa: F401
import emlib
from toolkit import composers as C
from toolkit import generator as G
from toolkit import loghook
from toolkit import primitives as P
from toolkit import probes as Q


def regular(points):
    e = P.edge_lengths(points)
    m = min(e)
    return sum(1 for x in e if abs(x - m) < 1e-9)


class Primitives(unittest.TestCase):
    def test_rotation_in_any_plane_preserves_norms(self):
        pts = P.hypercube(5)
        r = P.rotate_many(pts, [((0, 4), 0.3), ((2, 3), 1.9)])
        for a, b in zip(pts, r):
            self.assertAlmostEqual(P.norm(a), P.norm(b), places=12)

    def test_polytope_vertex_counts(self):
        for n in (3, 4, 5):
            self.assertEqual(len(P.simplex(n)), n + 1)
            self.assertEqual(len(P.hypercube(n)), 2 ** n)
            self.assertEqual(len(P.cross_polytope(n)), 2 * n)
        self.assertEqual(len(P.cell24()), 24)
        self.assertEqual(len(P.icosahedron()), 12)

    def test_edge_counts_at_minimum_distance(self):
        self.assertEqual(regular(P.simplex(4)), 10)        # complete graph K5
        self.assertEqual(regular(P.hypercube(4)), 32)
        self.assertEqual(regular(P.cross_polytope(4)), 24)
        self.assertEqual(regular(P.cell24()), 96)
        self.assertEqual(regular(P.icosahedron()), 30)

    def test_lattices(self):
        self.assertEqual(len(P.lattice_zn(3, 1)), 27)
        h = P.lattice_hex(1)
        self.assertEqual(len(h), 9)
        self.assertAlmostEqual(P.min_edge(h), 1.0)

    def test_projection_onto_basis(self):
        self.assertEqual(P.project([[1, 2, 3]], basis=[[0, 0, 1]]), [[3]])


class Composers(unittest.TestCase):
    def test_cut_and_project_golden_slope(self):
        par, perp = C.line_basis_2d(1 / P.PHI)
        pts = C.cut_and_project(P.lattice_zn(2, 40), par, perp, C.strip_window_1d(perp[0]))
        xs = sorted(p[0] for p in pts if abs(p[0]) < 15)
        gaps = sorted({round(b - a, 9) for a, b in zip(xs, xs[1:])})
        self.assertEqual(len(gaps), 2)
        self.assertAlmostEqual(gaps[1] / gaps[0], P.PHI, places=6)

    def test_phi_scale_and_stack(self):
        self.assertAlmostEqual(C.phi_scale([[1.0, 0.0]], 2)[0][0], P.PHI ** 2)
        st = C.phi_stack([[1.0, 0.0]], (0, 1))
        self.assertEqual([r["tag"] for r in st], ["phi^0", "phi^1"])

    def test_iterate_and_overlay(self):
        self.assertEqual(C.iterate(lambda x: x * 2, 1, 3), [1, 2, 4, 8])
        o = C.overlay([[0, 0]], [[1, 1]])
        self.assertEqual([r["tag"] for r in o], ["A", "B"])


class Probes(unittest.TestCase):
    def test_rotational_symmetry(self):
        self.assertEqual(Q.rotational_symmetry(P.regular_polygon(5)), [5])
        self.assertEqual(Q.rotational_symmetry(P.regular_polygon(6)), [2, 3, 6])
        self.assertEqual(Q.rotational_symmetry([[0, 0], [1, 0], [0, 2]]), [])

    def test_mirror_symmetry(self):
        self.assertEqual(len(Q.mirror_symmetry(P.regular_polygon(4))), 4)
        self.assertEqual(Q.mirror_symmetry([[0, 0], [1, 0], [0, 2], [3, 7]]), [])

    def test_periodicity(self):
        self.assertEqual(Q.smallest_period("abcabcab"), 3)
        self.assertEqual(Q.period_profile("ab" * 20)["32"], 2)

    def test_zlib_ratio_orders_repetitive_below_random(self):
        import random
        rng = random.Random(0)
        rep = "ab" * 500
        rnd = "".join(rng.choice("ab") for _ in range(1000))
        self.assertLess(Q.zlib_ratio(rep), Q.zlib_ratio(rnd))
        self.assertIsNone(Q.zlib_ratio(""))

    def test_angle_histogram_counts_everything_in_range(self):
        h = Q.angle_histogram([0, 5, 95, 180], width=90)
        self.assertEqual(h["counts"], [2, 2])


class LogHookAndGenerator(unittest.TestCase):
    def test_log_run_voids_interest_slots(self):
        with tempfile.TemporaryDirectory() as d:
            p = loghook.log_run("t", 3, {"a": 1}, {"x": 2}, symbols=list("abab"),
                                definition="test", out_dir=d)
            rec = json.loads(p.read_text())
        self.assertEqual(tuple(rec), emlib.SCHEMA)
        self.assertIsNone(rec["human_interest"])
        self.assertIsNone(rec["ai_interest"])
        self.assertIsNone(rec["divergence_note"])
        self.assertIn("zlib_ratio", rec["compressibility"])

    def test_logged_decorator(self):
        @loghook.logged("dec")
        def f(seed):
            return {"s": seed}, {"y": 1}, None, None, None
        with tempfile.TemporaryDirectory() as d:
            p = f(7, out_dir_override=d)
            self.assertTrue(p.name.endswith("dec_seed7.json"))

    def test_generator_deterministic_and_covers_all_pipelines(self):
        self.assertEqual(G.draw_spec(5), G.draw_spec(5))
        kinds = {G.draw_spec(s)["pipeline"] for s in range(30)}
        self.assertEqual(kinds, set(G.PIPELINES))
        for s in range(30):
            spec, outputs, angles, symbols, definition = G.run(s)
            self.assertEqual(spec, G.draw_spec(s))
            self.assertTrue(definition)

    def test_generator_writes_schema(self):
        with tempfile.TemporaryDirectory() as d:
            rec = json.loads(G.generate(1, out_dir=d).read_text())
        emlib.check_run(rec)
        self.assertIsNone(rec["ai_interest"])


if __name__ == "__main__":
    unittest.main()
