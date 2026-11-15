"""Tests for exact/ (ZPhi, Cyc5, inflation) and exact_inflation. CC0. Co-authored."""
import ast
import unittest
from fractions import Fraction
from pathlib import Path

import _path  # noqa: F401
import exact_inflation as E
from exact import Cyc5, ONE, PHI, ZERO, ZETA, ZPhi, cross, dot, norm2
from exact import inflation as I


class ZPhiOps(unittest.TestCase):
    def test_defining_relation_and_inverse(self):
        self.assertEqual(PHI * PHI, PHI + 1)
        self.assertEqual(PHI.inverse(), PHI - 1)
        self.assertEqual(PHI * PHI.inverse(), ONE)
        self.assertEqual(ZPhi(3, -2) / ZPhi(3, -2), ONE)
        with self.assertRaises(ZeroDivisionError):
            ZERO.inverse()

    def test_add_mul_against_fibonacci(self):
        # phi^n = F(n-1) + F(n) phi
        f0, f1, x = 0, 1, PHI
        for _ in range(60):
            self.assertEqual(x, ZPhi(f0, f1))
            x, f0, f1 = x * PHI, f1, f0 + f1
        self.assertEqual(PHI ** 10, ZPhi(34, 55))
        self.assertEqual(PHI ** -3 * PHI ** 3, ONE)

    def test_conjugate_and_norm(self):
        x = ZPhi(Fraction(3, 7), Fraction(-5, 2))
        self.assertEqual(x.conjugate().conjugate(), x)
        self.assertEqual(x * x.conjugate(), ZPhi(x.norm(), 0))
        self.assertEqual(PHI.conjugate(), 1 - PHI)
        self.assertEqual(PHI.norm(), -1)

    def test_compare_is_exact(self):
        self.assertGreater(PHI, ONE)
        self.assertLess(PHI, ZPhi(Fraction(1619, 1000), 0))
        self.assertGreater(PHI, ZPhi(Fraction(1618, 1000), 0))
        self.assertEqual(ZPhi(3, -2).sign(), -1)          # 3 - 2 phi < 0
        self.assertEqual(ZPhi(-3, 2).sign(), 1)
        self.assertEqual(ZERO.sign(), 0)
        # consecutive Fibonacci quotients straddle phi, however close they get
        f0, f1 = 1, 1
        for k in range(200):
            q = ZPhi(Fraction(f1, f0), 0)
            self.assertEqual((q - PHI).sign(), 1 if k % 2 else -1)
            f0, f1 = f1, f0 + f1
        self.assertEqual(sorted([PHI, ZERO, ONE, PHI - 1, -PHI]),
                         [-PHI, ZERO, PHI - 1, ONE, PHI])

    def test_bounds_are_rational_and_straddle(self):
        lo, hi = PHI.bounds()
        self.assertIsInstance(lo, Fraction)
        self.assertTrue(ZPhi(lo, 0) < PHI < ZPhi(hi, 0))

    def test_no_floats(self):
        for bad in (0.5, 1.0):
            with self.assertRaises(TypeError):
                ZPhi(bad, 0)
            with self.assertRaises(TypeError):
                PHI + bad
            with self.assertRaises(TypeError):
                PHI * bad
        with self.assertRaises(TypeError):
            float(PHI)
        with self.assertRaises(TypeError):
            float(ZETA)
        with self.assertRaises(TypeError):
            ZPhi(True, 0)

    def test_no_float_in_the_source_of_exact(self):
        # AST walk: no float literal, no call to float(), no import of math
        root = Path(__file__).resolve().parents[1] / "exact"
        for f in sorted(root.glob("*.py")):
            tree = ast.parse(f.read_text())
            for node in ast.walk(tree):
                if isinstance(node, ast.Constant):
                    self.assertNotIsInstance(node.value, float, f"{f.name}: float literal")
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                    self.assertNotEqual(node.func.id, "float", f"{f.name}: float() call")
                if isinstance(node, ast.Import):
                    self.assertNotIn("math", [a.name for a in node.names], f.name)
                if isinstance(node, ast.ImportFrom):
                    self.assertNotEqual(node.module, "math", f.name)

    def test_planted_float_is_caught_by_the_same_walk(self):
        tree = ast.parse("x = 0.5\ny = float(3)\n")
        hits = [n for n in ast.walk(tree)
                if (isinstance(n, ast.Constant) and isinstance(n.value, float))
                or (isinstance(n, ast.Call) and getattr(n.func, "id", "") == "float")]
        self.assertEqual(len(hits), 2)

    def test_immutable_and_hashable(self):
        with self.assertRaises(AttributeError):
            PHI.a = 3
        self.assertEqual(len({PHI, ZPhi(0, 1), ONE}), 2)


class Cyc5Ops(unittest.TestCase):
    def test_fifth_root_of_unity(self):
        z = Cyc5(1, 0)
        for _ in range(5):
            z = z * ZETA
        self.assertEqual(z, Cyc5(1, 0))
        self.assertEqual(norm2(ZETA), ONE)
        self.assertEqual(sum(I.pentagon()[1:], I.pentagon()[0]), Cyc5(0, 0))

    def test_pentagon_diagonal_over_side_is_phi(self):
        p = I.pentagon()
        self.assertEqual(norm2(p[2] - p[0]), norm2(p[1] - p[0]) * PHI * PHI)

    def test_cross_and_dot(self):
        p = I.pentagon()
        self.assertEqual(cross(p[1], p[1] * PHI), ZERO)
        self.assertNotEqual(cross(p[0], p[1]), ZERO)
        self.assertEqual(dot(p[0], p[0]), ONE)
        self.assertEqual(cross(p[0], p[1]), -cross(p[1], p[0]))


class Inflation(unittest.TestCase):
    def test_fibonacci_two_constructions_equal(self):
        word, pos = "L", I.fib_positions("L")
        for n in range(1, 13):
            gw, gp = I.fib_inflate(word, pos)
            sw = I.fib_substitute(word)
            self.assertEqual(gw, sw)
            self.assertEqual(gp, I.fib_positions(sw))
            self.assertEqual(gp[-1], PHI ** (n + 1))
            word, pos = gw, gp

    def test_penrose_area_and_edges(self):
        tris = I.penrose_seed()
        a0 = abs(I.area_over_w(tris[0]))
        for _ in range(5):
            tris = I.penrose_subdivide(tris)
        self.assertEqual(sum((abs(I.area_over_w(t)) for t in tris), ZERO), a0)
        for t in tris:
            for x in I.edge_norms(t):
                self.assertIsNotNone(I.phi_power(x))


class ExactInflationRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = E.run(levels={"fibonacci": 10, "pentagon": 10, "penrose": 5})

    def test_zero_drift_at_ten_levels(self):
        for key in ("fibonacci", "pentagon", "penrose"):
            self.assertTrue(self.out[key]["zero_drift"], key)
        self.assertEqual(sum(r["n_position_mismatches"] for r in self.out["fibonacci"]["rows"]), 0)

    def test_logged_run_has_ten_or_more_levels_and_zero_drift(self):
        rec = E.emlib.load_run(E.emlib.RUNS / "exact_inflation_seed0.json")
        o = rec["outputs"]
        for key in ("fibonacci", "pentagon", "penrose"):
            self.assertGreaterEqual(o[key]["levels"], 10)
            self.assertTrue(o[key]["zero_drift"])
        self.assertIsNone(rec["human_interest"])
        # OBSERVED, pinned: floats do drift on the same construction
        self.assertGreater(o["fibonacci"]["rows"][-1]["float_contrast_max_abs_difference"], 0.0)
        small = E.emlib.load_run(E.emlib.RUNS / "exact_inflation_small_seed0.json")
        self.assertEqual(set(small["outputs"]), set(o))


if __name__ == "__main__":
    unittest.main()
