"""Tests for exact/cyclo.py. CC0. Co-authored."""
import math
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from exact.cyclo import Field, cos_bounds, cyclotomic, pi_bounds, show, solve_exact  # noqa: E402


class TestCyclotomic(unittest.TestCase):
    def test_known_polynomials(self):
        self.assertEqual(cyclotomic(1), [-1, 1])
        self.assertEqual(cyclotomic(4), [1, 0, 1])
        self.assertEqual(cyclotomic(6), [1, -1, 1])
        self.assertEqual(cyclotomic(7), [1] * 7)
        self.assertEqual(cyclotomic(18), [1, 0, 0, -1, 0, 0, 1])

    def test_degree_is_euler_phi(self):
        for n in range(1, 41):
            phi = sum(1 for u in range(1, n + 1) if math.gcd(u, n) == 1)
            self.assertEqual(len(cyclotomic(n)) - 1, phi, n)


class TestField(unittest.TestCase):
    def test_cos_against_floats(self):
        for n in (5, 7, 18, 24, 28, 30):
            k = Field(n)
            for j in range(-3, n + 3):
                self.assertAlmostEqual(float(k.cos(j).approx()), math.cos(2 * math.pi * j / n),
                                       places=12)
                self.assertTrue(k.cos(j).is_real())
                self.assertEqual(k.sin_i(j) * k.sin_i(j), k.cos(j) * k.cos(j) - 1)

    def test_exact_values(self):
        k = Field(12)
        self.assertEqual(k.cos(2), F(1, 2))                 # cos 60
        self.assertEqual(k.cos(3), 0)                       # cos 90
        self.assertEqual(k.cos(1) * k.cos(1), F(3, 4))      # cos^2 30
        self.assertEqual(k.zeta(12), 1)
        self.assertEqual(k.zeta(6), -1)

    def test_inverse_and_division(self):
        k = Field(28)
        x = k.cos(1) * 3 - k.cos(5) + F(2, 7)
        self.assertEqual(x * x.inverse(), 1)
        self.assertEqual((x / x), k.one)
        self.assertEqual(1 / x, x.inverse())
        with self.assertRaises(ZeroDivisionError):
            k.zero.inverse()

    def test_galois(self):
        k = Field(28)
        for u in k.units():
            self.assertEqual(k.cos(3).galois(u), k.cos(3 * u))
            x, y = k.cos(1) + 2, k.cos(4) - k.cos(9)
            self.assertEqual((x * y).galois(u), x.galois(u) * y.galois(u))
        with self.assertRaises(ValueError):
            k.cos(1).galois(2)

    def test_sign(self):
        k = Field(18)
        self.assertEqual((k.cos(1) - k.cos(1)).sign(), 0)
        self.assertEqual(k.cos(1).sign(), 1)
        self.assertEqual(k.cos(8).sign(), -1)
        with self.assertRaises(ValueError):
            k.zeta(1).sign()

    def test_bounds_are_rational_and_bracket(self):
        lo, hi = pi_bounds(20)
        self.assertTrue(F(314159265358979323846, 10 ** 20) < lo < hi
                        < F(314159265358979323847, 10 ** 20))
        for k, n in ((1, 18), (5, 28), (13, 28), (0, 7), (9, 18)):
            lo, hi = cos_bounds(k, n, 16)
            self.assertIsInstance(lo, F)
            self.assertTrue(lo <= hi)
            self.assertTrue(float(lo) - 1e-12 <= math.cos(2 * math.pi * k / n) <= float(hi) + 1e-12)

    def test_sign_of_a_very_small_element(self):
        k = Field(28)
        tiny = (k.cos(1) - k.cos(1).approx(30))             # |tiny| < 1e-30, not zero
        self.assertIn(tiny.sign(), (-1, 1))
        self.assertEqual((tiny - tiny).sign(), 0)
        x = k.cos(3)
        self.assertLess(abs(x.approx(25) - x.approx(40)), F(2, 10 ** 25))

    def test_coordinates(self):
        k = Field(18)
        c = k.cos(2)                                        # cos(2 pi / 9)
        gens = [k.one, c, c * c]
        self.assertEqual(k.coordinates(c * c * c, gens), [F(1, -8), F(3, 4), 0])   # 4c^3 = 3c - 1/2
        self.assertIsNone(k.coordinates(k.zeta(1), gens))           # not real
        self.assertIsNone(k.coordinates(c, [k.one, k.one]))         # dependent generators
        k28 = Field(28)
        c7 = k28.cos(4)
        self.assertIsNone(k28.coordinates(k28.cos(1), [k28.one, c7, c7 * c7]))

    def test_other_field_refused(self):
        with self.assertRaises(ValueError):
            Field(18).of(Field(28).cos(1))

    def test_solve_exact(self):
        k = Field(7)
        a = [[k.cos(1), k.one], [k.one, k.cos(2)]]
        x = [k.cos(3), k.of(F(1, 3))]
        b = [a[0][0] * x[0] + a[0][1] * x[1], a[1][0] * x[0] + a[1][1] * x[1]]
        self.assertEqual(solve_exact(a, b, k.zero, k.one), x)
        self.assertIsNone(solve_exact([[k.one, k.one], [k.one, k.one]], b, k.zero, k.one))

    def test_show(self):
        self.assertEqual(show([F(1, 2), 0, F(-2)], ["1", "c", "c^2"]), "1/2 - 2*c^2")
        self.assertEqual(show([0, 0], ["1", "c"]), "0")


if __name__ == "__main__":
    unittest.main()
