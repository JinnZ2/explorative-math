import math
import unittest

import _path  # noqa: F401
import shell_concentration as E


class Shell(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = E.run(seed=0)

    def test_exact_formula(self):
        self.assertAlmostEqual(E.exact(1, 0.1), 0.1)
        self.assertAlmostEqual(E.exact(2, 0.5), 0.75)

    def test_ball_to_cube(self):
        self.assertAlmostEqual(math.exp(E.log_ball_to_cube(2)), math.pi / 4, places=12)
        self.assertAlmostEqual(math.exp(E.log_ball_to_cube(3)), math.pi / 6, places=12)
        E.log_ball_to_cube(500)  # no overflow

    def test_fraction_rises_to_one(self):
        for e in E.EPS:
            col = [r["exact"] for r in self.rows if r["eps"] == e]
            self.assertEqual(col, sorted(col))
        self.assertGreater(next(r["exact"] for r in self.rows if r["d"] == 500 and r["eps"] == 0.01), 0.99)

    def test_mc_void_above_cutoff_and_within_5_sigma_below(self):
        for r in self.rows:
            if r["d"] > E.MC_MAX_DIM:
                self.assertIsNone(r["monte_carlo"])
                self.assertIsNone(r["z"])
            else:
                self.assertLess(abs(r["z"]), 5.0)


if __name__ == "__main__":
    unittest.main()
