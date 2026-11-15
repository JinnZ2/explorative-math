import math
import unittest

import _path  # noqa: F401
import golden_angle as E


class GoldenAngle(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = E.run()

    def test_golden_value(self):
        self.assertAlmostEqual(E.GOLDEN, 137.50776405, places=6)

    def test_point_radii(self):
        pts = E.points(E.GOLDEN, n=10)
        for k, (x, y) in enumerate(pts, 1):
            self.assertAlmostEqual(math.hypot(x, y), math.sqrt(k), places=12)

    def test_nn_distance_known_case(self):
        self.assertEqual(E.nn_distances([(0, 0), (3, 4), (0, 1)], skip=0), [1.0, 4.242640687119285, 1.0])

    def test_rational_angles_make_spokes(self):
        g = self.rows["golden"]["min_nn"]
        for k in ("rational_2/5_turn", "rational_3/8_turn"):
            self.assertLess(self.rows[k]["min_nn"], g)

    def test_every_angle_recorded(self):
        self.assertEqual(set(self.rows), set(E.ANGLES))


if __name__ == "__main__":
    unittest.main()
