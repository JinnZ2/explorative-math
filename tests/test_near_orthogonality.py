import math
import random
import unittest

import _path  # noqa: F401
import emlib
import near_orthogonality as E


class NearOrthogonality(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.per_dim, cls.angles = E.run(seed=0)

    def test_unit_vectors(self):
        rng = random.Random(1)
        for d in (2, 50):
            self.assertAlmostEqual(emlib.norm(emlib.random_unit(rng, d)), 1.0, places=12)

    def test_all_angles_recorded(self):
        for d in E.DIMS:
            self.assertEqual(len(self.angles[str(d)]), E.PAIRS)

    def test_spread_shrinks_with_dimension(self):
        sds = [self.per_dim[str(d)]["sd_angle"] for d in E.DIMS]
        self.assertEqual(sds, sorted(sds, reverse=True))

    def test_sd_cos_tracks_one_over_sqrt_d(self):
        for d in E.DIMS:
            self.assertAlmostEqual(self.per_dim[str(d)]["sd_cos_times_sqrt_d"], 1.0, delta=0.15)

    def test_deterministic(self):
        self.assertEqual(E.run(seed=0)[1], self.angles)


if __name__ == "__main__":
    unittest.main()
