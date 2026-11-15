import itertools
import math
import random
import unittest

import _path  # noqa: F401
import simplex_4d as E


class Simplex4D(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = E.run(seed=0, trials=200)

    def test_regular(self):
        v4, _ = E.simplex_vertices()
        self.assertEqual(len(v4), 5)
        self.assertTrue(all(len(v) == 4 for v in v4))
        e = E.edge_lengths(v4)
        self.assertAlmostEqual(min(e), max(e), places=12)

    def test_rotation_preserves_distances(self):
        v4, _ = E.simplex_vertices()
        r = E.rotate(E.rotate(v4, 0, 3, 0.7), 1, 2, 2.1)
        for a, b in itertools.combinations(range(5), 2):
            self.assertAlmostEqual(math.dist(v4[a], v4[b]), math.dist(r[a], r[b]), places=12)

    def test_petrie_is_regular_pentagon(self):
        self.assertLess(self.out["petrie_score"], 1e-9)

    def test_score_detects_non_pentagon(self):
        sq = [(1, 0), (0, 1), (-1, 0), (0, -1), (0.5, 0.5)]
        self.assertGreater(E.pentagon_score(sq)[0], 1.0)

    def test_all_scores_recorded(self):
        self.assertEqual(len(self.out["scores"]), 200)
        self.assertEqual(self.out["n_logged_below"], len(self.out["logged"]))


if __name__ == "__main__":
    unittest.main()
