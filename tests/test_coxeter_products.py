"""Tests for experiments/coxeter_products.py. CC0. Co-authored."""
import unittest

import _path  # noqa: F401
import coxeter_products as X
import emlib


class Factors(unittest.TestCase):
    def test_computed_factors(self):
        self.assertEqual((X.factor("A", 1)["h"], X.factor("A", 1)["exponents"]), (2, [1]))
        self.assertEqual((X.factor("G", 2)["h"], X.factor("G", 2)["exponents"]), (6, [1, 5]))
        self.assertEqual((X.factor("F", 4)["h"], X.factor("F", 4)["exponents"]),
                         (12, [1, 5, 7, 11]))
        i = X.factor("I", 2, 7)
        self.assertEqual((i["name"], i["h"], i["exponents"]), ("I2(7)", 7, [1, 6]))

    def test_products_respect_the_rank_limit(self):
        fs = [X.factor("A", 1), X.factor("A", 4), X.factor("F", 4)]
        for parts in X.products(fs):
            self.assertLessEqual(sum(p["rank"] for p in parts), X.MAX_RANK)
            self.assertIn(len(parts), (2, 3))

    def test_logged_run(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = emlib.load_run(emlib.RUNS / "coxeter_products_seed0.json")
        self.assertIsNone(rec["ai_interest"])
        self.assertEqual(rec["outputs"]["summary"], {"chains": 4, "match": 0, "none": 4})
        for r in rec["outputs"]["chains"]:
            self.assertEqual(r["verdict"], "none")
            self.assertTrue(r["turn_fractions_no_factor_offers"])


if __name__ == "__main__":
    unittest.main()
