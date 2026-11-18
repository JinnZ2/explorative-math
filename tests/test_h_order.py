"""Tests for experiments/h_order.py. CC0. Co-authored."""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
import emlib  # noqa: E402
import h_order as X  # noqa: E402

OUT = X.run(0)
S = OUT["summary"]


class TestOrder(unittest.TestCase):
    def test_order_kept(self):
        self.assertTrue(X.order_kept(28, [1, 8, 13], 1))
        self.assertEqual(X.images(28, [1, 8, 13], 5), [5, 12, 9])
        self.assertFalse(X.order_kept(28, [1, 8, 13], 5))
        self.assertTrue(X.order_kept(28, [1, 8, 13], 3))          # (3, 4, 11)
        self.assertFalse(X.order_kept(28, [1, 8, 13], 13))        # (13, 8, 1): the set, reversed
        self.assertEqual(X.order_kept(30, [1, 11], 7), X.order_kept(30, [1, 11], 23))  # u, N - u

    def test_units_reaching(self):
        self.assertEqual(X.units_reaching(28, [1, 8, 13], [5, 9, 12]), [5, 9, 19, 23])
        self.assertEqual(X.units_reaching(28, [1, 8, 13], [1, 8, 13]), [1, 13, 15, 27])

    def test_variant_is_the_same_test(self):
        self.assertTrue(S["variant_is_the_same_test"])
        for n, p in ((28, [1, 8, 13]), (18, [1, 5, 8]), (20, [1, 4, 7, 10])):
            for u in range(1, n):
                if X.math.gcd(u, n) == 1:
                    self.assertEqual(X.order_kept(n, p, u), X.cosine_order_kept(n, p, u))


class TestCrossTab(unittest.TestCase):
    def test_rows(self):
        self.assertEqual((S["rows"], S["rows_real"], S["rows_negative"]), (5183, 1444, 3739))
        for how in ("first", "any", "every"):
            t = OUT["cross_tab"][how]["all_rows"]
            self.assertEqual(t["kept_real"] + t["kept_negative"] + t["broken_real"]
                             + t["broken_negative"], 5183)

    def test_h_order_does_not_hold(self):
        self.assertEqual(S["h_order_holds"], {"first": False, "any": False, "every": False})
        self.assertEqual(S["counterexamples"], {"first": 1079, "any": 1092, "every": 1126})
        t = OUT["cross_tab"]["any"]["all_rows"]
        self.assertEqual((t["kept_real"], t["kept_negative"], t["broken_real"], t["broken_negative"]),
                         (844, 492, 600, 3247))

    def test_named_counterexamples_in_both_cells(self):
        # order broken by every unit, and a found chain
        self.assertEqual(S["named_case_28_1_8_13_to_5_9_12"],
                         {"units_reaching": [5, 9, 19, 23],
                          "order_kept_by_unit": [False, False, False, False], "real": True})
        # order kept, and a negative share
        kept_neg = OUT["cross_tab"]["any"]["order_kept_and_a_negative_share_first_10"]
        self.assertIn([28, [1, 8, 13], [3, 4, 11]], [[r["n"], r["chain"], r["member"]] for r in kept_neg])
        self.assertEqual(len(kept_neg), 10)

    def test_a_table_can_come_out_clean(self):
        rows = [{"any": True, "real": True}, {"any": False, "real": False}]
        self.assertTrue(X.table(rows, "any")["h_order_holds_on_every_row"])
        rows.append({"any": True, "real": False})
        self.assertEqual(X.table(rows, "any")["counterexamples"], 1)

    def test_logged_run_matches_and_slots_are_void(self):
        rec = emlib.load_run(emlib.ROOT / "runs/h_order_seed0.json.gz")
        self.assertEqual(rec["outputs"]["summary"], S)
        self.assertIsNone(rec["human_interest"])
        self.assertIsNone(rec["ai_interest"])


if __name__ == "__main__":
    unittest.main()
