import unittest

import _path  # noqa: F401
import fibonacci_chain as E

FIB = {1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610}


class FibonacciChain(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = E.run()
        cls.g = cls.out["golden"]

    def test_two_lengths_ratio_phi(self):
        self.assertEqual(len(self.g["distinct_gap_lengths"]), 2)
        self.assertAlmostEqual(self.g["L_over_S_length"], E.PHI, places=9)

    def test_count_ratio_near_phi(self):
        self.assertAlmostEqual(self.g["count_ratio"], E.PHI, delta=0.02)

    def test_is_factor_of_substitution_word(self):
        self.assertTrue(self.g["factor_of_substitution_word"])

    def test_aperiodic_signature_vs_control(self):
        per = [self.g["prefix_smallest_periods"][str(k)] for k in E.PREFIX_LENGTHS
               if str(k) in self.g["prefix_smallest_periods"]]
        self.assertEqual(per, sorted(per))
        self.assertGreater(per[-1], per[0])
        ctrl = set(self.out["control_rational_1_2"]["prefix_smallest_periods"].values())
        self.assertEqual(len(ctrl), 1)

    def test_observed_prefix_periods_are_fibonacci(self):
        # OBSERVED on this run, pinned so a change shows up
        self.assertTrue(set(self.g["prefix_smallest_periods"].values()) <= FIB)

    def test_smallest_period(self):
        self.assertEqual(E.smallest_period("abab"), 2)
        self.assertEqual(E.smallest_period("abc"), 3)


if __name__ == "__main__":
    unittest.main()
