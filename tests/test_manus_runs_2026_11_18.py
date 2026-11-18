"""Tests for adversarial_probes (oracles and harness), q5_slopes (helpers), and the ported q12 / q11 / q5 / q9
checks. CC0. Co-authored.
Nothing here asserts that a probe is defective: a fixed probe must not turn the suite red. The tests pin the
independent oracles the harness judges by, the mechanism of the greedy-matching case, and the deterministic
readings the run files record."""
import math
import tempfile
import unittest
import _path  # noqa: F401
import adversarial_probes as A
import q12_q11_checks as K
import q5_q9_sweeps as W
import q5_slopes as S
from toolkit import primitives as P
from toolkit import probes as Q

PHI = P.PHI


class Oracles(unittest.TestCase):
    def test_kaspar_schuster_known_values(self):
        for s, v in {"0001101001000101": 6, "": 0, "a": 1, "aaaaaaaa": 2, "ab": 2}.items():
            self.assertEqual(A.ks_complexity(s), v, s)

    def test_border_period_known_values(self):
        for s, v in {"abcabcab": 3, "aaaa": 1, "abcd": 4, "a": 1, "": 0}.items():
            self.assertEqual(A.border_period(s), v, s)

    def test_exact_matching_finds_what_first_fit_misses(self):
        b = [[0.0, 0.0], [1.0, 0.0]]
        a = [[0.5, 0.0], [0.0, 0.0]]
        self.assertFalse(Q._matches(a, b, 0.5))      # first fit sends a[0] to b[0] and strands a[1]
        self.assertTrue(A.has_matching(a, b, 0.5))   # a[0] -> b[1], a[1] -> b[0]

    def test_finite_makes_strict_json(self):
        out = A.finite({"a": [float("nan"), float("inf"), -float("inf"), 1.5], "b": 2})
        self.assertEqual(out, {"a": ["NaN", "Infinity", "-Infinity", 1.5], "b": 2})

    def test_small_run_is_complete_and_classified(self):
        out = A.run(0, small=True)
        ids = [c["id"] for c in out["cases"]]
        self.assertEqual(ids, [f"A{i:02d}" for i in range(1, 21)])
        for c in out["cases"]:
            self.assertIn(c["class"], ("HOLDS", "DEFECT", "LIMIT", "UNDEFINED"))
            self.assertIn("observed", c)
        by = {c["id"]: c["class"] for c in out["cases"]}
        self.assertEqual(by["A14"], "HOLDS")   # lz76 against the independent algorithm
        self.assertEqual(by["A17"], "HOLDS")   # smallest_period against the border computation
        self.assertEqual(by["A20"], "HOLDS")   # the writer refuses bad records


class Q5Slopes(unittest.TestCase):
    def test_golden_canonical_width_gives_two_and_a_generic_width_three(self):
        s = 1 / PHI
        self.assertEqual(S.gaps_at(s, 1 + s)["n_gap_lengths"], 2)
        r = S.gaps_at(s, 0.5 * (1 + s))
        self.assertEqual(r["n_gap_lengths"], 3)
        self.assertEqual(r["three_gap_sum"]["verdict"], "pass")

    def test_records_of_the_golden_slope(self):
        rec = S.records(1 / PHI, 6)
        self.assertEqual(rec[(2, 1)][0], "record_above")
        self.assertEqual(rec[(3, 2)][0], "record_below")
        self.assertAlmostEqual(rec[(2, 1)][1], 2 / PHI - 1, places=12)
        self.assertEqual(rec[(0, 1)][0], "start")

    def test_continued_fraction_errors_shrink_by_phi(self):
        errs = S.cf_errors(1 / PHI, 13)
        self.assertEqual(errs[0], (0, 1, 1.0))
        d = [e[2] for e in errs]
        for a, b in zip(d[2:], d[3:]):
            self.assertAlmostEqual(a / b, PHI, places=6)

    def test_width_is_symmetric_in_the_sign_of_the_slope(self):
        a = S.gaps_at(1 / PHI, 1 + 1 / PHI)["n_gap_lengths"]
        b = S.gaps_at(-1 / PHI, 1 + 1 / PHI)["n_gap_lengths"]
        self.assertEqual((a, b), (2, 2))


class PortedChecks(unittest.TestCase):
    def test_q12_extended_matches_its_prediction(self):
        with tempfile.TemporaryDirectory() as d:
            K.OUT = d
            _, o = K.q12()
        self.assertEqual(o["prediction_status"], "CONFORMS")
        got = {(r["k"], r["k_max"]): r["got"] for r in o["kgons"]}
        self.assertEqual(got[(13, 12)], [])
        self.assertEqual(got[(13, 24)], [13])
        self.assertEqual(got[(24, 24)], [2, 3, 4, 6, 8, 12, 24])

    def test_q11_aligned_readings(self):
        with tempfile.TemporaryDirectory() as d:
            K.OUT = d
            _, o = K.q11()
        self.assertTrue(all(r["matches_within_1e-9"] for r in o["phi_composition"]))
        self.assertTrue(all(r["pairs_with_identical_vectors"] == 400 for r in o["orthogonality_aligned_stream"]))
        self.assertTrue(all(r["identical_after_unit_circumradius_and_sort"] is True
                            for r in o["polytopes"] if r["shape"] != "simplex"))
        self.assertTrue(o["constants"]["within_2_ulp_both"])

    def test_q9_fit_recovers_a_known_crossing(self):
        ds = list(range(30, 101, 5))
        ps = [0.5 + 0.004 * (d - 61.0) for d in ds]
        cross, slope = W.fit_cross(ds, ps)
        self.assertAlmostEqual(cross, 61.0, places=9)
        self.assertAlmostEqual(slope, 0.004, places=12)

    def test_q5_canonical_strip_gives_two_lengths(self):
        _, perp = W.C.line_basis_2d(W.SLOPE)
        lo, hi = W.C.strip_bounds_1d(perp[0])
        n, cl = W.gap_lengths((hi - lo) / 2.0)
        self.assertEqual(len(cl), 2)


if __name__ == "__main__":
    unittest.main()
