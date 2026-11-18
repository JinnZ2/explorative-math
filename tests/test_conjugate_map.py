"""Tests for experiments/conjugate_map.py. CC0. Co-authored."""
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
import conjugate_map as X  # noqa: E402
import emlib  # noqa: E402
from exact.cyclo import Field  # noqa: E402

OUT = X.run(0)
S = OUT["summary"]


class TestOrbit(unittest.TestCase):
    def test_named_checks(self):
        self.assertEqual(X.reduce_set(30, (1, 11), 7)[0], (7, 13))
        self.assertEqual(X.reduce_set(28, (1, 8, 13), 5)[0], (5, 9, 12))
        self.assertTrue(S["named_checks_hold"])

    def test_orbits_of_the_survivors(self):
        self.assertEqual(sorted(X.orbit_of(28, (1, 8, 13))), [(1, 8, 13), (3, 4, 11), (5, 9, 12)])
        self.assertEqual(sorted(X.orbit_of(18, (1, 5, 8))), [(1, 2, 7), (1, 5, 8), (4, 5, 7)])
        self.assertEqual(sorted(X.orbit_of(30, (1, 11))), [(1, 11), (7, 13)])
        self.assertEqual(sorted(X.orbit_of(24, (1, 7, 12))), [(1, 7, 12), (5, 11, 12)])

    def test_orbit_is_closed_and_the_same_from_any_member(self):
        for r in OUT["orbits"][:40]:
            sets = sorted(tuple(x["turns"]) for x in r["orbit"])
            for s in sets:
                self.assertEqual(sorted(X.orbit_of(r["n"], s)), sets)

    def test_every_chain_is_in_its_own_orbit_and_found(self):
        for r in OUT["orbits"]:
            own = [x for x in r["orbit"] if x["turns"] == r["turns"]]
            self.assertEqual((len(own), own[0]["verdict"]), (1, "FOUND"))

    def test_chains_read(self):
        self.assertEqual((S["chains_read"], S["of_them_from_the_r7_list"]), (782, 765))


class TestVerdicts(unittest.TestCase):
    def test_the_three_missing_sets(self):
        self.assertEqual(S["missing_three_verdicts"],
                         [[28, [3, 4, 11], "NEGATIVE_SHARE"], [18, [4, 5, 7], "NEGATIVE_SHARE"],
                          [18, [1, 2, 7], "NEGATIVE_SHARE"]])
        self.assertTrue(S["prediction_holds_on_the_three"])
        neg = {(m["n"], tuple(m["turns"])): m["negative_on_turns"] for m in OUT["missing"]}
        self.assertEqual(neg, {(28, (3, 4, 11)): [3], (18, (4, 5, 7)): [5],
                               (18, (1, 2, 7)): [2, 7]})

    def test_exact_and_float_solvers_agree_on_them(self):
        for m in OUT["missing"]:
            self.assertEqual(m["verdict"], m["float_solver_verdict"])
            f = X.float_shares(m["n"], m["turns"], 4)["shares"]
            for a, b in zip(m["shares"], f):
                self.assertAlmostEqual(a, b, places=9)
            self.assertAlmostEqual(sum(m["shares"]), 1.0, places=12)

    def test_no_member_is_without_a_solution_or_real_and_unlogged(self):
        self.assertFalse(S["any_no_solution_anywhere"])
        self.assertEqual(S["any_real_member_that_is_not_a_logged_chain"], [])
        self.assertEqual(S["not_found_members_by_verdict"],
                         {"R^7 list": {"NEGATIVE_SHARE": 2474},
                          "logged and R^5 / R^6": {"NEGATIVE_SHARE": 4}})

    def test_a_found_chain_solves_real(self):
        for n, t, m in ((18, (1, 5, 8), 4), (28, (5, 9, 12), 4), (30, (7, 13), 3), (7, (1, 2, 3), 6)):
            self.assertEqual(X.exact_shares(n, t, m)["verdict"], "REAL")
            self.assertEqual(X.float_shares(n, t, m)["verdict"], "REAL")

    def test_the_other_verdicts_are_reachable(self):
        # one turn, m = 2: a single plane cannot hold a triangle unless N = 3
        self.assertEqual(X.exact_shares(7, (1,), 2)["verdict"], "NO_SOLUTION")
        self.assertEqual(X.float_shares(7, (1,), 2)["verdict"], "NO_SOLUTION")
        self.assertEqual(X.exact_shares(3, (1,), 2)["verdict"], "REAL")
        # 12 (1, 3, 5) is a chain with shares 1/3 each; with a fourth turn the old chain
        # is met again with a share of exactly zero
        e = X.exact_shares(12, (1, 3, 5, 6), 4)
        self.assertEqual(e["verdict"], "ZERO_SHARE")
        self.assertEqual(e["shares"][3], 0)
        self.assertEqual(X.float_shares(12, (1, 3, 5, 6), 4)["verdict"], "ZERO_SHARE")


class TestExact(unittest.TestCase):
    def test_18_in_q_cos_2pi_9(self):
        e = OUT["exact"]["18_1_5_8"]
        self.assertTrue(e["every_share_in_the_smallest_basis"])
        self.assertEqual([s["form"] for s in e["shares"]],
                         ["-1/3 + 2/3*c + 2/3*c^2", "7/3 - 2/3*c - 8/3*c^2", "-1 + 2*c^2"])
        self.assertEqual(e["half_squared_edge"]["form"], "-1/2 + 2*c^2")
        c = math.cos(2 * math.pi / 9)
        self.assertAlmostEqual(e["shares"][0]["value"], -1 / 3 + 2 / 3 * c + 2 / 3 * c * c, places=12)
        self.assertTrue(e["gaps_1_to_m_exactly_equal"] and e["shares_sum_to_1_exactly"])

    def test_28_needs_sqrt7_on_two_turns(self):
        self.assertFalse(S["shares_28_1_8_13_all_in_Q_cos_2pi_7"])
        self.assertEqual(S["shares_28_1_8_13_in_Q_cos_2pi_7_by_turn"],
                         [[1, False], [8, True], [13, False]])
        self.assertEqual(S["shares_28_5_9_12_in_Q_cos_2pi_7_by_turn"],
                         [[5, False], [9, False], [12, True]])
        e = OUT["exact"]["28_1_8_13"]
        self.assertEqual(e["shares"][1]["form"], "-1*c + 2*c^2")
        self.assertEqual(e["half_squared_edge"]["form"], "1 - 1*c")
        c, r = math.cos(2 * math.pi / 7), 7 ** 0.5
        self.assertAlmostEqual(e["shares"][0]["value"],
                               1 / 2 + c / 2 - c * c + r / 14 + 5 * r * c / 14 - 3 * r * c * c / 7,
                               places=12)
        for s in e["shares"]:
            self.assertIsNotNone(s["form"])                 # every share has a closed form

    def test_sqrt7(self):
        k = Field(28)
        r = X.sqrt7(k)
        self.assertEqual(r * r, 7)
        self.assertAlmostEqual(float(r.approx()), 7 ** 0.5, places=12)

    def test_galois(self):
        g = OUT["galois"]
        for k in ("28_1_8_13_under_5", "28_1_8_13_under_3", "18_1_5_8_under_5",
                  "18_1_5_8_under_7", "30_1_11_under_7"):
            self.assertTrue(g[k]["equal"], k)
        self.assertEqual(g["28_1_8_13_under_5"]["turn_goes_to"], [[1, 5], [8, 12], [13, 9]])
        self.assertTrue(all(g["on_the_generators"].values()))
        self.assertTrue(g["unit_13_on_28_1_8_13"]["fixes_c"]
                        and g["unit_13_on_28_1_8_13"]["sqrt7_goes_to_minus"])
        self.assertEqual(S["galois_image_equals_solved_shares"], {"members_checked": 5, "equal": 5})

    def test_a_wrong_unit_does_not_pass_the_galois_check(self):
        k = Field(28)
        a = X.exact_shares(28, (1, 8, 13), 4, k)["shares"]
        b = X.exact_shares(28, (5, 9, 12), 4, k)["shares"]
        self.assertNotEqual(a[0].galois(3), b[0])           # 3 goes to (3, 4, 11), not here
        self.assertEqual(a[0].galois(5), b[0])


class TestCircles(unittest.TestCase):
    def test_agree_with_helix_isoclinic(self):
        self.assertTrue(S["circle_distances_agree_with_helix_isoclinic"])
        c = OUT["circles"]
        self.assertEqual(c["18_1_5_8"]["distinct_distances_deg"], [49.254243, 69.677963])
        self.assertEqual(c["28_1_8_13"]["distinct_distances_deg"], [46.20966])
        self.assertEqual(c["28_5_9_12"]["distinct_distances_deg"], [69.090317])

    def test_closed_forms(self):
        c = OUT["circles"]
        forms = {k: sorted({r["cos_squared"]["form"] for r in v["by_step"]}) for k, v in c.items()}
        self.assertEqual(forms, {"18_1_5_8": ["4 - 2*c - 4*c^2", "9 - 2*c - 12*c^2"],
                                 "28_1_8_13": ["-2 - 6*c + 16*c^2"],
                                 "28_5_9_12": ["16 - 8*c - 28*c^2"]})
        c9, c7 = math.cos(2 * math.pi / 9), math.cos(2 * math.pi / 7)
        for val, deg in ((4 - 2 * c9 - 4 * c9 * c9, 69.677963), (9 - 2 * c9 - 12 * c9 * c9, 49.254243),
                         (-2 - 6 * c7 + 16 * c7 * c7, 46.20966), (16 - 8 * c7 - 28 * c7 * c7, 69.090317)):
            self.assertAlmostEqual(math.degrees(math.acos(math.sqrt(val))), deg, places=5)
        self.assertEqual(S["circle_cos_squared_in_the_smallest_basis"],
                         {"18_1_5_8": True, "28_1_8_13": True, "28_5_9_12": True})


class TestRun(unittest.TestCase):
    def test_run_file_schema_and_void_slots(self):
        import tempfile
        with tempfile.TemporaryDirectory() as t:
            p = emlib.write_run(X.NAME, 0, {"tol": X.TOL}, OUT, None, None, None, None,
                                out_dir=t)
            rec = emlib.load_run(p)
            self.assertIsNone(rec["human_interest"])
            self.assertIsNone(rec["ai_interest"])
            self.assertEqual(rec["outputs"]["summary"]["chains_read"], 782)

    def test_logged_run_matches(self):
        rec = emlib.load_run(emlib.ROOT / "runs/conjugate_map_seed0.json.gz")
        self.assertEqual(rec["outputs"]["summary"], S)
        self.assertIsNone(rec["human_interest"])
        self.assertIsNone(rec["ai_interest"])


if __name__ == "__main__":
    unittest.main()
