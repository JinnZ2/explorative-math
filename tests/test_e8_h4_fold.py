"""Tests for experiments/e8_h4_fold.py. CC0. Co-authored."""
import unittest

import _path  # noqa: F401
import e8_h4_fold as E
import emlib
from exact import PHI, ZPhi, ONE, ZERO


class Pieces(unittest.TestCase):
    def test_gram_is_the_e8_shape_and_is_computed(self):
        g = E.e8_gram()
        shape = E.diagram_shape(g)
        self.assertTrue(shape["simply_laced_cartan"])
        self.assertEqual((shape["n_edges"], shape["branch_node_arms"]), (7, [1, 2, 4]))
        self.assertEqual(g[0][4], 0)          # alpha_1 and phi*alpha_1 are orthogonal
        self.assertEqual(g[0][1], 0)          # the 5-edge leaves no edge between alpha_1, alpha_2

    def test_roots(self):
        rs = E.roots(E.e8_gram())
        self.assertEqual(len(rs), 240)
        self.assertTrue(all(tuple(-c for c in r) in set(rs) for r in rs))
        a2 = E.roots([[2, -1], [-1, 2]])
        self.assertEqual(len(a2), 6)

    def test_h4_reflection_is_two_e8_reflections(self):
        x = (ONE, PHI, ZERO, ZPhi(2, -1))
        y = E.h4_reflect(x, 1)
        self.assertEqual(E.h4_reflect(y, 1), x)                    # an involution
        self.assertEqual(E.h4_form(y, y), E.h4_form(x, x))         # an isometry, exactly

    def test_coxeter_element_has_order_30_exactly(self):
        x = start = (ONE, ZERO, ZERO, ZERO)
        for k in range(1, 31):
            x = E.coxeter_h4(x)
            if x == start:
                break
        self.assertEqual(k, 30)


class LoggedRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = emlib.load_run(emlib.RUNS / "e8_h4_fold_seed0.json")
        cls.o = cls.rec["outputs"]

    def test_observed(self):
        # OBSERVED on the committed run, pinned so a change shows up
        self.assertIsNone(self.rec["ai_interest"])
        e, f, c, ch = (self.o[k] for k in ("e8", "fold", "coxeter_plane", "chain"))
        self.assertTrue(e["is_e8_shape"] and e["same_h_and_exponents"])
        self.assertEqual((e["n_roots"], e["coxeter_number"]), (240, 30))
        self.assertEqual(f["squared_lengths_text"], ["2", "2 + 2*phi"])
        self.assertTrue(f["two_sets_of_120"])
        self.assertTrue(f["squared_length_ratio_is_phi_squared"])
        self.assertTrue(f["larger_set_is_phi_times_the_smaller_exactly"])
        self.assertTrue(f["smaller_set"]["reads_as_600_cell_vertices"])
        self.assertEqual((c["n_rings"], c["every_ring_has_30"]), (8, True))
        self.assertEqual(ch["exponents_on_the_fold_space"], [1, 11, 19, 29])
        self.assertEqual(ch["exponents_on_the_conjugate_space"], [7, 13, 17, 23])
        self.assertTrue(ch["together_equal_e8_exponents"])
        self.assertEqual((ch["n_orbits_of_roots"], ch["orbit_sizes"]), (8, [30]))
        for key in ("chain_1_11", "chain_7_13"):
            self.assertEqual(ch[key]["verdict"], "match")
            self.assertEqual(ch[key]["orbits_that_are_the_logged_chain"], 2)
            self.assertEqual(ch[key]["orbits"], 8)


if __name__ == "__main__":
    unittest.main()
