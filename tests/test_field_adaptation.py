"""Tests for the rules, performance, inverse readout and field_adaptation. CC0. Co-authored."""
import ast
import json
import math
import unittest
from pathlib import Path

import _path  # noqa: F401
import field_adaptation as V
from field import hexsheet as H


class Rules(unittest.TestCase):
    def test_saturating_follows_linear_when_small_and_levels_off(self):
        lin = H.preferred_spacing(0.01, 0.01)
        sat = H.preferred_spacing(0.01, 0.01, "saturating", 1.5)
        self.assertAlmostEqual(lin, sat, places=8)
        self.assertLess(H.preferred_spacing(100.0, 0.01, "saturating", 1.5), 1.0 + 0.01 * 1.5 + 1e-12)
        with self.assertRaises(ValueError):
            H.preferred_spacing(1.0, 0.01, "saturating")
        with self.assertRaises(ValueError):
            H.preferred_spacing(1.0, 0.01, "no_such_rule")

    def test_inverse_of_each_rule(self):
        for rule, s in (("linear", None), ("saturating", 1.5)):
            for v in (-2.0, -0.3, 0.0, 0.7, 2.5):
                sp = H.preferred_spacing(v, 0.01, rule, s)
                self.assertAlmostEqual(H.field_from_spacing(sp, 0.01, rule, s), v, places=9)

    def test_outside_the_saturating_range_is_void_not_clipped(self):
        # clearly past the level-off value 1 + alpha * s = 1.015 (the exact edge is
        # a float boundary and is not probed)
        self.assertIsNone(H.field_from_spacing(1.0151, 0.01, "saturating", 1.5))
        self.assertIsNone(H.field_from_spacing(0.98, 0.01, "saturating", 1.5))
        self.assertIsNone(H.field_from_spacing(1.1, 0.01, "saturating", 1.5))


class Performance(unittest.TestCase):
    def test_zero_for_a_sheet_in_balance(self):
        sheet = H.cells(2)
        es = H.rest_lengths(sheet, H.uniform_field(sheet, 0.0), 0.01)
        p = H.performance(sheet, es, {c: H.ideal_xy(c) for c in sheet})
        self.assertLess(p["energy"], 1e-25)
        self.assertLess(p["max_net_pull"], 1e-12)

    def test_uniform_nonzero_field_strains_the_ideal_reference(self):
        sheet = H.cells(2)
        es = H.rest_lengths(sheet, H.uniform_field(sheet, 2.0), 0.01)
        p = H.performance(sheet, es, {c: H.ideal_xy(c) for c in sheet})
        self.assertAlmostEqual(p["max_abs_strain"], 0.02 / 1.02, places=9)
        self.assertGreater(p["energy"], 0.0)


class InverseReadout(unittest.TestCase):
    def test_signature_has_no_field_argument(self):
        # AST walk: the readout cannot be handed the field
        src = Path(H.__file__).read_text()
        fn = next(n for n in ast.walk(ast.parse(src))
                  if isinstance(n, ast.FunctionDef) and n.name == "inverse_readout")
        names = [a.arg for a in fn.args.args]
        self.assertNotIn("field", names)
        used = {n.id for n in ast.walk(fn) if isinstance(n, ast.Name)}
        self.assertNotIn("field", used)

    def test_recovers_a_gradient_on_interior_cells(self):
        sheet = H.cells(3)
        field = H.gradient_field(sheet, 0.5)
        pos, _ = H.relax(sheet, field, 0.01, 800)
        rec = H.inverse_readout(sheet, pos, 0.01)
        s = set(sheet)
        inner = [c for c in sheet
                 if all((c[0] + dq, c[1] + dr) in s for dq, dr in H.NEIGHBOURS)]
        for c in inner:
            self.assertAlmostEqual(rec[c], field[c], delta=0.02)


class FieldAdaptation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = V.run(radius=3, sweeps=600, slopes=(0.0, 0.5))

    def test_zero_field_control_both_sheets_equal(self):
        for rule in H.RULES:
            p = self.out["rules"][rule]["control_zero_field"]["performance"]
            self.assertTrue(all(p["equal"].values()), rule)
            self.assertFalse(any(p["adapted_lower"].values()), rule)

    def test_under_the_field_adapted_is_lower_on_every_measure(self):
        for rule in H.RULES:
            p = self.out["rules"][rule]["main"]["performance"]
            self.assertTrue(all(p["adapted_lower"].values()), rule)
            self.assertEqual(self.out["rules"][rule]["slopes_where_ideal_reference_is_not_worse"],
                             [0.0])

    def test_readout_logs_interior_and_edge_and_where_it_fails(self):
        for rule in H.RULES:
            r = self.out["rules"][rule]["main"]["inverse_readout"]
            for part in ("interior", "edge", "all"):
                self.assertIn("correlation", r[part])
            self.assertEqual(r["interior"]["n_cells"] + r["edge"]["n_cells"], r["all"]["n_cells"])
            self.assertEqual(r["where_it_fails"]["n"],
                             r["where_it_fails"]["interior"] + r["where_it_fails"]["edge"])
            self.assertEqual(len(r["cells"]), self.out["bound"]["n_cells"])
        self.assertIn("inverse_readout_linear_rule_assumed", self.out["rules"]["saturating"]["main"])
        self.assertNotIn("inverse_readout_linear_rule_assumed", self.out["rules"]["linear"]["main"])

    def test_naming(self):
        self.assertEqual(self.out["naming"]["ideal"], "the reference sheet, not a target")
        self.assertIn("relative to the environment", self.out["naming"]["optimum"])

        def keys(node):
            if isinstance(node, dict):
                for k, v in node.items():
                    yield k
                    yield from keys(v)
            elif isinstance(node, list):
                for v in node:
                    yield from keys(v)
        for k in set(keys(self.out)):
            self.assertNotIn("error", k.lower())
            self.assertNotIn("target", k.lower())
        doc = (H.__doc__ + __import__("field").__doc__).lower()
        self.assertIn("not a target", doc.replace("\n", " ").replace("  ", " "))
        self.assertIn("relative to the environment", " ".join(doc.split()))

    def test_sweep_must_hold_zero_and_the_main_slope(self):
        with self.assertRaises(ValueError):
            V.run(radius=2, sweeps=10, slopes=(0.5,))

    def test_logged_runs(self):
        # OBSERVED on the committed run, pinned so a change shows up
        rec = V.emlib.load_run(V.emlib.RUNS / "field_adaptation_seed0.json")
        o = rec["outputs"]
        self.assertIsNone(rec["human_interest"])
        lin, sat = o["rules"]["linear"], o["rules"]["saturating"]
        self.assertEqual(lin["main"]["inverse_readout"]["where_it_fails"]["n"], 0)
        self.assertEqual(sat["main"]["inverse_readout"]["where_it_fails"]["n"], 8)
        self.assertEqual(sat["main"]["inverse_readout"]["where_it_fails"]["interior"], 0)
        self.assertEqual(sat["main"]["inverse_readout_linear_rule_assumed"]["where_it_fails"]["n"], 30)
        for r in (lin, sat):
            self.assertEqual(r["slopes_where_ideal_reference_is_not_worse"], [0.0])
            self.assertTrue(all(r["control_zero_field"]["performance"]["equal"].values()))
        small = V.emlib.load_run(V.emlib.RUNS / "field_adaptation_small_seed0.json")
        self.assertTrue(small["params"]["small"])
        self.assertEqual(set(small["outputs"]), set(o))
        json.dumps(o)
        self.assertTrue(math.isfinite(lin["main"]["performance"]["adapted_sheet"]["energy"]))


if __name__ == "__main__":
    unittest.main()
