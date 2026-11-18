"""Opt-in options on three_gap_sum, angle_histogram, rotational_symmetry_bound and emlib.write_run. CC0. Co-authored.

Two kinds of test:
1. DEFAULT UNCHANGED. The functions as they were at commit 3162c99 are copied below, unedited, as frozen
   references. With the new options left off, the new functions must return exactly what the frozen ones return
   (same keys, same values) on many inputs, and write_run must write the same bytes as an independent json.dumps.
2. OPTIONS DO WHAT THEY SAY, on the inputs that the adversarial run recorded (A01, A02, A04, A05, A07, A19).
Nothing here pins the old behaviour on inputs where it is the finding (bare NaN, dropped angles): the
options are off by default, and the tests for the defaults use finite, in-range data.
"""
import gzip
import json
import math
import random
import tempfile
from decimal import Decimal, getcontext
import unittest
from pathlib import Path
import _path  # noqa: F401
import emlib
from toolkit import primitives as P
from toolkit import probes as Q


# ---- frozen references: toolkit/probes.py at 3162c99, bodies unedited --------------------------------
def old_three_gap_sum(lengths, tol=1e-6):
    ls = sorted(float(x) for x in lengths)
    if len(ls) != 3:
        return {"n_lengths": len(ls), "verdict": "not_applicable", "residual": None,
                "tol": tol}
    r = ls[2] - (ls[0] + ls[1])
    return {"n_lengths": 3, "verdict": "pass" if abs(r) <= tol else "fail",
            "residual": r, "tol": tol}


def old_angle_histogram(angles_deg, width=10.0, lo=0.0, hi=180.0):
    nb = int(round((hi - lo) / width))
    counts = [0] * nb
    for a in angles_deg:
        b = int((a - lo) // width)
        if 0 <= b < nb:
            counts[b] += 1
        elif a == hi:
            counts[-1] += 1
    return {"lo": lo, "hi": hi, "width": width, "counts": counts}


def old_rotational_symmetry_bound(points, k_max=12, rel_tol=1e-6):
    return {"n_points": len(points), "k_tested": [2, k_max], "rel_tol": rel_tol,
            "abs_tol": rel_tol * Q._scale_of(points) if points else None,
            "about": "centroid"}


def hexagon(t=0.0):
    return [[t + math.cos(math.pi * i / 3), t + math.sin(math.pi * i / 3)] for i in range(6)]


def a07_points():
    ring = [[1e6 * math.cos(2 * math.pi * i / 12), 1e6 * math.sin(2 * math.pi * i / 12)] for i in range(12)]
    return [[0.3, 0.0], [-0.1, 0.25], [-0.15, -0.2]] + ring


class DefaultsUnchanged(unittest.TestCase):
    def test_three_gap_sum(self):
        rnd = random.Random(1)
        for _ in range(2000):
            n = rnd.choice([0, 1, 2, 3, 3, 3, 4])
            ls = [rnd.uniform(1e-3, 5) for _ in range(n)]
            if n == 3 and rnd.random() < 0.5:
                ls[2] = ls[0] + ls[1] + rnd.choice([0.0, 1e-9, 1e-5])
            for tol in (1e-6, 1e-3):
                self.assertEqual(Q.three_gap_sum(ls, tol), old_three_gap_sum(ls, tol))
                self.assertEqual(Q.three_gap_sum(ls), old_three_gap_sum(ls))

    def test_angle_histogram(self):
        rnd = random.Random(2)
        for _ in range(500):
            width = rnd.choice([5.0, 10.0, 15.0, 7.0, 30.0, 0.5])
            lo, hi = rnd.choice([(0.0, 180.0), (0.0, 360.0), (-90.0, 90.0)])
            angs = [rnd.uniform(lo - 20, hi + 20) for _ in range(rnd.randint(0, 40))] + [lo, hi]
            self.assertEqual(Q.angle_histogram(angs, width, lo, hi), old_angle_histogram(angs, width, lo, hi))
        self.assertEqual(Q.angle_histogram([0, 5, 95, 180], width=90), old_angle_histogram([0, 5, 95, 180], width=90))

    def test_rotational_symmetry_bound(self):
        rnd = random.Random(3)
        for _ in range(300):
            pts = [[rnd.uniform(-3, 3), rnd.uniform(-3, 3)] for _ in range(rnd.randint(0, 12))]
            for k, t in ((12, 1e-6), (24, 1e-9)):
                self.assertEqual(Q.rotational_symmetry_bound(pts, k, t), old_rotational_symmetry_bound(pts, k, t))
        self.assertEqual(Q.rotational_symmetry_bound(hexagon()), old_rotational_symmetry_bound(hexagon()))

    def test_write_run_default_bytes_match_independent_dumps(self):
        outputs = {"x": [1, 2.5, None, "a"], "y": {"z": 0.1}}
        with tempfile.TemporaryDirectory() as d:
            p = emlib.write_run("optin_default", 3, {"k": 1}, outputs, None, None, None, None, out_dir=d)
            got = Path(p).read_bytes()
            q = emlib.write_run("optin_strict", 3, {"k": 1}, outputs, None, None, None, None, out_dir=d, strict_json=True)
            got_strict = Path(q).read_bytes()
        rec = {"params": {"k": 1, "seed": 3}, "outputs": outputs, "angles": None, "human_interest": None,
               "ai_interest": None, "compressibility": None, "nature_match": None, "divergence_note": None}
        want = (json.dumps(rec, indent=1, sort_keys=False) + "\n").encode("utf-8")
        self.assertEqual(got, want)
        self.assertEqual(got_strict, want)


class ThreeGapRel(unittest.TestCase):
    def test_a01_tiny_lengths(self):                       # A01: absolute 1e-6 passes a 57 % residual
        self.assertEqual(Q.three_gap_sum([1e-9, 2e-9, 7e-9])["verdict"], "pass")
        r = Q.three_gap_sum([1e-9, 2e-9, 7e-9], rel_tol=1e-6)
        self.assertEqual(r["verdict"], "fail")
        self.assertEqual((r["rel_tol"], r["tol"]), (1e-6, 1e-6))
        self.assertAlmostEqual(r["tol_used"], 7e-15, places=25)

    def test_a02_golden_triple_at_1e12(self):              # A02: one float step fails the absolute test
        getcontext().prec = 60                               # each value exact, then rounded to a float on its own
        phi = (1 + Decimal(5).sqrt()) / 2
        scale = Decimal(10) ** 12
        big = [float(scale), float(scale * phi), float(scale * phi * phi)]
        self.assertEqual(Q.three_gap_sum(big)["verdict"], "fail")
        self.assertEqual(Q.three_gap_sum(big, rel_tol=1e-9)["verdict"], "pass")

    def test_scale_invariance_of_the_relative_test(self):
        rnd = random.Random(4)
        for _ in range(300):
            a, b = rnd.uniform(0.1, 3), rnd.uniform(0.1, 3)
            for k in range(-12, 13, 3):
                c = 10.0 ** k
                ok = Q.three_gap_sum([a * c, b * c, (a + b) * c], rel_tol=1e-9)
                bad = Q.three_gap_sum([a * c, b * c, (a + b) * c * (1 + 1e-3)], rel_tol=1e-9)
                self.assertEqual((ok["verdict"], bad["verdict"]), ("pass", "fail"), (a, b, k))

    def test_other_counts_stay_not_applicable_and_record_rel_tol(self):
        r = Q.three_gap_sum([1.0, 2.0], rel_tol=1e-6)
        self.assertEqual((r["verdict"], r["residual"], r["rel_tol"]), ("not_applicable", None, 1e-6))


class HistogramAccount(unittest.TestCase):
    def test_a04_nothing_vanishes_without_a_count(self):
        h = Q.angle_histogram([-5.0, 185.0, 360.0, 90.0], width=10.0, lo=0.0, hi=180.0, account=True)
        self.assertEqual((h["n_input"], h["n_binned"], h["below_lo"], h["above_hi"]), (4, 1, 1, 2))
        self.assertEqual(h["n_input"], h["n_binned"] + h["below_lo"] + h["above_hi"])
        self.assertEqual(h["counts"], old_angle_histogram([-5.0, 185.0, 360.0, 90.0])["counts"])

    def test_a05_bins_cover_more_than_hi(self):
        h = Q.angle_histogram([181.0], width=7.0, lo=0.0, hi=180.0, account=True)
        self.assertEqual(h["bins_cover"], [0.0, 182.0])
        self.assertEqual((h["n_binned"], h["binned_above_hi"], h["above_hi"]), (1, 1, 0))
        self.assertEqual(sum(h["counts"]), 1)                # the count itself is not changed

    def test_identity_on_random_input(self):
        rnd = random.Random(5)
        for _ in range(300):
            width = rnd.choice([5.0, 7.0, 10.0, 15.0])
            angs = [rnd.uniform(-30, 400) for _ in range(rnd.randint(0, 50))] + [0.0, 360.0]
            h = Q.angle_histogram(angs, width, 0.0, 360.0, account=True)
            self.assertEqual(h["n_input"], h["n_binned"] + h["below_lo"] + h["above_hi"])
            self.assertEqual(h["counts"], old_angle_histogram(angs, width, 0.0, 360.0)["counts"])
            if (360.0 / width) == int(360.0 / width):
                self.assertEqual(h["binned_above_hi"], 0)

    def test_dividing_width_covers_exactly_hi(self):
        h = Q.angle_histogram([1.0], width=15.0, lo=0.0, hi=360.0, account=True)
        self.assertEqual(h["bins_cover"], [0.0, 360.0])


class BoundSpacing(unittest.TestCase):
    def test_a07_cluster_beside_a_far_ring(self):
        b = Q.rotational_symmetry_bound(a07_points(), 12, 1e-6, with_spacing=True)
        self.assertGreater(b["abs_tol"], b["min_pair_distance"])
        self.assertEqual(b["n_points_within_tol"], 3)
        self.assertEqual(Q.rotational_symmetry(a07_points(), 12), [2, 3, 4, 6, 12])   # the list itself is unchanged

    def test_clean_set_reports_none_within_tol(self):
        b = Q.rotational_symmetry_bound(hexagon(), 12, 1e-6, with_spacing=True)
        self.assertEqual(b["n_points_within_tol"], 0)
        self.assertAlmostEqual(b["min_pair_distance"], 1.0)

    def test_small_sets(self):
        self.assertIsNone(Q.rotational_symmetry_bound([], with_spacing=True)["min_pair_distance"])
        self.assertIsNone(Q.rotational_symmetry_bound([[0.0, 0.0]], with_spacing=True)["min_pair_distance"])
        self.assertEqual(Q.rotational_symmetry_bound([[0.0, 0.0], [0.0, 0.0]], with_spacing=True)["n_points_within_tol"], 2)

    def test_spacing_keys_are_a_superset(self):
        base = Q.rotational_symmetry_bound(hexagon())
        more = Q.rotational_symmetry_bound(hexagon(), with_spacing=True)
        self.assertTrue(all(more[k] == v for k, v in base.items()))


class StrictJson(unittest.TestCase):
    def rec(self, out_dir, **kw):
        return emlib.write_run("optin_nan", 0, {"k": 1}, {"x": float("nan"), "y": [1.0, float("inf")],
                                                         "z": {"w": -float("inf")}}, None, None, None, None,
                               out_dir=out_dir, **kw)

    def test_a19_strict_refuses_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(emlib.SchemaError) as cm:
                self.rec(d, strict_json=True)
            msg = str(cm.exception)
            for want in ("outputs/x", "outputs/y/1", "outputs/z/w"):
                self.assertIn(want, msg)
            self.assertEqual(list(Path(d).iterdir()), [])

    def test_nonfinite_paths(self):
        self.assertEqual(emlib.nonfinite_paths({"a": [1.0, float("nan")], "b": {"c": float("inf")}, "d": 2}),
                         ["a/1", "b/c"])
        self.assertEqual(emlib.nonfinite_paths({"a": [1, 2.5, None, "x"]}), [])
        self.assertEqual(len(emlib.nonfinite_paths([float("nan")] * 50, limit=7)), 7)

    def test_strict_output_parses_under_a_strict_parser(self):
        def refuse(t):
            raise ValueError(t)
        with tempfile.TemporaryDirectory() as d:
            p = emlib.write_run("optin_ok", 0, {"k": 1}, {"x": [1.5, 2.0]}, None, None, None, None, out_dir=d,
                                strict_json=True)
            text = Path(p).read_text()
        json.loads(text, parse_constant=refuse)

    def test_strict_with_gzip_form(self):
        with tempfile.TemporaryDirectory() as d:
            p = emlib.write_run("optin_gz", 0, {"k": 1}, {"x": [0.5] * 10}, None, None, None, None, out_dir=d,
                                gz=True, strict_json=True)
            self.assertTrue(str(p).endswith(".json.gz"))
            self.assertEqual(json.loads(gzip.open(p, "rt").read())["outputs"]["x"], [0.5] * 10)

    def test_off_by_default_the_writer_still_writes_what_it_always_did(self):
        with tempfile.TemporaryDirectory() as d:
            p = self.rec(d)                                  # default: no refusal (committed files depend on this)
            self.assertTrue(Path(p).exists())


if __name__ == "__main__":
    unittest.main()
