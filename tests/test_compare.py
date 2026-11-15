import json
import tempfile
import unittest

import _path  # noqa: F401
from compare import cross as X


class FieldVerdicts(unittest.TestCase):
    def test_numeric_exact_agree_and_differ_keep_diffs(self):
        self.assertEqual(X.field("x", 1.0, 1.0)["verdict"], X.AGREE)
        d = X.field("x", 1.0, 1.0 + 1e-16 * 3)
        self.assertEqual(d["verdict"], X.DIFFER)   # no tolerance filtering
        self.assertGreater(d["abs_diff"], 0.0)
        self.assertIn("rel_diff", d)

    def test_only_in(self):
        self.assertEqual(X.field("x", 1, None, True, False)["verdict"], X.ONLY_A)
        self.assertEqual(X.field("x", None, 1, False, True)["verdict"], X.ONLY_B)
        with self.assertRaises(ValueError):
            X.field("x", None, None, False, False)

    def test_vector_fields_length_mismatch(self):
        v = [f["verdict"] for f in X.vector_fields("v", [1, 2, 3], [1, 2])]
        self.assertEqual(v, [X.AGREE, X.AGREE, X.ONLY_A])


class Records(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = X.run_all()

    def test_union_schema_and_void_observer_slots(self):
        for recs in self.results.values():
            for r in recs:
                self.assertEqual(tuple(r), X.SCHEMA_UNION)
                for k in ("human_interest", "ai_interest", "compressibility",
                          "divergence_note", "nature_match", "angles"):
                    self.assertIsNone(r[k], (r["id"], k))
                self.assertEqual(r["status"], "OBSERVED")

    def _tally(self, name):
        t = {v: 0 for v in X.VERDICTS}
        for r in self.results[name]:
            for v in X.VERDICTS:
                t[v] += r["outputs"]["tally"][v]
        return t

    def test_observed_seeds_do_not_align(self):
        t = self._tally("seed_alignment")
        self.assertEqual(t[X.AGREE], 0)

    def test_observed_shared_stream_and_shared_rotations_agree(self):
        for name in ("shared_stream_sampler", "rotate_project"):
            t = self._tally(name)
            self.assertEqual(t[X.DIFFER], 0, name)
            self.assertGreater(t[X.AGREE], 0, name)

    def test_observed_phi_compositions_differ(self):
        self.assertGreater(self._tally("phi_composition")[X.DIFFER], 0)

    def test_write_all_jsonl(self):
        with tempfile.TemporaryDirectory() as d:
            paths = X.write_all(d)
            for p in paths:
                for line in p.read_text().splitlines():
                    json.loads(line)
        self.assertEqual(len(paths), len(X.COMPARISONS))


if __name__ == "__main__":
    unittest.main()
