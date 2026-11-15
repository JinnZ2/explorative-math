import json
import tempfile
import unittest

import _path  # noqa: F401
import emlib
import nature_catalog


class Schema(unittest.TestCase):
    def test_write_run_has_all_keys_and_void_human(self):
        with tempfile.TemporaryDirectory() as d:
            p = emlib.write_run("t", 0, {"a": 1}, {"x": 1}, None, emlib.NOT_INTERPRETABLE,
                                None, None, out_dir=d)
            rec = json.loads(p.read_text())
        self.assertEqual(tuple(rec), emlib.SCHEMA)
        self.assertIsNone(rec["human_interest"])
        self.assertIsNone(rec["divergence_note"])
        self.assertEqual(rec["params"]["seed"], 0)

    def test_refuses_missing_key_or_filled_human_interest(self):
        good = {k: None for k in emlib.SCHEMA}
        good["ai_interest"] = emlib.NOT_INTERPRETABLE
        emlib.check_run(good)
        bad = dict(good); del bad["angles"]
        with self.assertRaises(emlib.SchemaError):
            emlib.check_run(bad)
        bad = dict(good, human_interest="nice")
        with self.assertRaises(emlib.SchemaError):
            emlib.check_run(bad)
        bad = dict(good, ai_interest=None)
        with self.assertRaises(emlib.SchemaError):
            emlib.check_run(bad)


class LZ76(unittest.TestCase):
    def test_known_values(self):
        self.assertEqual(emlib.lz76(""), 0)
        self.assertEqual(emlib.lz76("aaaa"), 2)          # a | aaa
        self.assertEqual(emlib.lz76("0001101001000101"), 6)  # Kaspar-Schuster example
        self.assertLess(emlib.lz76("ab" * 50), emlib.lz76("abbabaabbaababba" * 6 + "abba"))


class Catalog(unittest.TestCase):
    def test_entries_carried_and_unknown_refused(self):
        self.assertTrue(all(v["status"] == "CARRIED" for v in nature_catalog.CATALOG.values()))
        self.assertIsNone(nature_catalog.match(None))
        with self.assertRaises(KeyError):
            nature_catalog.match("not-there")


if __name__ == "__main__":
    unittest.main()
