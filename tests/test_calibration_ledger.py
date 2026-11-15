"""Calibration ledger: every Y quote is verbatim in the file it cites. CC0. Co-authored."""
import json
import os
import subprocess
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class Ledger(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(ROOT, "docs", "calibration_ledger.json"), encoding="utf-8") as fh:
            self.d = json.load(fh)

    def _text(self, y):
        if y.get("ref"):
            return subprocess.run(["git", "-C", ROOT, "show", f"{y['ref']}:{y['file']}"],
                                  capture_output=True, text=True, check=True).stdout
        with open(os.path.join(ROOT, y["file"]), encoding="utf-8") as fh:
            return fh.read()

    def test_y_quotes_are_verbatim(self):
        for e in self.d["episodes"]:
            for y in e["Y"]:
                if "quote" in y:
                    self.assertIn(y["quote"], self._text(y), (e["id"], y["file"]))

    def test_y_counts_recompute(self):
        for e in self.d["episodes"]:
            for y in e["Y"]:
                if "count" not in y:
                    continue
                fields = [f for line in self._text(y).splitlines() if line.strip()
                          for f in json.loads(line)["outputs"].get("fields", [])]
                c = y["count"]
                self.assertEqual(len(fields), c["of"], (e["id"], y["file"]))
                self.assertEqual(sum(f["verdict"] == c["verdict"] for f in fields), c["n"],
                                 (e["id"], y["file"]))

    def test_corrections_keep_the_prior_entry(self):
        for c in self.d.get("corrections", []):
            self.assertTrue(c["was"] and c["why"] and c["by"], c["episode"])

    def test_fields_and_states(self):
        for e in self.d["episodes"]:
            for f in ("id", "H", "A", "C", "Y", "state", "relation", "reasoning"):
                self.assertIn(f, e)
            self.assertIn(e["state"], ("Y_PRESENT", "UNRESOLVED"))
            self.assertIn(e["relation"], ("SUPPORTS", "CONTRADICTS", "MIXED", "DOES_NOT_TEST"))
            if e["state"] == "UNRESOLVED":
                self.assertEqual(e["relation"], "DOES_NOT_TEST")

    def test_interest_declared(self):
        self.assertIn("party", self.d["interest"])


if __name__ == "__main__":
    unittest.main()
