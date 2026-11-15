"""Calibration ledger: every Y quote is verbatim in the file it cites. CC0. Co-authored."""
import json
import os
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class Ledger(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(ROOT, "docs", "calibration_ledger.json"), encoding="utf-8") as fh:
            self.d = json.load(fh)

    def test_y_quotes_are_verbatim(self):
        for e in self.d["episodes"]:
            for y in e["Y"]:
                with open(os.path.join(ROOT, y["file"]), encoding="utf-8") as fh:
                    self.assertIn(y["quote"], fh.read(), (e["id"], y["file"]))

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
