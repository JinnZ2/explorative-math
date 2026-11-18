"""Checks docs/probe_provenance_2026-11-18.json against the repo it cites. CC0. Co-authored.
A quote from a repo file must occur in that file (whitespace-normalised). Quotes from outside the repo
(RFC 8259) cannot be checked offline and are skipped by name. Facts marked COMPUTED are recomputed here
where the note states a number."""
import json
import math
import re
import unittest
from pathlib import Path
import _path  # noqa: F401
from toolkit import probes as Q

ROOT = Path(__file__).resolve().parents[1]
PROV = json.loads((ROOT / "docs" / "probe_provenance_2026-11-18.json").read_text(encoding="utf-8"))


def norm(t):
    return re.sub(r"\s+", " ", t).strip()


class Provenance(unittest.TestCase):
    def test_every_case_present_and_classed(self):
        ids = [f"A{i:02d}" for i in range(1, 21)]
        self.assertEqual(sorted(PROV["cases"]), ids)
        for i, c in PROV["cases"].items():
            self.assertIn(c["audited"], ("HOLDS", "DEFECT", "LIMIT", "UNDEFINED"), i)
            self.assertTrue(c["basis"], i)
            for b in c["basis"]:
                self.assertIn(b["kind"], PROV["kinds"], i)
                self.assertIn(b["strength"], PROV["strengths"], i)

    def test_repo_quotes_occur_in_their_files(self):
        checked = 0
        for i, c in PROV["cases"].items():
            for b in c["basis"]:
                if not b["quote"] or b["strength"] == "EXTERNAL":
                    continue
                text = norm((ROOT / b["source"]).read_text(encoding="utf-8"))
                self.assertIn(norm(b["quote"]), text, f"{i}: quote not found in {b['source']}")
                checked += 1
        self.assertGreaterEqual(checked, 20)

    def test_external_quotes_are_dated(self):
        for i, c in PROV["cases"].items():
            for b in c["basis"]:
                if b["strength"] == "EXTERNAL":
                    self.assertIn("read 2026-11-18", b["source"], i)

    def test_computed_facts(self):
        self.assertAlmostEqual(math.ulp(1e12), 1.220703125e-04)
        self.assertAlmostEqual(math.ulp(2.618033988749895e12), 4.8828125e-04)
        h = Q.angle_histogram([181.0], width=7.0, lo=0.0, hi=180.0)      # A05
        self.assertEqual((h["hi"], len(h["counts"]), sum(h["counts"])), (180.0, 26, 1))
        self.assertGreater(h["lo"] + len(h["counts"]) * h["width"], h["hi"])
        self.assertEqual((-1e-15) % 360.0, 360.0)                        # A06 b
        self.assertEqual(math.atan2(0.0, 0.0), 0.0)                      # A06 a


if __name__ == "__main__":
    unittest.main()
