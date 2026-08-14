from __future__ import annotations

import re
import unittest

from ontseq_platform.smoke import synthetic_sam_text


def _query_length(cigar: str) -> int:
    return sum(
        int(length)
        for length, operation in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if operation in {"M", "I", "S", "=", "X"}
    )


class SyntheticSmokeFixtureTests(unittest.TestCase):
    def test_sam_fixture_is_structurally_consistent_and_synthetic(self) -> None:
        text = synthetic_sam_text()
        records = [line.split("\t") for line in text.splitlines() if not line.startswith("@")]

        self.assertEqual(len(records), 24)
        self.assertEqual(sum(fields[5] == "5000M200D5000M" for fields in records), 12)
        self.assertTrue(all(len(fields[9]) >= 10_000 for fields in records))
        for fields in records:
            self.assertTrue(fields[0].startswith("SYNTH_"))
            self.assertEqual(len(fields[9]), _query_length(fields[5]))
            self.assertEqual(len(fields[10]), len(fields[9]))
            self.assertIn("RG:Z:SYNTHETIC_RG", fields)
        self.assertNotIn("patient", text.lower())


if __name__ == "__main__":
    unittest.main()
