from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_sea3_riccati_metric_p3 as gate  # noqa: E402


class Sea3FullWordSourceIdentityTest(unittest.TestCase):
    def test_same_complete_source_drives_both_modes(self):
        d = gate.build()
        self.assertEqual(gate.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(d["same_complete_SEA3_word_used_for_H18_and_A21"])
        self.assertTrue(d["same_xs_lambda_drives_entire_word"])
        self.assertFalse(d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"])
        self.assertFalse(d["P3_FULL_WORD_ENCLOSED"])
        self.assertFalse(d["P3_FULL_MATRIX_COMPARISON_CLOSED"])
        self.assertFalse(d["P3_CANONICAL_PASS"])


if __name__ == "__main__":
    unittest.main()
