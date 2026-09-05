from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_complete_source as source  # noqa: E402
import ou3_sea3_riccati_metric_p3 as gate  # noqa: E402


class Sea3SourceFamilyFailClosedTest(unittest.TestCase):
    def test_contract_readiness_does_not_equal_materialization(self):
        s = source.build()
        self.assertEqual(source.validate(s), [])
        self.assertTrue(s["P3_source_contract_ready"])
        self.assertFalse(s["P3_source_family_materialized"])

        d = gate.build()
        self.assertEqual(gate.validate(d), [])
        self.assertFalse(d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"])
        self.assertFalse(d["P3_CANONICAL_PASS"])
        self.assertTrue(d["P3_CANONICAL_FAIL_REASONS"])


if __name__ == "__main__":
    unittest.main()
