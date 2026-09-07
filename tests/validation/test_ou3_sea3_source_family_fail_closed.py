from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

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
        # The universal conditional chain does not claim a materialized
        # physical source family. Fail-closed deployment is the distinct gate.
        self.assertTrue(d["P3_CONDITIONAL_SEA3_PASS"])
        self.assertEqual(d["P3_CONDITIONAL_SEA3_FAIL_REASONS"], [])
        self.assertFalse(d["P3_DEPLOYMENT_PASS"])
        self.assertEqual(
            d["P3_DEPLOYMENT_FAIL_REASONS"],
            ["physical SEA0->SEA3 left inclusion remains open"],
        )
        self.assertFalse(d["global_physical_deployment_left_inclusion_closed_here"])


if __name__ == "__main__":
    unittest.main()
