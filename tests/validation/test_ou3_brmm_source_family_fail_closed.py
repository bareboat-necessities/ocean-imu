from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_complete_source as source  # noqa: E402
import ou3_brmm_riccati_metric_p3 as gate  # noqa: E402


class BrmmSourceFamilyFailClosedTest(unittest.TestCase):
    def test_universal_chain_does_not_claim_finite_materialization(self):
        s = source.build()
        self.assertEqual(source.validate(s), [])
        self.assertTrue(s["P3_source_contract_ready"])
        self.assertFalse(s["P3_source_family_materialized"])

        d = gate.build()
        self.assertEqual(gate.validate(d), [])
        self.assertFalse(d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"])
        self.assertFalse(d["finite_source_family_materialization_required"])
        self.assertTrue(d["UNIVERSAL_COMPLETE_BRMM_CERTIFICATE_CHAIN_CLOSED"])
        self.assertTrue(d["P3_CONDITIONAL_BRMM_PASS"])
        self.assertFalse(d["P3_DEPLOYMENT_PASS"])


if __name__ == "__main__":
    unittest.main()
