from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_riccati_metric_p3 as gate  # noqa: E402


class BrmmFullWordSourceIdentityTest(unittest.TestCase):
    def test_same_complete_source_drives_both_modes(self):
        d = gate.build()
        self.assertEqual(gate.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_BRMM_NORMAL_LIVE_WORD")
        self.assertTrue(d["same_complete_BRMM_execution_continues_across_H_to_A"])
        self.assertTrue(d["same_primitive_root_drives_entire_execution"])
        self.assertFalse(d["same_three_second_same_mode_word_used_for_H18_and_A21"])
        self.assertTrue(d["H_to_A_is_separate_dimension_changing_hybrid_event"])
        self.assertFalse(d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"])
        self.assertTrue(d["universal_complete_BRMM_certificate_chain_used_instead_of_finite_materialization"])
        self.assertTrue(d["P3_FULL_WORD_ENCLOSED"])
        self.assertTrue(d["P3_FULL_MATRIX_COMPARISON_CLOSED"])
        self.assertTrue(d["P3_CONDITIONAL_BRMM_PASS"])
        self.assertFalse(d["P3_DEPLOYMENT_PASS"])


if __name__ == "__main__":
    unittest.main()
