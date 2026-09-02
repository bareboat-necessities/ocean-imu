#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p3_canonical_gate as G


class CanonicalP3GateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build()

    def test_definition_is_frozen_and_gate_is_unchanged(self):
        self.assertEqual(G.validate(self.d), [])
        self.assertTrue(self.d["canonical_definition_frozen"])
        self.assertTrue(self.d["only_this_module_may_promote_P3_for_P4"])
        self.assertEqual(self.d["useful_gate"], 1.0e-18)

    def test_current_state_cannot_promote_without_refined_p2_interface(self):
        self.assertTrue(self.d["P2_interface"]["timing_pass"])
        self.assertFalse(self.d["P2_interface"]["correlation_ready"])
        self.assertFalse(self.d["P2_interface"]["ready_for_canonical_P3"])
        self.assertFalse(self.d["P3_CANONICAL_PASS"])
        self.assertFalse(self.d["P4_MAY_CONSUME_P3"])
        self.assertTrue(any("P2" in x for x in self.d["P3_CANONICAL_FAIL_REASONS"]))

    def test_both_H_and_A_margins_are_part_of_same_verdict(self):
        self.assertIn("H", self.d["mode_margins"])
        self.assertIn("A", self.d["mode_margins"])
        self.assertGreater(self.d["mode_margins"]["H"], 0.0)
        self.assertGreater(self.d["mode_margins"]["A"], 0.0)
        self.assertLess(self.d["worst_H_A_margin"], self.d["useful_gate"])

    def test_prerequisite_success_cannot_be_relabelled_as_canonical_pass(self):
        self.assertFalse(self.d["P3_CANONICAL_PASS"])
        required = self.d["required_properties"]
        self.assertTrue(required["time_varying_tuner_over_word"])
        self.assertTrue(required["interleaved_accelerometer_and_S_measurements"])
        self.assertTrue(required["H_mode"])
        self.assertTrue(required["A_mode"])


if __name__ == "__main__":
    unittest.main()
