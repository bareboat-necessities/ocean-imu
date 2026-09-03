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
        self.assertTrue(self.d["proof_mechanism_not_part_of_canonical_definition"])
        self.assertTrue(self.d["only_this_module_may_promote_P3_for_P4"])
        self.assertEqual(self.d["useful_gate"], 1.0e-18)

    def test_refined_p2_is_ready_but_current_p3_has_not_consumed_it(self):
        self.assertTrue(self.d["P2_interface"]["timing_pass"])
        self.assertTrue(self.d["P2_interface"]["correlation_ready"])
        self.assertTrue(self.d["P2_interface"]["ready_for_canonical_P3"])
        self.assertFalse(self.d["P3_CANONICAL_PASS"])
        self.assertFalse(self.d["P4_MAY_CONSUME_P3"])
        reasons = self.d["P3_CANONICAL_FAIL_REASONS"]
        self.assertTrue(any("has not consumed" in x for x in reasons))
        self.assertTrue(any("source-history" in x for x in reasons))

    def test_both_H_and_A_margins_are_part_of_same_verdict(self):
        self.assertIn("H", self.d["mode_margins"])
        self.assertIn("A", self.d["mode_margins"])
        self.assertGreater(self.d["mode_margins"]["H"], 0.0)
        self.assertGreater(self.d["mode_margins"]["A"], 0.0)
        self.assertLess(self.d["worst_H_A_margin"], self.d["useful_gate"])

    def test_prerequisite_success_cannot_be_relabelled_as_canonical_pass(self):
        self.assertFalse(self.d["P3_CANONICAL_PASS"])
        required = self.d["required_properties"]
        self.assertTrue(required["same_source_history_for_process_covariance_measurement"])
        self.assertTrue(required["time_varying_tuner_over_word"])
        self.assertTrue(required["interleaved_accelerometer_and_S_measurements"])
        self.assertTrue(required["H_mode"])
        self.assertTrue(required["A_mode"])

    def test_equivalent_full_matrix_proof_can_satisfy_frozen_semantic_gate(self):
        p2 = G.P2I.build()
        version = p2["P3_required_correlation_interface_version"]
        cand = {
            "qualification": "SYNTHETIC_EQUIVALENT_FULL_MATRIX_P3_FOR_GATE_TEST",
            "P2_correlation_interface_consumed": True,
            "P2_correlation_interface_version": version,
            "process_covariance_measurement_bounds_same_source_history": True,
            "independent_cartesian_tau_sigma_RS_extrema_used": False,
            "source_generated_not_trajectory_fit": True,
            "trajectory_replay_used": False,
            "filter_changed": False,
            "zero_lever_arm_branch": True,
            "dormant_transparent_vibration_guard_branch": True,
            "time_varying_tuner_over_word_covered": True,
            "interleaved_accelerometer_and_S_measurements_covered": True,
            "modes": {
                "H": {"relative_Riccati_injection_margin_lower": 2.0e-18},
                "A": {"relative_Riccati_injection_margin_lower": 3.0e-18},
            },
        }
        d = G.build(p3_candidate=cand)
        self.assertEqual(G.validate(d), [])
        self.assertTrue(d["P3_CANONICAL_PASS"])
        self.assertTrue(d["P4_MAY_CONSUME_P3"])
        self.assertEqual(d["worst_H_A_margin"], 2.0e-18)
        self.assertEqual(
            d["semantic_coverage"]["time_varying_tuner_proof_advertisement"],
            "semantic",
        )
        self.assertEqual(
            d["semantic_coverage"]["interleaved_measurement_proof_advertisement"],
            "semantic",
        )

    def test_retained_mechanism_is_only_a_compatibility_advertisement(self):
        self.assertEqual(
            G._semantic_coverage({}, "semantic", "legacy"),
            (False, None),
        )
        self.assertEqual(
            G._semantic_coverage({"legacy": True}, "semantic", "legacy"),
            (True, "retained_mechanism"),
        )
        self.assertEqual(
            G._semantic_coverage({"semantic": True}, "semantic", "legacy"),
            (True, "semantic"),
        )


if __name__ == "__main__":
    unittest.main()
