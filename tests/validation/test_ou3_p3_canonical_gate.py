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
        cls.p2 = G.P2I.build()
        cls.version = cls.p2["P3_required_correlation_interface_version"]

    def candidate(self, h=2.0e-19, a=3.0e-19, *, consume=True):
        return {
            "qualification": "SYNTHETIC_SEMANTIC_P3_GATE_TEST",
            "P2_correlation_interface_consumed": consume,
            "P2_correlation_interface_version": self.version if consume else None,
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
                "H": {"relative_Riccati_injection_margin_lower": h},
                "A": {"relative_Riccati_injection_margin_lower": a},
            },
        }

    def test_definition_is_frozen_and_gate_is_unchanged(self):
        d = G.build(p3_candidate=self.candidate())
        self.assertEqual(G.validate(d), [])
        self.assertTrue(d["canonical_definition_frozen"])
        self.assertTrue(d["proof_mechanism_not_part_of_canonical_definition"])
        self.assertTrue(d["only_this_module_may_promote_P3_for_P4"])
        self.assertEqual(d["useful_gate"], 1.0e-18)

    def test_refined_p2_is_ready(self):
        d = G.build(p3_candidate=self.candidate())
        self.assertTrue(d["P2_interface"]["timing_pass"])
        self.assertTrue(d["P2_interface"]["correlation_ready"])
        self.assertTrue(d["P2_interface"]["ready_for_canonical_P3"])
        self.assertEqual(d["P2_interface"]["correlation_version"], self.version)

    def test_both_H_and_A_margins_are_part_of_same_fail_closed_verdict(self):
        d = G.build(p3_candidate=self.candidate())
        self.assertFalse(d["P3_CANONICAL_PASS"])
        self.assertFalse(d["P4_MAY_CONSUME_P3"])
        self.assertEqual(d["mode_margins"]["H"], 2.0e-19)
        self.assertEqual(d["mode_margins"]["A"], 3.0e-19)
        self.assertEqual(d["worst_H_A_margin"], 2.0e-19)
        reasons = d["P3_CANONICAL_FAIL_REASONS"]
        self.assertTrue(any("H: margin below" in x for x in reasons))
        self.assertTrue(any("A: margin below" in x for x in reasons))

    def test_missing_refined_source_interface_stays_blocked(self):
        d = G.build(p3_candidate=self.candidate(h=2e-18, a=3e-18, consume=False))
        self.assertFalse(d["P3_CANONICAL_PASS"])
        self.assertFalse(d["P4_MAY_CONSUME_P3"])
        reasons = d["P3_CANONICAL_FAIL_REASONS"]
        self.assertTrue(any("has not consumed" in x for x in reasons))
        self.assertTrue(any("required P2 correlation-interface version" in x for x in reasons))

    def test_equivalent_full_matrix_proof_can_satisfy_frozen_semantic_gate(self):
        d = G.build(p3_candidate=self.candidate(h=2.0e-18, a=3.0e-18))
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
