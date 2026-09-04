#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p3_p2_v1_full_state_join as J
import ou3_source_reachable_matrix_p3 as BASE


class FullStatePrecisionJoinTests(unittest.TestCase):
    def test_join_factor_and_canonical_gate_are_frozen(self):
        self.assertEqual(J.JOIN_FACTOR, 0.5)
        self.assertEqual(BASE.MIN_USEFUL_DELTA, 1.0e-18)

    def test_scalar_two_block_precision_inequality(self):
        # For J=[[a,c],[c,b]] > 0, verify numerically that
        # J^-1 - 0.5*diag(1/a,1/b) is positive semidefinite.  This is the
        # scalar-block instance of the operator inequality used by the join.
        cases = (
            (2.0, 3.0, 0.0),
            (2.0, 3.0, 1.0),
            (5.0, 7.0, 4.0),
            (1.0, 1.0, 0.999),
        )
        for a, b, c in cases:
            det = a * b - c * c
            self.assertGreater(det, 0.0)
            x11 = b / det - 0.5 / a
            x22 = a / det - 0.5 / b
            x12 = -c / det
            self.assertGreaterEqual(x11, -1e-14)
            self.assertGreaterEqual(x22, -1e-14)
            self.assertGreaterEqual(x11 * x22 - x12 * x12, -1e-14)

    def test_conditional_bias_floor_is_strict_and_source_independent(self):
        domain = J.json.loads(J.DEFAULT_DOMAIN.read_text(encoding="utf-8"))
        blocks = J._common_blocks(domain)
        self.assertTrue(blocks["vector_PE_recurrence_quantifier_consumed"])
        self.assertEqual(
            blocks["final_pe_recurrence_s"],
            float(domain["normal_live"]["vector_pe_recurrence_window_s"]),
        )
        for mode in ("H", "A"):
            row = J._conditional_bias_floor(mode, domain, blocks)
            self.assertGreater(row["conditional_measurement_attenuation_lower"], 0.0)
            self.assertGreater(row["attitude_conditional_posterior_lower"], 0.0)
            self.assertGreater(row["gyro_bias_conditional_posterior_lower"], 0.0)
            self.assertTrue(row["translation_columns_conditioned_known"])
            self.assertTrue(row["fresh_final_prediction_modes_only"])
            if mode == "A":
                self.assertGreater(row["accel_bias_conditional_posterior_lower"], 0.0)
            else:
                self.assertIsNone(row["accel_bias_conditional_posterior_lower"])

    def test_history_bias_upper_uses_final_pe_window_not_whole_word(self):
        domain = J.json.loads(J.DEFAULT_DOMAIN.read_text(encoding="utf-8"))
        blocks = J._common_blocks(domain)
        summary = {
            "all_statistics_from_one_legal_P2_history": True,
            "independent_global_source_extrema_used": False,
            "history_duration_s": [10.0, 10.0],
            "pseudo_update_cadence_s": [0.1, 0.1],
            "q_c_upper": 1.0,
        }
        for mode in ("H", "A"):
            row = J._history_bias_upper(summary, mode, domain, blocks)
            self.assertTrue(row["final_PE_selected_from_recurrence_window"])
            self.assertFalse(row["whole_word_horizon_used_for_post_PE_propagation"])
            self.assertTrue(row["vector_PE_recurrence_quantifier_consumed"])
            self.assertEqual(
                row["post_vector_PE_to_endpoint_horizon_s_upper"],
                J.BASE.up(blocks["final_pe_recurrence_s"]),
            )
            self.assertLess(
                row["post_vector_PE_to_endpoint_horizon_s_upper"],
                row["whole_word_horizon_s_upper"],
            )

    def test_validation_is_fail_closed_on_numeric_flags(self):
        d = {
            "schema": J.SCHEMA,
            "qualification": "OU3_P3_P2_V1_SOURCE_COMPLETE_HA_PRECISION_BLOCK_JOIN",
            "source_generated_not_trajectory_fit": True,
            "trajectory_replay_used": False,
            "filter_changed": False,
            "declared_domain_changed": False,
            "zero_lever_arm_branch": True,
            "dormant_transparent_vibration_guard_branch": True,
            "P2_correlation_interface_consumed": True,
            "P2_correlation_interface_version": J.CORR.INTERFACE_VERSION,
            "process_covariance_measurement_bounds_same_source_history": True,
            "independent_cartesian_tau_sigma_RS_extrema_used": False,
            "independent_cartesian_tau_sigma_R_S_extrema_used": False,
            "time_varying_tuner_over_word_covered": True,
            "interleaved_accelerometer_and_S_measurements_covered": True,
            "finite_clock_13_26_stage_language_covered": True,
            "frozen_clock_absorbing_hold_covered": True,
            "translation_full_matrix_samplewise_measurements_consumed": True,
            "translation_whole_word_lift_charged_again": False,
            "attitude_bias_fresh_final_prediction_modes_used": True,
            "conditional_precision_block_theorem_used": True,
            "same_history_bias_upper_evaluated_before_uniform_envelope": True,
            "vector_PE_word_language_certificate_consumed": True,
            "final_vector_PE_recurrence_window_used_for_attitude_bias_propagation": True,
            "whole_covariance_word_used_for_post_PE_propagation": False,
            "final_vector_PE_recurrence_window_s": 1.0,
            "precision_block_join_factor": 0.5,
            "useful_gate": 1.0e-18,
            "modes": {
                "H": {
                    "relative_Riccati_injection_margin_lower": 2.0e-18,
                    "useful_margin_established": True,
                    "same_history_bias_covariance_upper": {
                        "final_PE_selected_from_recurrence_window": True,
                        "whole_word_horizon_used_for_post_PE_propagation": False,
                    },
                },
                "A": {
                    "relative_Riccati_injection_margin_lower": 3.0e-18,
                    "useful_margin_established": True,
                    "same_history_bias_covariance_upper": {
                        "final_PE_selected_from_recurrence_window": True,
                        "whole_word_horizon_used_for_post_PE_propagation": False,
                    },
                },
            },
            "P3_PRODUCER_NUMERIC_PASS": True,
            "P3_CANONICAL_PROMOTED": False,
            "P4_PROMOTED": False,
            "failures": [],
        }
        self.assertEqual(J.validate(d), [])
        d["P3_PRODUCER_NUMERIC_PASS"] = False
        self.assertIn(
            "producer numeric pass flag does not match H/A margins",
            J.validate(d),
        )
        d["P3_PRODUCER_NUMERIC_PASS"] = True
        d["whole_covariance_word_used_for_post_PE_propagation"] = True
        self.assertIn(
            "whole_covariance_word_used_for_post_PE_propagation is not false",
            J.validate(d),
        )


if __name__ == "__main__":
    unittest.main()
