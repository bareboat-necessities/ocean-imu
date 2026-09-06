#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_a21_prior_free_completion as A21


class Sea3A21PriorFreeCompletionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = A21.build()
        cls.failures = A21.validate(cls.d)

    def test_hybrid_full_matrix_certificate_closes(self):
        d = self.d
        self.assertEqual([], self.failures)
        self.assertEqual("COMPLETE_SEA3_NORMAL_LIVE_WORD", d["canonical_source"])
        self.assertEqual(1e-18, d["useful_gate"])
        self.assertTrue(d["H18_full_matrix_margin_inherited_before_release"])
        self.assertTrue(d["H18_full_18x18_interval_LDLT_consumed"])
        self.assertTrue(d["full_21x21_Omega_minus_delta_P_closed"])
        self.assertTrue(d["A21_prior_free_completion_closed"])
        self.assertTrue(d["exact_direct_sum_full_matrix_hybrid_proof_used"])
        self.assertGreater(d["H18_worst_LDLT_pivot_lower"], 0.0)

    def test_shipping_hold_forces_H18_to_finish_before_A_release(self):
        d = self.d
        self.assertTrue(d["shipping_outer_mag_refinement_hold_consumed"])
        self.assertTrue(d["H_to_A_is_separate_dimension_changing_hybrid_event"])
        self.assertTrue(d["H18_word_finishes_before_A_release"])
        self.assertGreaterEqual(
            d["shipping_H_mode_minimum_duration_before_A_release_s"],
            d["canonical_H18_word_horizon_s"],
        )
        self.assertTrue(d["held_ba_cross_covariances_zero_at_release"])
        self.assertEqual(
            d["held_ba_covariance_variance_at_release"],
            d["release_ba_floor_variance"],
        )

    def test_first_active_ba_prediction_closes_new_directions(self):
        d = self.d
        self.assertTrue(d["first_active_prediction_block_diagonal_from_zero_release_cross_covariance"])
        self.assertTrue(d["first_active_ba_block_strictly_positive"])
        self.assertGreater(d["shipping_active_ba_Q_lambda_min_lower"], 0.0)
        self.assertGreater(d["first_active_ba_M_delta_margin_lower"], 0.0)
        self.assertGreater(
            d["shipping_active_ba_Q_lambda_min_lower"],
            d["delta_times_release_ba_variance_upper"],
        )
        self.assertEqual(1.0, d["phi_b_squared_upper"])
        self.assertEqual([], d["prediction_source_parity_failures"])
        self.assertTrue(all(d["prediction_source_parity"].values()))

    def test_no_disallowed_A21_shortcuts(self):
        d = self.d
        self.assertEqual(
            "ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION",
            d["paper_active_bias_route"],
        )
        self.assertTrue(d["A21_detectability_certificate_retained_as_independent_support"])
        self.assertFalse(d["eta9_point_packet_shortcut_used"])
        self.assertFalse(d["full_A21_prior_free_D_inverse_identity_used"])
        self.assertFalse(d["finite_full_A21_linear_estimator_constructed"])
        self.assertFalse(d["blockwise_minimum_contraction_used"])
        self.assertFalse(d["old_one_step_Euclidean_full_state_Q_min_used"])
        self.assertFalse(d["scalar_beta_contraction_used"])

    def test_source_and_live_measurement_semantics_are_retained(self):
        d = self.d
        self.assertTrue(d["event_algebra_preserves_margin_after_first_active_prediction"])
        self.assertTrue(d["actual_applied_SpectralMSE_R_S_retained_through_inherited_H18_margin"])
        self.assertTrue(d["all_Normal_Live_accelerometer_updates_retained"])
        self.assertFalse(d["accelerometer_rejection_after_certified_Normal_Live_allowed"])
        for key in (
            "source_family_replaced",
            "trajectory_replay_used",
            "independent_tau_sigma_RS_source_created",
            "filter_changed_for_A21_proof",
            "P3_promoted",
        ):
            self.assertFalse(d[key], key)


if __name__ == "__main__":
    unittest.main()
