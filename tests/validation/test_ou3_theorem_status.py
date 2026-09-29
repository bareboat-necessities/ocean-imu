import json,sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.theorem_status import status_report

class TheoremStatusTests(unittest.TestCase):
    def test_single_architecture_is_fail_closed(self):
        r=status_report()
        self.assertEqual(r["principal_assumptions"],["MARINE MOTION","IMU BIAS","MAGNETIC SERVICE"])
        self.assertFalse(r["parallel_no_magnetometer_stability_path"])
        self.assertFalse(r["theorem_closed"])
        self.assertFalse(r["regional_practical_stability_claimed"])
        self.assertIn("H18-to-A21 release",r["proof_path"])
    def test_full_state_claims_require_more_than_partial_lin_and_service(self):
        o=status_report()["obligations"]
        self.assertTrue(o["LIN_endpoint_matrix_path_action"])
        self.assertTrue(o["complete_word_covariance_energy_identity"])
        self.assertTrue(o["historical_AG_readout_covariance_implication"])
        self.assertFalse(o["uniform_historical_AG_readout_action"])
        self.assertTrue(o["source_covariance_projection_guard"])
        self.assertFalse(o["certified_tail_prefix_retention"])
        for key in ("full_state_magnetic_information_lifting",
                    "constructive_full_A21_mu_rho_enclosure",
                    "source_uniform_A21_linear_dissipativity",
                    "constructive_fixed_coordinate_A21_mu_enclosure"):
            self.assertFalse(o[key])
    def test_regime_and_pivot_algebra_does_not_promote_stability(self):
        r=status_report()
        self.assertEqual(r["physical_regimes"],["STILL","TRANSITION","MOVING"])
        self.assertFalse(r["shipping_mode_switch_enabled"])
        self.assertTrue(r["exact_rest_detector_obstruction"]["identical_rest_motion_IMU_histories"])
        o=r["obligations"]
        for key in ("stationary_observability_structure", "stationary_gyro_information_bound",
                    "finite_transition_prefix_composition", "conditional_two_group_six_column_bound",
                    "singular_floor_to_all_six_historical_pivots"):
            self.assertTrue(o[key])
        for key in ("literal_reset_inverse_nonexpansion", "conditional_inverse_frame_gyro_transport_bound",
                    "same_prediction_cell_actual_row_factorization",
                    "same_cell_reset_pulled_field_geometry_implication",
                    "quiet_zero_residual_nominal_historical_covariance_upper",
                    "quiet_zero_residual_nominal_homogeneous_linear_loss_exists"):
            self.assertTrue(o[key])
        for key in ("stationary_compatible_class_practical_stability", "sound_runtime_regime_certification",
                    "finite_regime_detection_and_transition_retention", "recurring_regime_storage_budget",
                    "temporal_margins_to_six_historical_pivots", "full_21_covariance_upper"):
            self.assertFalse(o[key])

    def test_world_frame_rows_do_not_promote_aggregate_or_aw_bounds(self):
        r=status_report()
        o=r["obligations"]
        for key in ("world_frame_historical_row_factorization", "attitude_invariant_same_cell_geometry",
                    "literal_injection_covariance_NIS_budget",
                    "quiet_nominal_six_column_floor_any_attitude_admitted_reference",
                    "nominal_attitude_column_physical_transfer_given_AW_tracking",
                    "sharp_isotropic_sync_AW_covariance_ceiling"):
            self.assertTrue(o[key])
        for key in ("uniform_AW_tracking_error_bound", "aggregate_world_frame_six_column_floor",
                    "uniform_historical_AG_readout_action"):
            self.assertFalse(o[key])
        w=r["world_frame_rows"]
        self.assertFalse(w["attitude_estimate_or_error_in_six_column_geometry"])
        self.assertTrue(w["aggregate_rows_avoid_magnetic_phase_premise"])
        self.assertTrue(w["same_cell_floor_requires_magnetic_cadence_coupling"])
        self.assertTrue(w["aw_covariance_ceiling_tight_at_sync"])
        self.assertTrue(w["uniform_storage_route_ratio_above_one_on_carried_collinear_words"])
        for key in ("same_cell_floor_from_motion_and_bias_bounds_alone",
                    "collinear_1Hz_witness_magnetic_service_admitted",
                    "same_cell_two_group_uniform_floor_established",
                    "signed_world_injection_transport_budget"):
            self.assertFalse(w[key])

    def test_moving_continuation_keeps_its_open_premises_open(self):
        r=status_report()
        o=r["obligations"]
        for key in ("nominal_signed_mean_attitude_columns",
                    "pointwise_physical_AW_tracking_premise_refuted_on_admitted_history",
                    "signed_injection_rotation_identity", "third_order_literal_reset_factor",
                    "magnetic_service_tube_monotone_field_axis",
                    "injection_free_aggregate_six_column_floor_given_nominal_window_premises"):
            self.assertTrue(o[key])
        for key in ("source_uniform_nominal_AW_window_statistics",
                    "source_uniform_half_angle_injection_transport",
                    "aggregate_world_frame_six_column_floor",
                    "practical_rho0_margin_from_six_column_floor",
                    "source_uniform_A21_linear_dissipativity"):
            self.assertFalse(o[key])
        w=r["world_frame_rows"]
        self.assertTrue(w["attitude_columns_need_only_nominal_signed_mean"])
        self.assertIn("A_tilde=I", w["injection_free_floor_premises"])

    def test_committed_status_matches_code(self):
        p=ROOT/"reports/results/ou3_stability/theorem-status.json"
        self.assertEqual(json.loads(p.read_text()),status_report())

if __name__=="__main__": unittest.main()
