from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))
import ou3_p4_complete_brmm_accelerometer_operation_coordinate as acccoord  # noqa: E402
import ou3_p4_complete_brmm_finite_angle_information as finfo  # noqa: E402
import ou3_p4_complete_brmm_phi_differential_metric as diffmetric  # noqa: E402
import ou3_p4_complete_brmm_vector_remainder_geometry as remgeom  # noqa: E402
import ou3_p4_complete_word_endpoint_transport as endpoint  # noqa: E402
import ou3_brmm_riccati_metric_p4 as mod  # noqa: E402


class BrmmFiniteStateP4Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()
        cls.ep = endpoint.build()
        cls.dm = diffmetric.build()
        cls.i = finfo.build()
        cls.c = acccoord.build()
        cls.g = remgeom.build()
        print("P4_FINITE_ANGLE_CELLS", [
            (x["attitude_angle_deg"], x["full_H18_information_lambda_min_lower"])
            for x in cls.i["candidate_cells"]
        ])
        print("P4_ENDPOINT_MASTER", {
            "qualification": cls.ep["qualification"],
            "actual_RS_in_suffix": cls.ep["actual_RS_regularization_enters_every_applicable_suffix"],
            "endpoint_closed": cls.ep["source_uniform_master_endpoint_domination_closed"],
        })

    def test_p4_uses_paper_finite_state_endpoint_and_prefix_architecture(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(endpoint.validate(self.ep), [])
        self.assertEqual(
            self.d["canonical_P4_architecture"],
            "FINITE_STATE_COMPLETE_BRMM_QUADRATIC_ENDPOINT_AND_PREFIX",
        )
        self.assertEqual(self.d["canonical_source"], "COMPLETE_BRMM_NORMAL_LIVE_WORD")
        self.assertEqual(self.d["paper_Lyapunov_function"], "V(e,zeta)=e^T M(zeta)e")
        self.assertTrue(self.d["finite_state_endpoint_dissipation_required"])
        self.assertTrue(self.d["finite_state_prefix_gain_required"])
        self.assertTrue(self.d["prefix_chart_and_source_domain_retention_required"])
        self.assertGreaterEqual(self.d["outer_angle_rad"], 0.80)
        self.assertEqual(self.d["candidate_angles_deg"], [30.0, 25.0, 20.0, 15.0])

    def test_frozen_conditional_p3_is_consumed_without_modification(self):
        self.assertTrue(self.d["P3_CONDITIONAL_BRMM_PASS_consumed"])
        self.assertTrue(self.d["P3_frozen_not_modified"])
        self.assertEqual(self.d["P3_H18_delta_consumed"], 1.0e-18)
        self.assertEqual(self.d["P3_A21_delta_consumed"], 1.0e-18)
        self.assertFalse(self.d["P3_DEPLOYMENT_PASS_consumed_as_if_closed"])

    def test_primary_motion_target_does_not_relabel_full_state_P4(self):
        self.assertEqual(self.d["primary_proof_target"], mod.MOTION.TARGET)
        motion = self.d["bounded_bias_motion_contract"]
        self.assertEqual(mod.MOTION.validate(motion), [])
        self.assertFalse(motion["bias_error_convergence_required"])
        self.assertFalse(motion["zero_error_floor_required"])
        self.assertFalse(self.d["P4_MOTION_PASS"])
        self.assertFalse(self.d["P5_MOTION_MAY_START"])
        self.assertFalse(self.d["P4_CANONICAL_PASS"])

    def test_endpoint_master_retains_complete_word_and_actual_rs(self):
        for key in (
            "same_complete_BRMM_word_required",
            "same_frontend_tuner_covariance_history_required",
            "all_due_S_updates_with_actual_applied_RS_required",
            "all_valid_accelerometer_updates_required",
            "all_process_Q_floor_reset_events_required",
            "asynchronous_vector_events_required",
            "H_to_A_rectangular_hybrid_event_required",
            "full_prediction_F_Eaw_rows_retained",
            "full_epsilon_aw_retained",
            "endpoint_master_object_emitted",
            "actual_RS_regularization_enters_suffix_maps",
        ):
            self.assertTrue(self.d[key], key)
        self.assertEqual(self.d["H_to_A_homogeneous_lift"], "[I18;0]")
        self.assertTrue(self.d["H_to_A_held_ba_error_retained_as_separate_forcing"])
        self.assertTrue(self.d["H_to_A_covariance_floor_retained_as_separate_metric_event"])
        self.assertIn("B_W", self.d["joint_accelerometer_suffix_operator"])
        self.assertFalse(self.d["independent_RS_schedule_used"])
        self.assertFalse(self.d["point_word_rho_used_to_promote"])
        self.assertFalse(self.d["longer_point_window_optimization_used_to_promote"])

    def test_full_state_finite_taub_joint_sector_bridge_is_retained_but_open(self):
        self.assertTrue(self.d["full_state_joint_sector_master_available"])
        self.assertTrue(self.d["A21_finite_taub_full_matrix_P3_required"])
        self.assertTrue(self.d["A21_finite_taub_detectability_consumed"])
        self.assertTrue(self.d["joint_sector_all_21_state_cross_terms_retained"])
        self.assertTrue(self.d["joint_sector_actual_RS_provenance_preserved"])
        self.assertFalse(self.d["source_uniform_same_history_joint_sector_closed"])
        self.assertFalse(self.d["source_uniform_full_augmented_LDLT_closed"])

    def test_differential_ad_is_machinery_not_replacement_p4(self):
        self.assertEqual(diffmetric.validate(self.dm), [])
        self.assertTrue(self.d["differential_AD_used_only_for_finite_map_enclosure"])
        self.assertFalse(self.d["differential_pullback_used_as_replacement_P4"])
        self.assertFalse(self.d["differential_finite_distance_bridge_closed"])
        self.assertTrue(self.d["same_source_omega_h_tau_prediction_required"])
        self.assertTrue(self.d["prediction_independent_F_forbidden"])
        self.assertTrue(self.d["same_P_H_R_cell_required_for_Joseph"])
        self.assertTrue(self.d["independent_K_forbidden"])
        self.assertTrue(self.d["A21_bias_projection_generalized_Jacobian_machinery_retained"])
        self.assertTrue(self.d["outward_interval_AD_machinery_retained"])

    def test_accelerometer_operation_coordinate_retains_full_shift_and_rs(self):
        self.assertEqual(acccoord.validate(self.c), [])
        self.assertEqual(self.c["canonical_source"], "COMPLETE_BRMM_NORMAL_LIVE_WORD")
        self.assertTrue(self.c["aw_error_exactly_linear_in_accelerometer_operation_coordinate"])
        self.assertTrue(self.c["accelerometer_bias_error_exactly_linear"])
        self.assertGreater(self.c["latent_aw_nonlinear_eta_coefficient"], 0.0)
        self.assertTrue(self.c["mixed_aw_shipping_tangent_remainder_retained"])
        self.assertFalse(self.c["nonlinear_Phi_storage_is_original_metric_isometry"])
        self.assertTrue(self.c["actual_RS_regularizer_not_removed_by_coordinate_change"])
        self.assertFalse(self.c["source_family_replaced"])

    def test_vector_remainder_geometry_remains_nonpromoting_input(self):
        self.assertEqual(remgeom.validate(self.g), [])
        self.assertEqual(self.g["canonical_source"], "COMPLETE_BRMM_NORMAL_LIVE_WORD")
        self.assertTrue(self.g["mixed_aw_shipping_tangent_remainder_retained"])
        self.assertTrue(self.g["all_due_S_updates_and_actual_RS_remain_in_complete_word"])
        self.assertFalse(self.g["packet_count_multiplier_used"])
        self.assertFalse(self.g["standalone_eta_disturbance_budget_used"])

    def test_finite_angle_information_keeps_declared_candidates_without_promotion(self):
        self.assertEqual(finfo.validate(self.i), [])
        self.assertEqual(self.i["canonical_source"], "COMPLETE_BRMM_NORMAL_LIVE_WORD")
        self.assertTrue(self.i["P3_frozen_not_modified"])
        self.assertTrue(self.i["actual_applied_SpectralMSE_R_S_consumed_through_frozen_H18_information"])
        self.assertTrue(self.i["all_due_S_updates_remain_in_complete_word"])
        self.assertTrue(self.i["all_valid_accelerometer_updates_remain_in_complete_word"])
        self.assertFalse(self.i["selected_PE_or_four_S_replace_complete_word"])
        self.assertEqual(self.i["widest_information_cell_deg"], 30.0)
        self.assertTrue(self.i["information_headroom_closed"])
        self.assertFalse(self.i["P4_promoted_here"])

    def test_only_canonical_blockers_are_endpoint_and_prefix_theorems(self):
        self.assertFalse(self.d["source_uniform_joint_BW_epsilon_enclosure_closed"])
        self.assertFalse(self.d["source_uniform_r_word_enclosure_closed"])
        self.assertFalse(self.d["source_uniform_master_endpoint_domination_closed"])
        self.assertFalse(self.d["source_uniform_prefix_gain_closed"])
        self.assertFalse(self.d["source_uniform_prefix_domain_retention_closed"])
        self.assertFalse(self.d["P4_CANONICAL_PASS"])
        self.assertFalse(self.d["P5_MAY_START"])
        self.assertEqual(2, len(self.d["P4_CANONICAL_FAIL_REASONS"]))
        self.assertTrue(any("endpoint" in x.lower() for x in self.d["P4_CANONICAL_FAIL_REASONS"]))
        self.assertTrue(any("prefix" in x.lower() for x in self.d["P4_CANONICAL_FAIL_REASONS"]))
        for key in (
            "packet_count_remainder_budget_used",
            "packetwise_remainder_norm_sum_used",
            "state_elimination_used",
            "a_w_Schur_final_certificate_used",
            "correction_radius_claim_used",
            "inverse_metric_floor_claim_used",
        ):
            self.assertFalse(self.d[key], key)


if __name__ == "__main__":
    unittest.main()
