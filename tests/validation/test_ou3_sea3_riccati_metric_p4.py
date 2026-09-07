from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))
import ou3_p4_complete_sea3_accelerometer_operation_coordinate as acccoord  # noqa: E402
import ou3_p4_complete_sea3_finite_angle_information as finfo  # noqa: E402
import ou3_p4_complete_sea3_phi_differential_metric as diffmetric  # noqa: E402
import ou3_p4_complete_sea3_vector_remainder_geometry as remgeom  # noqa: E402
import ou3_sea3_riccati_metric_p4 as mod  # noqa: E402


class Sea3DifferentialP4Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()
        cls.dm = diffmetric.build()
        cls.i = finfo.build()
        cls.c = acccoord.build()
        cls.g = remgeom.build()
        print("P4_FINITE_ANGLE_CELLS", [
            (x["attitude_angle_deg"], x["full_H18_information_lambda_min_lower"])
            for x in cls.i["candidate_cells"]
        ])
        print("P4_DIFFERENTIAL_METRIC", {
            "type": cls.dm["metric_type"],
            "det_DPhi": cls.dm["Phi_jacobian_determinant_exact"],
            "full_rank_H18": cls.dm["H18_full_rank"],
            "full_rank_A21": cls.dm["A21_full_rank"],
            "actual_RS_required": cls.dm["all_due_S_updates_with_actual_applied_RS_required"],
        })

    def test_p4_uses_complete_sea3_full_state_differential_architecture(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(diffmetric.validate(self.dm), [])
        self.assertEqual(
            self.d["canonical_P4_architecture"],
            "FULL_STATE_COMPLETE_SEA3_DIFFERENTIAL_PULLBACK",
        )
        self.assertEqual(self.d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertGreaterEqual(self.d["outer_angle_rad"], 0.80)
        self.assertEqual(self.d["candidate_angles_deg"], [30.0, 25.0, 20.0, 15.0])

    def test_differential_metric_is_full_rank_and_reduces_to_frozen_p3(self):
        self.assertTrue(self.d["P3_CONDITIONAL_SEA3_PASS_consumed"])
        self.assertTrue(self.d["P3_frozen_not_modified"])
        self.assertEqual(self.d["P3_H18_delta_consumed"], 1.0e-18)
        self.assertEqual(self.d["P3_A21_delta_consumed"], 1.0e-18)
        self.assertTrue(self.d["all_active_states_retained"])
        self.assertTrue(self.d["H18_full_rank"])
        self.assertTrue(self.d["A21_full_rank"])
        self.assertEqual(self.d["Phi_jacobian_determinant_exact"], 1.0)
        self.assertTrue(self.d["metric_reduces_exactly_to_P3_at_zero_error"])
        self.assertFalse(self.d["finite_Phi_storage_used_as_Lyapunov_function"])
        self.assertFalse(self.d["finite_raw_endpoint_storage_used_as_P4_certificate"])

    def test_complete_shipping_word_and_actual_rs_are_mandatory(self):
        for key in (
            "same_complete_SEA3_word_required",
            "same_frontend_tuner_covariance_history_required",
            "all_due_S_updates_with_actual_applied_RS_required",
            "all_valid_accelerometer_updates_required",
            "all_process_Q_floor_reset_events_required",
            "H_to_A_rectangular_differential_event_required",
        ):
            self.assertTrue(self.d[key], key)
        self.assertFalse(self.d["independent_RS_schedule_used"])
        self.assertFalse(self.d["point_word_rho_used_to_promote"])
        self.assertFalse(self.d["longer_point_window_optimization_used_to_promote"])

    def test_accelerometer_operation_coordinate_retains_full_shift_and_rs(self):
        self.assertEqual(acccoord.validate(self.c), [])
        self.assertEqual(self.c["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.c["aw_error_exactly_linear_in_accelerometer_operation_coordinate"])
        self.assertTrue(self.c["accelerometer_bias_error_exactly_linear"])
        self.assertGreater(self.c["latent_aw_nonlinear_eta_coefficient"], 0.0)
        self.assertTrue(self.c["mixed_aw_shipping_tangent_remainder_retained"])
        self.assertFalse(self.c["nonlinear_Phi_storage_is_original_metric_isometry"])
        self.assertTrue(self.c["actual_RS_regularizer_not_removed_by_coordinate_change"])
        self.assertFalse(self.c["source_family_replaced"])

    def test_vector_remainder_geometry_remains_nonpromoting_input(self):
        self.assertEqual(remgeom.validate(self.g), [])
        self.assertEqual(self.g["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.g["mixed_aw_shipping_tangent_remainder_retained"])
        self.assertTrue(self.g["all_due_S_updates_and_actual_RS_remain_in_complete_word"])
        self.assertFalse(self.g["packet_count_multiplier_used"])
        self.assertFalse(self.g["standalone_eta_disturbance_budget_used"])

    def test_finite_angle_information_keeps_declared_candidates_without_promotion(self):
        self.assertEqual(finfo.validate(self.i), [])
        self.assertEqual(self.i["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.i["P3_frozen_not_modified"])
        self.assertTrue(self.i["actual_applied_SpectralMSE_R_S_consumed_through_frozen_H18_information"])
        self.assertTrue(self.i["all_due_S_updates_remain_in_complete_word"])
        self.assertTrue(self.i["all_valid_accelerometer_updates_remain_in_complete_word"])
        self.assertFalse(self.i["selected_PE_or_four_S_replace_complete_word"])
        self.assertEqual(self.i["widest_information_cell_deg"], 30.0)
        self.assertTrue(self.i["information_headroom_closed"])
        self.assertFalse(self.i["P4_promoted_here"])

    def test_only_canonical_blockers_are_full_jacobian_and_differential_ldlt(self):
        self.assertFalse(self.d["source_uniform_complete_word_Jacobian_enclosed"])
        self.assertFalse(self.d["source_uniform_pullback_differential_contraction_closed"])
        self.assertFalse(self.d["P4_CANONICAL_PASS"])
        self.assertFalse(self.d["P5_MAY_START"])
        self.assertEqual(2, len(self.d["P4_CANONICAL_FAIL_REASONS"]))
        self.assertTrue(any("jacobian" in x.lower() for x in self.d["P4_CANONICAL_FAIL_REASONS"]))
        self.assertTrue(any("ldlt" in x.lower() or "contraction" in x.lower()
                            for x in self.d["P4_CANONICAL_FAIL_REASONS"]))
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
