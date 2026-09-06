from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p4_complete_sea3_accelerometer_operation_coordinate as acccoord  # noqa: E402
import ou3_p4_complete_sea3_finite_angle_information as finfo  # noqa: E402
import ou3_p4_complete_sea3_vector_remainder_geometry as remgeom  # noqa: E402
import ou3_p4_moving_metric_rebind as rebind  # noqa: E402
import ou3_sea3_riccati_metric_p4 as mod  # noqa: E402


class Sea3MovingRiccatiP4Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()
        cls.r = rebind.build()
        cls.i = finfo.build()
        cls.c = acccoord.build()
        cls.g = remgeom.build()
        print("P4_FINITE_ANGLE_CELLS", [
            (x["attitude_angle_deg"], x["full_H18_information_lambda_min_lower"])
            for x in cls.i["candidate_cells"]
        ])
        print("P4_ACCEL_OPERATION_COORDINATE", {
            "aw_eta": cls.c["latent_aw_nonlinear_eta_coefficient"],
            "moving_metric": cls.c["state_coordinate_transform"]["moving_metric_energy_invariant"],
            "actual_RS_retained": cls.c["actual_RS_regularizer_not_removed_by_coordinate_change"],
        })
        print("P4_VECTOR_REMAINDER_CELLS", [
            (x["attitude_angle_deg"], x["eta_squared_over_linear_tangent_squared_upper"])
            for x in cls.g["candidate_cells"]
        ])

    def test_p4_uses_same_moving_metric_and_not_800_endpoint_scan(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(rebind.validate(self.r), [])
        self.assertEqual(
            self.d["canonical_P4_architecture"],
            "NONLINEAR_WORD_IN_MOVING_SHIPPING_RICCATI_METRIC",
        )
        self.assertFalse(self.d["old_800_endpoint_signed_Joseph_scan_consumed"])
        self.assertFalse(self.d["old_terminal_source_phase_metric_attachment_consumed"])
        self.assertFalse(self.d["old_group_isotropic_P3_P4_metric_assumed"])
        self.assertGreaterEqual(self.d["outer_angle_rad"], 0.80)

    def test_exact_moving_metric_rebind_is_closed(self):
        self.assertTrue(self.d["P3_CANONICAL_PASS_consumed"])
        self.assertTrue(self.d["full_nonlinear_measurement_metric_rebind_closed"])
        self.assertFalse(self.d["exact_vector_accelerometer_congruence_rebind_pending"])
        self.assertTrue(self.d["moving_metric_coordinate_congruence_exact"])
        self.assertTrue(self.d["Joseph_nonlinear_injection_metric_closed"])
        self.assertTrue(self.r["prediction_linear_map_nonexpansive"])
        self.assertTrue(self.r["Joseph_linear_map_nonexpansive"])
        self.assertTrue(self.r["left_error_reset_exact_metric_isometry"])
        self.assertFalse(self.r["group_isotropic_metric_attachment_used"])
        self.assertFalse(self.r["endpoint_source_word_scan_used"])

    def test_accelerometer_operation_coordinate_is_exact_in_moving_metric(self):
        self.assertEqual(acccoord.validate(self.c), [])
        self.assertEqual(self.c["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.c["aw_error_exactly_linear_in_accelerometer_operation_coordinate"])
        self.assertTrue(self.c["accelerometer_bias_error_exactly_linear"])
        self.assertEqual(self.c["latent_aw_nonlinear_eta_coefficient"], 0.0)
        self.assertEqual(self.c["accelerometer_bias_nonlinear_eta_coefficient"], 0.0)
        self.assertTrue(self.c["state_coordinate_transform"]["moving_metric_energy_invariant"])
        self.assertTrue(self.c["state_coordinate_transform"]["innovation_covariance_S_invariant"])
        self.assertFalse(self.c["state_coordinate_transform"]["group_isotropic_metric_assumption_required"])
        self.assertTrue(self.c["actual_RS_regularizer_not_removed_by_coordinate_change"])
        self.assertFalse(self.c["source_family_replaced"])
        self.assertFalse(self.c["retired_endpoint_attachment_module_reintroduced"])

    def test_pure_vector_remainder_is_homogeneous_and_keeps_RS(self):
        self.assertEqual(remgeom.validate(self.g), [])
        self.assertEqual(self.g["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertEqual(self.g["accelerometer_aw_nonlinear_eta_coefficient"], 0.0)
        self.assertEqual(self.g["accelerometer_bias_nonlinear_eta_coefficient"], 0.0)
        self.assertTrue(self.g["accelerometer_eta_is_pure_force_rotation"])
        self.assertTrue(self.g["magnetometer_radial_Joseph_cancellation_consumed"])
        self.assertTrue(self.g["R_inverse_whitening_preserves_pure_vector_sector"])
        self.assertTrue(self.g["all_due_S_updates_and_actual_RS_remain_in_complete_word"])
        self.assertFalse(self.g["packet_count_multiplier_used"])
        self.assertFalse(self.g["standalone_eta_disturbance_budget_used"])
        self.assertFalse(self.g["retired_source_grid_remainder_route_used"])

    def test_complete_sea3_finite_angle_information_is_a_p4_prerequisite(self):
        self.assertEqual(finfo.validate(self.i), [])
        self.assertEqual(self.i["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.i["P3_frozen_not_modified"])
        self.assertTrue(self.i["actual_applied_SpectralMSE_R_S_consumed_through_frozen_H18_information"])
        self.assertTrue(self.i["all_due_S_updates_remain_in_complete_word"])
        self.assertTrue(self.i["all_valid_accelerometer_updates_remain_in_complete_word"])
        self.assertFalse(self.i["selected_PE_or_four_S_replace_complete_word"])
        self.assertFalse(self.i["ordinary_libm_trigonometric_used_in_pass_decision"])
        self.assertTrue(self.i["validated_transcendental_backend_used"])
        self.assertEqual(self.i["widest_information_cell_deg"], 30.0)
        self.assertTrue(self.i["information_headroom_closed"])
        self.assertFalse(self.i["P4_promoted_here"])

    def test_only_live_blocker_is_finite_nonlinear_remainder(self):
        self.assertGreaterEqual(self.d["P3_H_delta_consumed"], 1.0e-18)
        self.assertGreaterEqual(self.d["P3_A_delta_consumed"], 1.0e-18)
        self.assertFalse(self.d["P4_CANONICAL_PASS"])
        self.assertFalse(self.d["P5_MAY_START"])
        self.assertFalse(self.d["nonlinear_remainder_dominated_on_full_sector"])
        self.assertEqual(1, len(self.d["P4_CANONICAL_FAIL_REASONS"]))
        self.assertIn(
            "nonlinear remainder",
            self.d["P4_CANONICAL_FAIL_REASONS"][0].lower(),
        )


if __name__ == "__main__":
    unittest.main()
