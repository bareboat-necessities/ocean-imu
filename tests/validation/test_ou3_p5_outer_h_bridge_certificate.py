import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_outer_h_bridge_certificate as BRIDGE


class Ou3P5OuterHBridgeCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = BRIDGE.build()

    def test_bridge_validates_and_keeps_first_s_widening(self):
        d = self.d
        self.assertEqual(BRIDGE.validate(d), [])
        self.assertEqual(d["goLive_covariance_stage"]["status"], "PASS")
        self.assertEqual(d["first_due_S_gain_stage"]["status"], "PASS")
        self.assertEqual(d["first_due_S_state_prefix_stage"]["status"], "PASS_CONDITIONAL")
        sx = d["first_due_S_exact_prefix_stage"]
        self.assertEqual(sx["status"], "PASS_WIDENED_CHART")
        self.assertFalse(sx["diagnostic_q_lt_1_is_promotion_gate"])
        self.assertGreater(sx["antipodal_margin_lower"], 0.0)

    def test_effective_vector_input_is_active_route(self):
        d = self.d
        eta = d["exact_eta_geometry_support_stage"]
        self.assertEqual(eta["status"], "PASS")
        self.assertFalse(eta["standalone_eta_penalty_is_active_P5_route"])
        self.assertFalse(eta["global_packet_count_times_Lipschitz_defect_used"])

        e = d["effective_vector_input_stage"]
        self.assertEqual(e["status"], "PASS")
        self.assertTrue(e["standalone_vector_eta_penalty_retired"])
        self.assertTrue(e["magnetometer_radial_gain_action_exact_zero"])
        self.assertTrue(e["magnetometer_effective_coordinate_nonexpansive"])
        self.assertTrue(e["accelerometer_effective_aw_norm_preserved"])
        self.assertTrue(e["gravity_quotient_uses_effective_accelerometer_input"])
        self.assertIn("K_m y_m=K_m H_m d_m", e["magnetometer_exact_state_correction_identity"])
        self.assertIn("z+E_aw e_eta", e["accelerometer_exact_state_correction_identity"])
        self.assertTrue(d["standalone_vector_eta_penalty_retired_as_P5_promotion_route"])

    def test_exact_transport_and_signed_cayley_remain_required(self):
        d = self.d
        tr = d["complete_word_transport_stage"]
        self.assertEqual(tr["status"], "PASS")
        self.assertTrue(tr["full_S_to_attitude_gain_retained"])
        self.assertTrue(tr["sequential_immediate_quaternion_resets_retained"])
        self.assertFalse(tr["standalone_vector_eta_penalty_used"])
        self.assertEqual(tr["gauged_numerical_status"], "NOT_ESTABLISHED")
        self.assertEqual(tr["quotient_numerical_status"], "NOT_ESTABLISHED")

        signed = d["signed_cayley_cell_stage"]
        self.assertEqual(signed["status"], "PASS")
        self.assertTrue(signed["signed_a_dot_c_retained"])
        self.assertFalse(signed["independent_abs_a_abs_c_denominator_used"])

    def test_corrected_outer_routes_are_retained(self):
        d = self.d
        self.assertEqual(
            d["raw_V_R_large_angle_sector_audit"]["status"],
            "DISPROVED_ON_DECLARED_SOURCE_FAMILY",
        )
        self.assertTrue(d["raw_V_R_large_angle_sector_retired_as_P5_promotion_route"])
        self.assertEqual(d["finite_angle_information_geometry"]["status"], "PASS")
        dq = d["detectable_gravity_quotient"]
        self.assertEqual(dq["status"], "PASS")
        self.assertIn("NEUTRAL_BOUNDED", dq["axial_gyro_bias_role"])
        self.assertTrue(d["yaw_only_full_bias_quotient_retired"])
        for node in d["gauged_full_heading_nodes"].values():
            self.assertEqual(node["first_due_S_exact_prefix_certificate"], "PASS_WIDENED_CHART")
            self.assertTrue(node["inside_widened_first_S_chart"])
            self.assertGreater(node["exact_pair_residual_information_vs_goLive_attitude_metric_lower"], 0.0)

    def test_remaining_obligations_are_cell_propagation(self):
        d = self.d
        self.assertEqual(
            d["gauged_full_heading_first_failure"],
            "COMPLETE_WORD_EFFECTIVE_VECTOR_INPUT_RESET_PREFIX_BUDGET_NOT_CERTIFIED",
        )
        self.assertEqual(
            d["first_failure"],
            "GRAVITY_QUOTIENT_EFFECTIVE_ACCEL_INPUT_RESET_PREFIX_BUDGET_NOT_CERTIFIED",
        )
        self.assertEqual(
            d["ungauged_timeout_route"]["current_numerical_obligation"],
            "GRAVITY_QUOTIENT_EFFECTIVE_ACCEL_INPUT_RESET_PREFIX_BUDGET_NOT_CERTIFIED",
        )
        self.assertIn("P,H,R,K,r,d_eff", d["next_full_heading_numerical_certificate"])
        self.assertIn("b_g_parallel", d["next_complete_startup_family_certificate"])
        self.assertEqual(d["P5_OUTER_H_BRIDGE_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertIsNone(d["N_H_words"])


if __name__ == "__main__":
    unittest.main()
