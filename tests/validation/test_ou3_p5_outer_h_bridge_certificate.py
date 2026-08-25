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

    def test_bridge_keeps_source_staged_first_S_work(self):
        d = self.d
        self.assertEqual(BRIDGE.validate(d), [])
        self.assertTrue(d["global_normal_live_P3_covariance_bound_retired_as_outer_S_prefix_gate"])
        self.assertTrue(d["finite_angle_perturbation_vs_tiny_P3_gap_retired_as_P5_promotion_route"])
        self.assertEqual(d["goLive_covariance_stage"]["status"], "PASS")
        self.assertTrue(d["goLive_covariance_stage"]["S_to_attitude_gain_exact_zero"])
        self.assertEqual(d["first_due_S_gain_stage"]["status"], "PASS")
        p = d["first_due_S_state_prefix_stage"]
        self.assertEqual(p["status"], "PASS_CONDITIONAL")
        self.assertLess(p["first_due_S_induced_attitude_correction_norm_upper_rad"], p["group_helper_limit_rad"])
        x = d["first_due_S_exact_prefix_stage"]
        self.assertEqual(x["status"], "PASS_WIDENED_CHART")
        self.assertFalse(x["diagnostic_q_lt_1_is_promotion_gate"])
        self.assertGreater(x["antipodal_margin_lower"], 0.0)
        self.assertGreater(x["vector_information_vs_goLive_metric_lower"], 0.0)

    def test_bridge_consumes_validated_raw_VR_counterexample(self):
        d = self.d
        a = d["raw_V_R_large_angle_sector_audit"]
        self.assertEqual(a["status"], "DISPROVED_ON_DECLARED_SOURCE_FAMILY")
        self.assertTrue(a["beta_cannot_repair_xi_zero_counterexample"])
        self.assertLess(a["counterexample"]["D_R_deployed_interval"][1], 0.0)
        self.assertTrue(d["raw_V_R_large_angle_sector_retired_as_P5_promotion_route"])
        self.assertTrue(d["source_shaped_Cayley_information_is_primary_full_heading_outer_route"])
        for n in d["gauged_full_heading_nodes"].values():
            self.assertTrue(n["raw_V_R_sector_witness_inside_node"])
            self.assertFalse(n["raw_V_R_sector_is_P5_promotion_route"])
            self.assertTrue(n["source_shaped_Cayley_information_outer_sector_required"])

    def test_finite_angle_information_and_exact_transport_prerequisites_are_closed(self):
        d = self.d
        fi = d["finite_angle_information_geometry"]
        self.assertEqual(fi["status"], "PASS")
        self.assertEqual(fi["complete_word_sector"], "NOT_ESTABLISHED")
        self.assertIn("eta^T R^-1 eta", fi["exact_joseph_tangent_information_identity"])
        tr = d["complete_word_transport_stage"]
        self.assertEqual(tr["status"], "PASS")
        self.assertTrue(tr["full_S_to_attitude_gain_retained"])
        self.assertTrue(tr["sequential_immediate_quaternion_resets_retained"])
        self.assertEqual(tr["gauged_numerical_status"], "NOT_ESTABLISHED")
        eta = d["exact_eta_subdivision_stage"]
        self.assertEqual(eta["status"], "PASS")
        self.assertGreaterEqual(eta["annular_subdivision_cell_count"], 1)
        self.assertFalse(eta["global_packet_count_times_Lipschitz_defect_used"])
        signed = d["signed_cayley_cell_stage"]
        self.assertEqual(signed["status"], "PASS")
        self.assertTrue(signed["signed_a_dot_c_retained"])
        self.assertFalse(signed["independent_abs_a_abs_c_denominator_used"])
        for n in d["gauged_full_heading_nodes"].values():
            self.assertEqual(n["finite_angle_information_geometry_status"], "PASS")
            self.assertGreater(n["exact_cayley_residual_factor_lower"], 0.0)
            self.assertGreater(n["exact_pair_residual_information_per_cayley_norm_sq_lower"], 0.0)
            self.assertGreater(n["exact_pair_residual_information_vs_goLive_attitude_metric_lower"], 0.0)
            self.assertTrue(n["source_correlated_Joseph_information_identity_retained"])
            self.assertEqual(n["first_due_S_exact_prefix_certificate"], "PASS_WIDENED_CHART")
            self.assertTrue(n["inside_widened_first_S_chart"])

    def test_bridge_uses_full_attitude_gauged_nodes_not_tilt_only_cosines(self):
        d = self.d
        h = d["heading_handoff_contract"]
        self.assertTrue(h["P1_gravity_cosines_are_tilt_only"])
        self.assertGreater(h["gauged_quality_full_cayley_norm_upper"], 0.20)
        self.assertGreater(h["gauged_timeout_full_cayley_norm_upper"], h["gauged_quality_full_cayley_norm_upper"])
        for n in d["gauged_full_heading_nodes"].values():
            self.assertTrue(n["P1_gravity_tilt_cosine_not_used_as_full_attitude_cosine"])
            self.assertTrue(n["inside_candidate_outer_cayley_bootstrap"])
            self.assertTrue(n["S_induced_correction_inside_group_helper"])

    def test_detectable_gravity_quotient_replaces_false_yaw_only_route(self):
        d = self.d
        q = d["yaw_only_quotient_audit"]
        self.assertEqual(q["obstruction_identified"], "PASS")
        self.assertEqual(q["status"], "NOT_ESTABLISHED")
        dq = d["detectable_gravity_quotient"]
        self.assertEqual(dq["status"], "PASS")
        self.assertEqual(dq["complete_word_status"], "NOT_ESTABLISHED")
        self.assertGreater(dq["reduced_attitude_bias_information"]["alpha_4_quotient_information_lower"], 0.0)
        self.assertTrue(dq["translation_word"]["source_complete"])
        self.assertIn("NEUTRAL_BOUNDED", dq["axial_gyro_bias_role"])
        u = d["ungauged_timeout_route"]
        self.assertTrue(u["yaw_only_quotient_disproved"])
        self.assertEqual(u["reduced_detectability_certificate"], "PASS")
        self.assertFalse(u["full_heading_cayley_bound_available"])
        self.assertIn("AXIAL_GYRO_BIAS", u["required_route"])
        self.assertEqual(u["current_numerical_obligation"], "GRAVITY_QUOTIENT_EXACT_ETA_RESET_PREFIX_BUDGET_NOT_CERTIFIED")

    def test_remaining_obligations_are_later_numerical_subdivision_only(self):
        d = self.d
        self.assertEqual(
            d["gauged_full_heading_first_failure"],
            "COMPLETE_WORD_ETA_RESET_INFORMATION_BUDGET_NOT_CERTIFIED",
        )
        self.assertEqual(
            d["first_failure"],
            "GRAVITY_QUOTIENT_EXACT_ETA_RESET_PREFIX_BUDGET_NOT_CERTIFIED",
        )
        self.assertIn("source-correlated covariance/gain/residual K*r cells", d["next_full_heading_numerical_certificate"])
        self.assertIn("b_g_parallel", d["next_complete_startup_family_certificate"])
        self.assertEqual(d["P5_OUTER_H_BRIDGE_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertIsNone(d["N_H_words"])


if __name__ == "__main__":
    unittest.main()
