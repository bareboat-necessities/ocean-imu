import math
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

    def test_bridge_consumes_exact_goLive_covariance_not_global_P3_box(self):
        d = self.d
        self.assertEqual(BRIDGE.validate(d), [])
        self.assertTrue(d["global_normal_live_P3_covariance_bound_retired_as_outer_S_prefix_gate"])
        self.assertTrue(d["old_global_P3_S_induced_attitude_bound_is_diagnostic_only"])
        g = d["goLive_covariance_stage"]
        self.assertEqual(g["status"], "PASS")
        self.assertEqual(g["theta_S_cross_covariance_operator_norm_upper"], 0.0)
        self.assertTrue(g["S_to_attitude_gain_exact_zero"])
        self.assertEqual(g["pseudo_elapsed_s"], 0.0)

    def test_bridge_consumes_first_due_gain_and_conditional_S_prefix(self):
        d = self.d
        s = d["first_due_S_gain_stage"]
        self.assertEqual(s["status"], "PASS")
        self.assertGreater(s["K_thetaS_operator_norm_upper"], 0.0)
        self.assertLess(s["theta_S_canonical_correlation_upper"], 1e-2)
        self.assertGreater(s["gain_widening_factor_vs_global_P3_bound_lower"], 1e4)
        p = d["first_due_S_state_prefix_stage"]
        self.assertEqual(p["status"], "PASS_CONDITIONAL")
        self.assertTrue(p["conditional_on_outer_node_bootstrap"])
        self.assertGreater(p["first_due_S_error_norm_upper_m_s"], 300.0)
        self.assertLess(p["first_due_S_induced_attitude_correction_norm_upper_rad"], p["group_helper_limit_rad"])
        for n in d["gauged_full_heading_nodes"].values():
            self.assertEqual(n["first_due_S_gain_certificate"], "PASS")
            self.assertEqual(n["first_due_S_state_prefix_certificate"], "PASS_CONDITIONAL")
            self.assertTrue(n["S_induced_correction_inside_group_helper"])

    def test_bridge_retires_tiny_P3_perturbation_comparison_as_outer_route(self):
        d = self.d
        self.assertTrue(d["finite_angle_perturbation_vs_tiny_P3_gap_retired_as_P5_promotion_route"])
        self.assertTrue(d["exact_large_angle_vector_dissipation_is_primary_full_heading_outer_route"])
        for n in d["gauged_full_heading_nodes"].values():
            self.assertTrue(n["exact_large_angle_sector_required"])
            self.assertFalse(n["perturbation_vs_P3_gap_is_P5_promotion_route"])
            self.assertGreater(n["perturbation_over_P3_sqrt_gap"], 1.0)
            self.assertGreaterEqual(n["handoff_group_energy_V_R_upper"], 0.0)
            self.assertLess(n["handoff_group_energy_V_R_upper"], 2.0)

    def test_bridge_uses_full_attitude_gauged_nodes_not_tilt_only_cosines(self):
        d = self.d
        h = d["heading_handoff_contract"]
        self.assertTrue(h["P1_gravity_cosines_are_tilt_only"])
        self.assertGreater(h["gauged_quality_full_cayley_norm_upper"], 0.20)
        self.assertGreater(h["gauged_timeout_full_cayley_norm_upper"], h["gauged_quality_full_cayley_norm_upper"])
        for n in d["gauged_full_heading_nodes"].values():
            self.assertTrue(n["P1_gravity_tilt_cosine_not_used_as_full_attitude_cosine"])
            self.assertTrue(n["inside_candidate_outer_cayley_bootstrap"])

    def test_ungauged_timeout_is_not_faked_as_full_heading_node(self):
        d = self.d
        q = d["ungauged_timeout_route"]
        self.assertFalse(q["full_heading_cayley_bound_available"])
        self.assertIn("YAW_QUOTIENT", q["required_route"])
        self.assertEqual(d["first_failure"], "UNGAUGED_TIMEOUT_YAW_QUOTIENT_CAPTURE_NOT_CERTIFIED")
        self.assertIn("yaw-quotient", d["next_complete_startup_family_certificate"])

    def test_next_gauged_obstruction_is_exact_large_angle_sector(self):
        d = self.d
        self.assertEqual(
            d["gauged_full_heading_first_failure"],
            "EXACT_LARGE_ANGLE_VECTOR_DISSIPATION_SECTOR_NOT_CERTIFIED",
        )
        self.assertIn("large-angle", d["next_full_heading_numerical_certificate"])
        self.assertEqual(d["P5_OUTER_H_BRIDGE_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertIsNone(d["N_H_words"])


if __name__ == "__main__":
    unittest.main()
