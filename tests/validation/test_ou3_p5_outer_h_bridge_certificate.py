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

    def test_bridge_consumes_first_due_source_staged_S_gain(self):
        d = self.d
        s = d["first_due_S_gain_stage"]
        self.assertEqual(s["status"], "PASS")
        self.assertGreater(s["K_thetaS_operator_norm_upper"], 0.0)
        self.assertLess(s["theta_S_canonical_correlation_upper"], 1e-2)
        self.assertGreater(s["gain_widening_factor_vs_global_P3_bound_lower"], 1e4)
        for name in ("normal", "timeout"):
            n = d["nodes"][name]
            self.assertEqual(n["first_due_S_gain_certificate"], "PASS")
            self.assertFalse(n["handoff_radius_is_certified_prefix_bound"])
            self.assertGreater(n["first_due_injection_if_S_never_exceeded_handoff_radius_diagnostic"], 0.0)

    def test_bridge_retires_tiny_P3_perturbation_comparison_as_outer_route(self):
        d = self.d
        self.assertTrue(d["finite_angle_perturbation_vs_tiny_P3_gap_retired_as_P5_promotion_route"])
        self.assertTrue(d["exact_large_angle_vector_dissipation_is_primary_outer_route"])
        for name in ("normal", "timeout"):
            n = d["nodes"][name]
            self.assertTrue(n["exact_large_angle_sector_required"])
            self.assertFalse(n["perturbation_vs_P3_gap_is_P5_promotion_route"])
            self.assertGreater(n["perturbation_over_P3_sqrt_gap"], 1.0)
            self.assertGreaterEqual(n["handoff_group_energy_V_R_upper"], 0.0)
            self.assertLess(n["handoff_group_energy_V_R_upper"], 2.0)

    def test_bridge_stops_at_outer_S_state_prefix_not_covariance_gain(self):
        d = self.d
        self.assertEqual(d["first_failure"], "OUTER_S_STATE_PREFIX_NOT_CERTIFIED_TO_FIRST_DUE_S")
        self.assertEqual(d["P5_OUTER_H_BRIDGE_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertIsNone(d["N_H_words"])
        self.assertIn("S-state prefix", d["first_required_numerical_certificate"])
        self.assertIn("large-angle", d["second_required_numerical_certificate"])


if __name__ == "__main__":
    unittest.main()
