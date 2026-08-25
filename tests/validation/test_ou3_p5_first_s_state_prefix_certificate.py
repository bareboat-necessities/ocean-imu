import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_first_s_state_prefix_certificate as PREFIX


class Ou3P5FirstSStatePrefixCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = PREFIX.build()

    def test_conditional_prefix_is_source_bound_and_nonpromoting(self):
        d = self.d
        self.assertEqual(PREFIX.validate(d), [])
        self.assertEqual(d["P5_FIRST_S_STATE_PREFIX_CERTIFICATE"], "PASS_CONDITIONAL")
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertFalse(d["outer_node_bootstrap_supplied_here"])
        self.assertFalse(d["exact_large_angle_dissipation_supplied_here"])
        self.assertFalse(d["candidate_outer_bootstrap"]["bootstrap_proved_here"])

    def test_prefix_keeps_full_S_cross_gain_and_all_possible_physical_corrections(self):
        d = self.d
        g = d["source_correlated_final_S_gain"]
        self.assertTrue(g["full_S_cross_gain_retained"])
        self.assertGreater(g["L_K_acc_operator_norm_upper"], 0.0)
        self.assertGreater(g["L_K_mag_operator_norm_upper"], 0.0)
        timing = d["first_due_timing"]
        self.assertGreaterEqual(timing["accelerometer_packets_charged_upper"], 1)
        self.assertGreaterEqual(timing["magnetometer_packets_charged_upper"], 1)
        self.assertIn("accepted", d["accepted_rejected_branch_policy"])
        self.assertIn("rejected", d["accepted_rejected_branch_policy"])

    def test_first_due_S_bound_is_finite_and_makes_staged_S_injection_finite(self):
        d = self.d
        S = d["first_due_S_error_norm_upper_m_s"]
        inj = d["first_due_S_induced_attitude_correction_norm_upper_rad"]
        self.assertTrue(math.isfinite(S))
        self.assertGreater(S, 300.0)
        self.assertLess(S, 600.0)
        self.assertTrue(math.isfinite(inj))
        self.assertGreater(inj, 0.0)
        self.assertLess(inj, d["deployed_group_helper_correction_limit_rad"])
        self.assertTrue(d["S_induced_correction_inside_group_helper"])

    def test_candidate_chart_is_wide_but_explicitly_only_a_bootstrap(self):
        b = self.d["candidate_outer_bootstrap"]
        self.assertEqual(b["cayley_norm_upper"], 1.0)
        self.assertGreater(b["rotation_R_minus_I_norm_upper"], 0.0)
        self.assertLess(b["rotation_R_minus_I_norm_upper"], 1.0)
        self.assertEqual(b["latent_acceleration_error_norm_upper_mps2"], 10.0)
        self.assertIn("large-angle", b["required_next"])


if __name__ == "__main__":
    unittest.main()
