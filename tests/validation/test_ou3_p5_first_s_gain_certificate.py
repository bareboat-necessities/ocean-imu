import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_first_s_gain_certificate as FIRST


class Ou3P5FirstSGainCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = FIRST.build()

    def test_first_s_gain_is_source_staged_and_passes(self):
        d = self.d
        self.assertEqual(FIRST.validate(d), [])
        self.assertEqual(d["P5_FIRST_DUE_S_GAIN_CERTIFICATE"], "PASS")
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertTrue(d["full_S_to_attitude_gain_retained"])
        self.assertTrue(d["accepted_and_rejected_physical_prefixes_covered"])

    def test_canonical_correlation_uses_persistent_constructor_covariance(self):
        d = self.d
        D = d["persistent_S_conditional_covariance_lambda_min_lower"]
        U = d["P_SS_lambda_max_upper_before_first_due"]
        rho = d["theta_S_canonical_correlation_upper"]
        self.assertGreaterEqual(D, 2500.0)
        self.assertGreaterEqual(U, D)
        self.assertGreater(rho, 0.0)
        self.assertLess(rho, 1.0)
        self.assertLess(rho, 1e-2)

    def test_directional_theta_bound_replaces_translation_dominated_global_bound(self):
        d = self.d
        self.assertFalse(d["global_translation_dominated_covariance_used_for_theta"])
        knew = d["K_thetaS_operator_norm_upper_first_due"]
        kold = d["old_global_P3_K_thetaS_operator_norm_upper"]
        self.assertGreater(knew, 0.0)
        self.assertLess(knew, kold)
        self.assertGreater(d["gain_widening_factor_vs_global_P3_bound_lower"], 1e4)

    def test_gain_stage_does_not_claim_state_prefix_or_finite_injection(self):
        d = self.d
        self.assertFalse(d["S_state_error_prefix_bound_supplied_here"])
        self.assertFalse(d["first_due_attitude_injection_bound_supplied_here"])
        self.assertIn("outer S-state", d["next_obligation"])
        self.assertGreater(d["timing"]["first_due_samples_upper"], 0)
        self.assertGreaterEqual(d["timing"]["aw_sync_count_before_first_due_upper"], 0)


if __name__ == "__main__":
    unittest.main()
