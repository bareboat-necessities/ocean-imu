from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_two_word_packet_null_lift as LIFT
from ou3_proof_module_state import preserve_module_bindings


class TwoWordPacketNullLiftTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with preserve_module_bindings():
            cls.d = LIFT.build()

    def test_structural_null_lift_is_source_uniform_and_fail_closed(self):
        self.assertEqual(LIFT.validate(self.d), [])
        self.assertTrue(self.d["two_word_packet_null_is_structurally_observed"])
        self.assertTrue(self.d["shared_H18_differential_operations_used"])
        self.assertGreater(self.d["packet_null_aw_per_theta_norm_lower_mps2_per_rad"], 0.0)
        self.assertGreater(self.d["aw_survival_factor_to_following_word_lower"], 0.0)
        self.assertGreater(self.d["following_word_aw_per_theta_norm_lower"], 0.0)
        self.assertGreater(self.d["following_word_four_S_information_gramian_lambda_min_lower"], 0.0)
        self.assertGreater(self.d["packet_null_following_word_raw_information_per_theta2_lower"], 0.0)
        self.assertEqual(
            self.d["packet_null_following_word_raw_information_bound_formula"],
            "lambda_min(G_S)*(a_w_next/theta)^2",
        )
        self.assertFalse(self.d["raw_coordinate_product_is_final_P4_metric_margin"])
        self.assertFalse(self.d["full_H18_metric_directional_credit_established_here"])
        self.assertFalse(self.d["P4_USABLE_CERTIFICATE_PROMOTED"])

    def test_PE_lower_bound_is_not_weakened(self):
        # Declared shipping theorem values currently imply f_min*sin(phi_min)=0.5.
        self.assertGreaterEqual(
            self.d["packet_null_aw_per_theta_norm_lower_mps2_per_rad"],
            0.5 - 1e-15,
        )
        self.assertEqual(self.d["following_word_four_S_firing_count"], 4)

    def test_raw_information_credit_is_derived_from_existing_four_S_gramian(self):
        expected = (
            self.d["following_word_four_S_information_gramian_lambda_min_lower"]
            * self.d["following_word_aw_per_theta_norm_lower"] ** 2
        )
        actual = self.d["packet_null_following_word_raw_information_per_theta2_lower"]
        self.assertGreater(actual, 0.0)
        self.assertLessEqual(actual, expected)
        self.assertGreater(actual, expected * (1.0 - 1e-12))

    def test_next_step_is_metric_pullback_not_new_domain_assumption(self):
        text = self.d["next_obligation"]
        self.assertIn("shared H18 interval-AD", text)
        self.assertIn("pull back", text)
        self.assertIn("before scalarization", text)
        self.assertIn("actual source-correlated covariance/information metric", text)


if __name__ == "__main__":
    unittest.main()
