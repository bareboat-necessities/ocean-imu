import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_yaw_quotient_word_certificate as QUOT


class Ou3P5YawQuotientWordCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = QUOT.build()

    def test_exact_gravity_only_zero_dynamics_witness_is_source_valid(self):
        d = self.d
        self.assertEqual(QUOT.validate(d), [])
        w = d["witness"]
        self.assertTrue(w["zero_dynamics_source_word_valid"])
        self.assertEqual(w["body_rate_rad_s"], 0.0)
        self.assertEqual(w["a_w_error_norm_mps2"], 0.0)
        self.assertEqual(w["post_word_quotient_attitude_error"], 0.0)
        self.assertGreater(w["created_yaw_before_quotient_interval_rad"][0], 0.0)
        self.assertFalse(w["magnetometer_correction_available"])

    def test_yaw_only_quotient_does_not_contract_axial_gyro_bias(self):
        d = self.d
        r = d["witness"]["parallel_bias_physical_contraction_ratio_interval"]
        self.assertLessEqual(r[0], 1.0)
        self.assertGreaterEqual(r[1], 1.0)
        self.assertFalse(d["strict_lambda_less_than_one_possible_on_yaw_only_quotient"])
        self.assertEqual(d["P5_YAW_ONLY_QUOTIENT_OBSTRUCTION_IDENTIFIED"], "PASS")
        self.assertEqual(d["P5_UNGAUGED_TIMEOUT_QUOTIENT_WORD_CERTIFICATE"], "NOT_ESTABLISHED")

    def test_unobserved_axial_bias_covariance_grows(self):
        c = self.d["covariance_obstruction"]
        self.assertGreater(c["one_word_unobserved_parallel_bias_covariance_growth_interval"][0], 0.0)
        self.assertGreater(c["repeated_unobserved_covariance_growth_lower"], 0.0)
        self.assertFalse(c["uniform_compact_information_metric_possible_if_parallel_bias_retained"])

    def test_corrected_quotient_must_not_penalize_neutral_direction(self):
        d = self.d
        q = d["required_quotient_correction"]
        self.assertTrue(q["remove_horizontal_attitude_gauge"])
        self.assertTrue(q["do_not_require_strict_contraction_of_instantaneous_gravity_parallel_gyro_bias"])
        self.assertIn("bounded input", q["parallel_gyro_bias_role"])
        self.assertTrue(d["paper_theorem_requires_revision"])


if __name__ == "__main__":
    unittest.main()
