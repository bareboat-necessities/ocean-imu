"""Exact one-cell gyro transport certificate; no word-level promotion."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_theorem.gyro_bias_projection import certificate, source_radius  # noqa: E402


class GyroProjectionProofTests(unittest.TestCase):
    def test_source_radius_is_engineering_margin_not_truth_bound(self):
        c = certificate()
        self.assertEqual(source_radius(), F(1,2))
        self.assertEqual(F(c["estimator_radius_rad_s"])/F(c["physical_slow_bias_rad_s"]), 25)
        self.assertFalse(c["physical_amplitude_numbers_changed"])

    def test_prediction_charges_both_residual_terms(self):
        c = certificate()
        rate = sum(F(c[k]) for k in ("physical_angular_rate_rad_s", "physical_slow_bias_rad_s", "fast_measurement_amplitude_rad_s", "estimator_radius_rad_s"))
        self.assertEqual(F(c["prediction_argument_upper_rad"]), F(3,500)*rate)
        self.assertLess(F(c["prediction_argument_upper_rad"])+F(c["quaternion_polynomial_angle_defect_upper_rad"]), F(7,1000))
        self.assertGreater(F(c["margin_from_pi_lower_rad"]), 3)
        self.assertGreater(F(c["margin_from_two_pi_lower_rad"]), 6)

    def test_both_literal_transport_branches_have_positive_floor(self):
        c = certificate()
        floor = F(c["gyro_transport_singular_floor_s"])
        self.assertEqual(floor, F(1,250)*(1-F(7,1000)**2/24))
        self.assertGreater(floor, F(399999,100000000))
        self.assertGreater(F(c["small_rate_polynomial_singular_floor_s"]), floor)

    def test_changed_radius_cannot_silently_reuse_angle_certificate(self):
        with mock.patch("tools.stability.ou3_theorem.gyro_bias_projection.source_radius", return_value=F(1000)):
            with self.assertRaises(ValueError):
                certificate()

    def test_one_cell_result_does_not_promote_signed_word_or_device_timing(self):
        c = certificate()
        self.assertTrue(c["complete_turn_nominal_bias_alias_excluded"])
        for key in ("signed_temporal_Delta_gyr_closed", "force_field_collinearity_excluded", "uniform_historical_AG_action_closed", "all_positive_device_timesteps_covered", "float32_word_totality_certified"):
            self.assertFalse(c[key])

    def test_committed_certificate_matches_exact_reproduction(self):
        artifact = ROOT/"reports/results/ou3_stability/gyro-bias-projection.json"
        self.assertEqual(json.loads(artifact.read_text()), certificate())


if __name__ == "__main__":
    unittest.main()
