import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_go_live_covariance_stage as STAGE


class Ou3P5GoLiveCovarianceStageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = STAGE.build()

    def test_stage_is_source_bound_and_nonpromoting(self):
        d = self.d
        self.assertEqual(STAGE.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertFalse(d["global_normal_live_P3_covariance_used_at_goLive"])
        self.assertEqual(d["P5_GOLIVE_COVARIANCE_STAGE_CERTIFICATE"], "PASS")
        self.assertEqual(d["P5_FIRST_DUE_CROSS_COVARIANCE_CERTIFICATE"], "NOT_ESTABLISHED")

    def test_goLive_seed_has_exact_zero_cross_covariances(self):
        s = self.d["goLive_H_covariance_seed"]
        self.assertTrue(s["attitude_linear_cross_covariance_exact_zero"])
        self.assertTrue(s["P_awaw_reset_to_current_stationary_covariance"])
        self.assertTrue(s["S_to_attitude_gain_at_goLive_exact_zero"])
        self.assertEqual(s["pseudo_update_elapsed_s_at_goLive"], 0.0)
        for key in (
            "theta_S_cross_covariance_operator_norm_upper",
            "theta_aw_cross_covariance_operator_norm_upper",
            "bg_S_cross_covariance_operator_norm_upper",
            "aw_S_cross_covariance_operator_norm_upper",
        ):
            self.assertEqual(s[key], 0.0)
        self.assertTrue(math.isclose(s["P_SS_variance_per_axis"], 2500.0, rel_tol=0.0, abs_tol=1e-12))

    def test_first_pseudo_stage_is_finite_and_keeps_all_source_branches(self):
        s = self.d["pre_first_S_stage"]
        self.assertGreaterEqual(s["first_due_prediction_samples_upper"], 1)
        self.assertGreater(s["first_due_time_upper_s"], 0.0)
        self.assertTrue(s["source_branch_enumeration_may_not_assume_accelerometer_rejection"])
        self.assertTrue(s["accepted_accelerometer_can_create_attitude_linear_cross_via_aw"])
        self.assertTrue(s["prediction_only_preserves_zero_attitude_linear_cross_covariance"])
        self.assertFalse(s["first_due_S_to_attitude_gain_enclosed"])
        self.assertIn("PENDING_ACCEPTED_ACCEL_SOURCE_STAGE", s["pre_first_due_theta_S_cross_covariance_enclosure"])


if __name__ == "__main__":
    unittest.main()
