import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_first_accel_exact_source_v2 as G


class Ou3P5FirstAccelExactSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(source_pieces=2)

    def test_source_exact_first_accel_stage_validates(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertTrue(d["first_linear_prediction_homogeneous_zero_mean_certified"])
        self.assertTrue(d["first_S_pseudo_zero_residual_mean_identity_certified"])

    def test_first_force_and_covariance_axis_are_source_reachable_not_generic(self):
        d = self.d
        self.assertTrue(d["first_accel_aw_mean_exact_zero_before_measurement"])
        self.assertTrue(d["first_due_S_mean_correction_exact_zero"])
        self.assertTrue(d["first_accel_specific_force_magnitude_exact_gravity"])
        self.assertTrue(d["first_accel_yaw_covariance_axis_aligned_with_force_axis"])
        self.assertTrue(d["yaw_alignment_x_equals_zero"])
        self.assertGreater(d["first_accel_specific_force_magnitude_mps2"], 9.0)
        self.assertLess(d["first_accel_specific_force_magnitude_mps2"], 10.5)

    def test_physical_gravity_error_only_uses_tangent_cayley_component(self):
        d = self.d
        self.assertTrue(d["exact_gravity_cayley_tangent_identity_used"])
        self.assertTrue(d["exact_rotational_residual_tangent_bound_used"])
        self.assertGreater(d["post_prediction_full_cayley_norm_upper"], 0.0)
        self.assertGreaterEqual(d["post_prediction_cayley_tangent_norm_upper"], 0.0)
        self.assertLessEqual(
            d["post_prediction_cayley_tangent_norm_upper"],
            d["post_prediction_full_cayley_norm_upper"],
        )

    def test_latent_rotation_cross_is_combined_exactly(self):
        d = self.d
        self.assertTrue(d["latent_linear_plus_rotation_cross_combined_before_norm"])
        self.assertFalse(d["independent_latent_rotation_cross_norm_added"])
        self.assertIn("||R^T e||=||e||", d["latent_combined_norm_identity"])

    def test_deployed_range_and_whole_word_guards_are_unchanged(self):
        d = self.d
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])
        self.assertGreater(d["evaluated_source_phase_cells"], 0)
        if d["P5_FIRST_ACCEL_EXACT_SOURCE_CERTIFICATE"] == "PASS":
            self.assertTrue(d["all_first_accelerometer_source_cells_inside_validated_correction_range"])
            self.assertIsNone(d["first_unclosed_child"])
        else:
            self.assertIsNotNone(d["first_unclosed_child"])


if __name__ == "__main__":
    unittest.main()
