import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_first_accel_rotation_gauge_v2 as G


class Ou3P5FirstAccelRotationGaugeV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(
            source_pieces=2,
            yaw_axis_face_pieces=4,
            force_magnitude_pieces=4,
        )

    def test_v2_hulls_only_equivalent_axis_interval_representations(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertTrue(d["axis_isotropic_source_intervals_hulled_across_equivalent_axes"])
        self.assertFalse(d["bit_identical_axis_interval_endpoints_required"])
        self.assertTrue(d["cross_axis_covariance_still_required_exact_zero"])
        self.assertTrue(d["theta_aw_covariance_still_required_exact_zero"])

    def test_source_and_deployed_range_contracts_remain_unchanged(self):
        d = self.d
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertTrue(d["rotation_gauge_sets_J_aw_to_identity"])
        self.assertTrue(d["rotation_gauge_sets_specific_force_direction_to_e3"])

    def test_innovation_inverse_must_close_without_spectral_fallback(self):
        d = self.d
        self.assertGreater(d["evaluated_child_count"], 0)
        self.assertGreater(d["fixed_pivot_inverse_count"], 0)
        self.assertEqual(d["spectral_fallback_inverse_count"], 0)

    def test_stage_remains_fail_closed_and_never_sets_NH(self):
        d = self.d
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])
        if d["P5_FIRST_ACCEL_ROTATION_GAUGED_CERTIFICATE"] == "PASS":
            self.assertTrue(d["all_first_accelerometer_children_inside_validated_correction_range"])
            self.assertIsNone(d["first_unclosed_child"])
            self.assertEqual(
                d["next_obligation"],
                "PROPAGATE_ROTATION_GAUGED_CHILDREN_THROUGH_ACCEL_JOSEPH_RESET_AND_LATER_PREFIXES",
            )
        else:
            self.assertIsNotNone(d["first_unclosed_child"])
            self.assertEqual(
                d["next_obligation"],
                "REFINE_FIRST_ACCEL_ATTITUDE_COVARIANCE_AND_EFFECTIVE_AW_DIRECTION_COUPLING",
            )


if __name__ == "__main__":
    unittest.main()
