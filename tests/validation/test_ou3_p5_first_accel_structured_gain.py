import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_first_accel_structured_gain as G


class Ou3P5FirstAccelStructuredGainTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(source_pieces=2, alignment_pieces=16, force_magnitude_pieces=4)

    def test_structured_gain_stage_validates_and_is_source_bound(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertTrue(d["first_prefix_source_sparsity_certified"])

    def test_rank_two_gain_avoids_loose_matrix_inverse(self):
        d = self.d
        self.assertTrue(d["rotation_gauge_sets_J_aw_to_identity"])
        self.assertTrue(d["specific_force_direction_gauged_to_e3"])
        self.assertTrue(d["yaw_axis_reduced_by_axial_symmetry_to_scalar_alignment"])
        self.assertTrue(d["analytic_rank_two_gain_channels_used"])
        self.assertTrue(d["attitude_PSD_remainder_retained_by_resolvent_bound"])
        self.assertFalse(d["matrix_inverse_used_for_first_accel_gain"])
        self.assertFalse(d["loose_spectral_inverse_fallback_used"])

    def test_deployed_range_and_promotion_guards_remain_unchanged(self):
        d = self.d
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])
        self.assertGreater(d["evaluated_child_count"], 0)

    def test_status_is_fail_closed(self):
        d = self.d
        if d["P5_FIRST_ACCEL_STRUCTURED_GAIN_CERTIFICATE"] == "PASS":
            self.assertTrue(d["all_first_accelerometer_children_inside_validated_correction_range"])
            self.assertIsNone(d["first_unclosed_child"])
            self.assertEqual(
                d["next_obligation"],
                "PROPAGATE_STRUCTURED_FIRST_ACCEL_CHILDREN_THROUGH_JOSEPH_RESET_AND_LATER_PREFIXES",
            )
        else:
            self.assertIsNotNone(d["first_unclosed_child"])
            self.assertEqual(
                d["next_obligation"],
                "REFINE_ACCEPTED_ACCEL_EFFECTIVE_AW_RESIDUAL_DIRECTION_COUPLING",
            )


if __name__ == "__main__":
    unittest.main()
