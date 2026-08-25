import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_first_accel_rotation_gauge as G


class Ou3P5FirstAccelRotationGaugeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(
            source_pieces=2,
            yaw_axis_face_pieces=4,
            force_magnitude_pieces=4,
        )

    def test_stage_is_source_bound_and_does_not_change_filter_or_range(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])

    def test_rotation_gauge_consumes_exact_first_prefix_structure(self):
        d = self.d
        self.assertTrue(d["rotation_gauge_sets_J_aw_to_identity"])
        self.assertTrue(d["rotation_gauge_sets_specific_force_direction_to_e3"])
        self.assertTrue(d["gauge_requires_first_prefix_theta_aw_cross_zero"])
        self.assertTrue(d["gauge_requires_first_prefix_aw_axis_isotropy"])
        self.assertTrue(d["first_prefix_source_structure_checked"])
        self.assertTrue(d["attitude_seed_rank_one_yaw_axis_retained"])
        self.assertTrue(d["yaw_axis_cube_face_cover_complete"])
        self.assertTrue(d["pseudo_phase_coupled_to_tau_before_branching"])

    def test_loose_innovation_inverse_fallback_is_retired_for_this_stage(self):
        d = self.d
        self.assertGreater(d["evaluated_child_count"], 0)
        self.assertGreater(d["fixed_pivot_inverse_count"], 0)
        self.assertEqual(d["spectral_fallback_inverse_count"], 0)

    def test_first_stage_never_promotes_complete_word(self):
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
