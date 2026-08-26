import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_first_accel_rotation_gauge_v3 as G


class Ou3P5FirstAccelRotationGaugeV3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(
            source_pieces=2,
            yaw_axis_face_pieces=4,
            force_magnitude_pieces=4,
        )

    def test_structural_zeros_are_source_certified_not_thresholded(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertTrue(d["first_prefix_source_sparsity_certified_before_interval_product"])
        self.assertGreater(d["source_sparsity_cells_checked"], 0)
        self.assertFalse(d["structural_zero_canonicalization_uses_numeric_threshold"])
        self.assertFalse(d["cross_axis_interval_dust_treated_as_physical_covariance"])
        self.assertTrue(d["axis_equivalent_same_axis_intervals_hulled"])
        self.assertGreaterEqual(d["arithmetic_zero_dust_abs_upper_max"], 0.0)

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

    def test_stage_never_promotes_complete_word(self):
        d = self.d
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])
        if d["P5_FIRST_ACCEL_ROTATION_GAUGED_CERTIFICATE"] == "PASS":
            self.assertTrue(d["all_first_accelerometer_children_inside_validated_correction_range"])
            self.assertIsNone(d["first_unclosed_child"])
        else:
            self.assertIsNotNone(d["first_unclosed_child"])


if __name__ == "__main__":
    unittest.main()
