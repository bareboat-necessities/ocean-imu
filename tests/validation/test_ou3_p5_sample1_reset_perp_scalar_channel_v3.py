import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p5_sample1_reset_perp_scalar_channel_v3 as V3


class ResetPerpScalarChannelV3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = V3.build(p_pieces=4, d_pieces=4, axial_pieces=4)

    def test_source_state_residual_semantics(self):
        self.assertTrue(self.d["canonical_first_perpendicular_residual_zero_used"])
        self.assertTrue(self.d["first_perpendicular_aw_error_inferred_from_zero_residual"])
        self.assertTrue(self.d["reset_gauge_preserves_perpendicular_aw_x_component"])
        self.assertTrue(self.d["one_step_OU_decay_included"])
        self.assertTrue(self.d["one_step_body_rotation_aw_mixing_included"])
        self.assertFalse(self.d["global_raw_30p5_residual_multiplier_used"])
        self.assertFalse(self.d["temporal_force_slew_assumed"])

    def test_finite_fail_closed_result(self):
        self.assertGreater(self.d["evaluated_joint_cells"], 0)
        self.assertTrue(math.isfinite(self.d["max_source_correlated_perp_residual_abs_upper_mps2"]))
        self.assertTrue(math.isfinite(self.d["max_correction_norm_upper_rad"]))
        status = self.d["P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL_V3"]
        self.assertIn(status, ("PASS", "NOT_ESTABLISHED"))
        self.assertEqual(self.d["first_unclosed_joint_cell"] is None, status == "PASS")

    def test_no_overclaim(self):
        self.assertFalse(self.d["first_attitude_PSD_cross_axis_remainder_included"])
        self.assertFalse(self.d["sample1_tangent_force_perturbation_included"])
        self.assertFalse(self.d["sample1_S_attitude_correction_included"])
        self.assertFalse(self.d["complete_sample1_branch_closed_here"])
        self.assertFalse(self.d["whole_word_promoted_here"])
        self.assertFalse(self.d["N_H_words_set_here"])

    def test_validate(self):
        self.assertEqual(V3.validate(self.d), [])


if __name__ == "__main__":
    unittest.main()
