import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p5_sample1_reset_perp_scalar_channel_v2 as V2


class ResetPerpScalarChannelV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = V2.build(p_pieces=4, d_pieces=4, axial_pieces=4)

    def test_positive_dependency_preserving_semantics(self):
        for key in (
            "first_posterior_determinant_identity_used",
            "reset_determinant_identity_used",
            "one_step_process_determinant_identity_used",
            "positive_innovation_identity_used",
            "direct_scalar_ratio_maximization_used",
        ):
            self.assertTrue(self.d[key])
        self.assertGreater(self.d["minimum_A_lower"], 0.0)
        self.assertGreaterEqual(self.d["minimum_predicted_determinant_lower"], 0.0)

    def test_finite_result(self):
        self.assertGreater(self.d["evaluated_joint_cells"], 0)
        self.assertTrue(math.isfinite(self.d["max_Ktheta_abs_upper"]))
        self.assertTrue(math.isfinite(self.d["max_correction_norm_upper_rad"]))

    def test_fail_closed_and_no_promotion(self):
        self.assertFalse(self.d["source_replay_used"])
        self.assertFalse(self.d["filter_changed"])
        self.assertFalse(self.d["source_attitude_remainder_cross_terms_included"])
        self.assertFalse(self.d["sample1_body_rate_rotation_perturbation_included"])
        self.assertFalse(self.d["sample1_tangent_force_perturbation_included"])
        self.assertFalse(self.d["sample1_S_attitude_correction_included"])
        self.assertFalse(self.d["complete_sample1_branch_closed_here"])
        self.assertFalse(self.d["whole_word_promoted_here"])
        self.assertFalse(self.d["N_H_words_set_here"])
        status = self.d["P5_SAMPLE1_RESET_PERP_SCALAR_CHANNEL_V2"]
        witness = self.d["first_unclosed_joint_cell"]
        self.assertIn(status, ("PASS", "NOT_ESTABLISHED"))
        self.assertEqual(witness is None, status == "PASS")

    def test_validate(self):
        self.assertEqual(V2.validate(self.d), [])


if __name__ == "__main__":
    unittest.main()
