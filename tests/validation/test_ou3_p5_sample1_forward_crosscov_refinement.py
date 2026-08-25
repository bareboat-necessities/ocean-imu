import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_forward_crosscov_refinement as G


class Ou3P5Sample1ForwardCrosscovRefinementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(delta_pieces=8, axial_pieces=8)

    def test_validates_fail_closed(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertGreater(d["evaluated_joint_cells"], 0)
        self.assertIn(
            d["P5_SAMPLE1_FORWARD_CROSSCOV_WITNESS_REFINEMENT"],
            ("PASS", "NOT_ESTABLISHED"),
        )
        if d["P5_SAMPLE1_FORWARD_CROSSCOV_WITNESS_REFINEMENT"] == "PASS":
            self.assertIsNone(d["first_unclosed_joint_cell"])
        else:
            self.assertIsNotNone(d["first_unclosed_joint_cell"])

    def test_forward_identity_semantics(self):
        d = self.d
        self.assertTrue(d["first_posterior_identity_P_Ht_equals_KR_used"])
        self.assertTrue(d["actual_H1_pushed_through_forward_reset_prediction_map"])
        self.assertTrue(d["mismatch_E_formed_before_covariance_multiplication"])
        self.assertTrue(d["source_reachable_force_cone_used"])
        self.assertTrue(d["sample1_S_identity_subbranch_only"])
        self.assertLessEqual(d["first_aw_tangent_gain_norm_upper"], 1.0)
        self.assertLessEqual(d["first_aw_axial_gain_abs_upper"], 1.0)

    def test_finite_diagnostics(self):
        d = self.d
        for key in (
            "max_E_theta_norm_upper",
            "max_E_bg_norm_upper",
            "max_E_aw_norm_upper",
            "max_Pj_Et_first6_norm_upper",
            "max_actual_Ctheta_norm_upper",
            "max_sample1_acc_correction_norm_upper_rad",
        ):
            self.assertTrue(math.isfinite(d[key]))
            self.assertGreaterEqual(d[key], 0.0)

    def test_no_gate_relaxation_or_promotion(self):
        d = self.d
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertFalse(d["complete_sample1_branch_refined_here"])
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])


if __name__ == "__main__":
    unittest.main()
