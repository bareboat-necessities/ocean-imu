from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_structured_full_gain_v5 as V5


class StructuredFullGainV5Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = V5.build(source_pieces=4, source_cell_index=0,
                         p_pieces=4, d_pieces=4, axial_pieces=4)

    def test_validate(self):
        self.assertEqual(V5.validate(self.d), [])

    def test_residual_coupling_semantics(self):
        for k in (
            "first_residual_to_tangent_correction_identity_used",
            "first_residual_to_axial_aw_correction_identity_used",
            "first_residual_norm_ball_enforced_on_joint_d_axial_cells",
            "first_posterior_aw_error_residual_identity_used",
            "posterior_aw_bound_uses_min_of_residual_identity_and_triangle",
            "sample1_nonaxial_force_included",
            "complete_3x3_accelerometer_innovation_used",
            "complete_3x3_attitude_gain_used",
            "full_propagated_aw_error_norm_used",
        ):
            self.assertIs(self.d[k], True)
        self.assertIs(self.d["independent_prior_plus_correction_aw_bound_used_as_only_bound"], False)
        self.assertGreater(self.d["evaluated_joint_cells"], 0)
        self.assertGreaterEqual(self.d["candidate_joint_cells_before_residual_ball"], self.d["evaluated_joint_cells"])

    def test_fail_closed_numeric_result(self):
        self.assertIn(self.d["P5_SAMPLE1_FIRST_RESIDUAL_COUPLED_FULL_GAIN_V5"],
                      ("PASS", "NOT_ESTABLISHED"))
        self.assertTrue(math.isfinite(self.d["max_post_prediction_aw_error_norm_upper_mps2"]))
        self.assertTrue(math.isfinite(self.d["max_correction_norm_upper_rad"]))
        if self.d["P5_SAMPLE1_FIRST_RESIDUAL_COUPLED_FULL_GAIN_V5"] == "PASS":
            self.assertIsNone(self.d["first_unclosed_joint_cell"])
            self.assertLess(self.d["max_correction_norm_upper_rad"], 9.0)
        else:
            self.assertIsNotNone(self.d["first_unclosed_joint_cell"])

    def test_no_promotion(self):
        for k in (
            "first_attitude_PSD_cross_axis_remainder_included",
            "sample1_S_covariance_update_included",
            "sample1_S_attitude_injection_included",
            "complete_sample1_branch_closed_here",
            "q8_word_promoted_here",
            "whole_word_promoted_here",
            "N_H_words_set_here",
            "filter_changed",
            "source_replay_used",
        ):
            self.assertIs(self.d[k], False)


if __name__ == "__main__":
    unittest.main()
