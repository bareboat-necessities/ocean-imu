from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_structured_full_gain_v7 as V7


class StructuredFullGainV7Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=V7.build(source_pieces=4,source_cell_index=0,p_pieces=4,tangent_pieces=4,axial_pieces=4)

    def test_validate(self):
        self.assertEqual(V7.validate(self.d),[])

    def test_analytic_block_semantics(self):
        for k in (
            "direct_first_residual_coordinate_family_retained",
            "simultaneous_Rx_inverse_gauge_is_orthogonal",
            "Ktheta_operator_norm_invariant_under_block_gauge",
            "sample1_innovation_exactly_one_plus_two_block_diagonal",
            "scalar_innovation_completed_square_positive_identity_used",
            "two_by_two_innovation_positive_determinant_identity_used",
            "sample1_nonaxial_force_included",
            "full_propagated_aw_error_norm_used",
        ):
            self.assertIs(self.d[k],True)
        for k in ("three_by_three_interval_inverse_used","fixed_pivot_inverse_used","spectral_inverse_fallback_used"):
            self.assertIs(self.d[k],False)
        self.assertGreater(self.d["minimum_scalar_innovation_lower"],0.0)
        self.assertGreater(self.d["minimum_two_by_two_determinant_lower"],0.0)

    def test_fail_closed_numeric_result(self):
        st=self.d["P5_SAMPLE1_ANALYTIC_BLOCK_GAIN_V7"]
        self.assertIn(st,("PASS","NOT_ESTABLISHED"))
        self.assertTrue(math.isfinite(self.d["max_Ktheta_operator_norm_upper"]))
        self.assertTrue(math.isfinite(self.d["max_correction_norm_upper_rad"]))
        self.assertEqual(self.d["first_unclosed_joint_cell"] is None,st=="PASS")
        if st=="PASS": self.assertLess(self.d["max_correction_norm_upper_rad"],9.0)

    def test_no_promotion(self):
        for k in (
            "source_replay_used","filter_changed","temporal_force_slew_assumed",
            "first_attitude_PSD_cross_axis_remainder_included","sample1_S_covariance_update_included",
            "sample1_S_attitude_injection_included","complete_sample1_branch_closed_here",
            "q8_word_promoted_here","whole_word_promoted_here","N_H_words_set_here",
        ):
            self.assertIs(self.d[k],False)


if __name__=="__main__": unittest.main()
