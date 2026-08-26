from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
from ou3_interval import Interval
import ou3_p5_sample1_structured_full_gain_v8 as V8


class StructuredFullGainV8Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=V8.build(source_pieces=4,source_cell_index=0,p_pieces=4,tangent_pieces=4,axial_pieces=4)

    def test_validate(self):
        self.assertEqual(V8.validate(self.d),[])

    def test_positive_ratio_stationary_line(self):
        # f=(x+z)/(x+z+1)^2 has its maximum 1/4 on x+z=1.
        val,detail=V8._linear_over_square_rect_max(
            1.0,1.0,Interval.outward_bounds(0.0,2.0),Interval.outward_bounds(0.0,2.0),1.0)
        self.assertGreaterEqual(val,0.25)
        self.assertLess(val,0.250000000001)
        self.assertAlmostEqual(sum(detail["maximizer"]),1.0,places=12)

    def test_positive_ratio_corner_maximum(self):
        # f=x/(x+z+1)^2 is decreasing for x>=2,z>=0 here; max=2/9.
        val,detail=V8._linear_over_square_rect_max(
            1.0,0.0,Interval.outward_bounds(2.0,3.0),Interval.outward_bounds(0.0,1.0),1.0)
        self.assertGreaterEqual(val,2.0/9.0)
        self.assertLess(val,2.0/9.0+1e-12)
        self.assertEqual(detail["maximizer"],[2.0,0.0])

    def test_ratio_dependency_semantics(self):
        for k in (
            "direct_first_residual_coordinate_family_retained",
            "analytic_one_plus_two_block_structure_retained",
            "block_numerator_denominator_dependency_preserved_by_positive_ratio_maximization",
            "rectangle_boundary_stationary_maximization_used",
            "sample1_nonaxial_force_included",
            "full_propagated_aw_error_norm_used",
        ):
            self.assertIs(self.d[k],True)
        self.assertIs(self.d["three_by_three_interval_inverse_used"],False)
        self.assertIs(self.d["spectral_inverse_fallback_used"],False)
        self.assertGreater(self.d["minimum_scalar_innovation_lower"],0.0)
        self.assertGreater(self.d["minimum_two_by_two_determinant_lower"],0.0)

    def test_fail_closed_numeric_result(self):
        st=self.d["P5_SAMPLE1_POSITIVE_RATIO_BLOCK_GAIN_V8"]
        self.assertIn(st,("PASS","NOT_ESTABLISHED"))
        self.assertTrue(math.isfinite(self.d["max_Ktheta_operator_norm_upper"]))
        self.assertTrue(math.isfinite(self.d["max_correction_norm_upper_rad"]))
        self.assertEqual(self.d["first_unclosed_joint_cell"] is None,st=="PASS")
        if st=="PASS": self.assertLess(self.d["max_correction_norm_upper_rad"],9.0)

    def test_no_promotion(self):
        for k in (
            "source_replay_used","filter_changed","three_by_three_interval_inverse_used","spectral_inverse_fallback_used",
            "first_attitude_PSD_cross_axis_remainder_included","sample1_S_covariance_update_included",
            "sample1_S_attitude_injection_included","complete_sample1_branch_closed_here",
            "q8_word_promoted_here","whole_word_promoted_here","N_H_words_set_here",
        ):
            self.assertIs(self.d[k],False)


if __name__=="__main__": unittest.main()
