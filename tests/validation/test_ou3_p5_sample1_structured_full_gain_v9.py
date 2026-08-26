from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_structured_full_gain_v9 as V9


class StructuredFullGainV9Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=V9.build(source_pieces=4,source_cell_index=0,p_pieces=4,tangent_pieces=4,axial_pieces=4)

    def test_two_block_allocation_when_perp_dominates(self):
        v=V9._two_block_correction_upper(2.0,1.0,3.0,5.0)
        self.assertGreaterEqual(v,math.sqrt(52.0))
        self.assertLess(v,math.sqrt(52.0)+1e-12)

    def test_two_block_allocation_when_parallel_dominates(self):
        v=V9._two_block_correction_upper(1.0,2.0,3.0,5.0)
        self.assertGreaterEqual(v,10.0)
        self.assertLess(v,10.0+1e-12)

    def test_validate(self):
        self.assertEqual(V9.validate(self.d),[])

    def test_directional_semantics(self):
        for k in (
            "V8_positive_ratio_gain_parent_used",
            "analytic_one_plus_two_block_structure_retained",
            "canonical_first_perpendicular_residual_zero_used",
            "first_exact_gravity_chord_residual_used_in_perpendicular_component",
            "reset_gauge_preserves_perpendicular_aw_component",
            "one_step_body_rotation_aw_mixing_included",
            "orthogonal_two_block_residual_ball_optimization_used",
        ):
            self.assertIs(self.d[k],True)
        self.assertIs(self.d["large_scalar_gain_multiplied_by_full_residual_norm"],False)
        self.assertIs(self.d["same_full_residual_norm_double_counted_across_blocks"],False)
        self.assertGreater(self.d["evaluated_joint_cells"],0)
        self.assertLessEqual(
            self.d["max_directional_correction_norm_upper_rad"],
            self.d["max_parent_isotropic_correction_upper_rad"]+1e-10)
        for row in self.d["rows"]:
            self.assertLessEqual(
                row["directional_correction_norm_upper_rad"],
                row["parent_isotropic_correction_upper_rad"]+1e-10)

    def test_fail_closed_numeric_result(self):
        st=self.d["P5_SAMPLE1_DIRECTIONAL_RESIDUAL_BLOCK_GAIN_V9"]
        self.assertIn(st,("PASS","NOT_ESTABLISHED"))
        self.assertTrue(math.isfinite(self.d["max_directional_correction_norm_upper_rad"]))
        if st=="PASS":
            self.assertEqual(self.d["unclosed_joint_cells"],0)
            self.assertIsNone(self.d["first_unclosed_joint_cell"])
            self.assertLess(self.d["max_directional_correction_norm_upper_rad"],9.0)
        else:
            self.assertGreater(self.d["unclosed_joint_cells"],0)
            self.assertIsNotNone(self.d["first_unclosed_joint_cell"])

    def test_no_promotion(self):
        for k in (
            "source_replay_used","filter_changed","temporal_force_slew_assumed",
            "first_attitude_PSD_cross_axis_remainder_included","sample1_S_covariance_update_included",
            "sample1_S_attitude_injection_included","complete_sample1_branch_closed_here",
            "q8_word_promoted_here","whole_word_promoted_here","N_H_words_set_here",
        ):
            self.assertIs(self.d[k],False)


if __name__=="__main__": unittest.main()
