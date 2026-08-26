from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_structured_full_gain_v10 as V10


class StructuredFullGainV10Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=V10.build(source_pieces=4,source_cell_index=0,p_pieces=4,tangent_pieces=4,axial_pieces=4)

    def test_validate(self):
        self.assertEqual(V10.validate(self.d),[])

    def test_exact_combined_residual_decay_without_remainders(self):
        r=V10._combined_x_residual_upper(
            alpha_lo=0.9,alpha_hi=0.9,first_rot_x_upper=2.0,bias_upper=0.0,
            error_transport_rotation_norm_upper=0.0,series_rotation_mismatch_upper=0.0,
            pre_first_aw_error_norm_upper=10.0,gravity=9.8)
        self.assertGreaterEqual(r["combined_x_residual_upper_mps2"],0.2)
        self.assertLess(r["combined_x_residual_upper_mps2"],0.200000000001)
        self.assertEqual(r["rotation_mismatch_residual_upper_mps2"],0.0)

    def test_bias_difference_uses_one_plus_alpha(self):
        r=V10._combined_x_residual_upper(
            alpha_lo=0.95,alpha_hi=0.95,first_rot_x_upper=0.0,bias_upper=0.5,
            error_transport_rotation_norm_upper=0.0,series_rotation_mismatch_upper=0.0,
            pre_first_aw_error_norm_upper=0.0,gravity=9.8)
        self.assertGreaterEqual(r["bias_difference_upper_mps2"],0.975)
        self.assertLess(r["bias_difference_upper_mps2"],0.975000000001)

    def test_series_axis_mismatch_only_below_source_threshold(self):
        self.assertEqual(V10._series_vs_axis_rotation_mismatch_upper(0.02,0.03),0.0)
        self.assertGreaterEqual(V10._series_vs_axis_rotation_mismatch_upper(0.0,0.005),0.0)
        self.assertGreaterEqual(V10._series_vs_axis_rotation_mismatch_upper(0.005,0.02),0.0)

    def test_directional_result_never_exceeds_v8_parent(self):
        self.assertTrue(self.d["rows"])
        for r in self.d["rows"]:
            self.assertLessEqual(
                r["combined_directional_correction_norm_upper_rad"],
                r["parent_isotropic_correction_upper_rad"]+1e-10)

    def test_combined_residual_semantics(self):
        for k in (
            "source_generated_not_trajectory_fit","V8_positive_ratio_gain_parent_used",
            "analytic_one_plus_two_block_structure_retained","first_perpendicular_residual_exact_zero_in_ideal_SO2_gauge",
            "first_perpendicular_aw_estimator_correction_exact_zero_in_ideal_block",
            "latent_term_carried_as_rotated_E_eaw_component",
            "first_force_and_post_update_aw_correction_cancel_in_sample1_x_residual",
            "nominal_body_rotation_removed_by_V7_simultaneous_orthogonal_gauge",
            "same_one_step_error_transport_bound_reused_for_sample1",
            "deployed_series_vs_axis_gauge_mismatch_included","above_series_threshold_axis_gauge_mismatch_exact_zero",
            "orthogonal_two_block_residual_ball_optimization_used",
        ):
            self.assertIs(self.d[k],True)
        self.assertLessEqual(
            self.d["max_combined_source_x_residual_upper_mps2"],
            self.d["max_full_residual_norm_upper_mps2"]+1e-12)

    def test_fail_closed_numeric_result(self):
        st=self.d["P5_SAMPLE1_COMBINED_PERPENDICULAR_RESIDUAL_V10"]
        self.assertIn(st,("PASS","NOT_ESTABLISHED"))
        self.assertTrue(math.isfinite(self.d["max_combined_directional_correction_norm_upper_rad"]))
        self.assertEqual(self.d["first_unclosed_joint_cell"] is None,st=="PASS")
        self.assertEqual(self.d["unclosed_joint_cells"]==0,st=="PASS")
        if st=="PASS":
            self.assertLess(self.d["max_combined_directional_correction_norm_upper_rad"],9.0)

    def test_no_promotion(self):
        for k in (
            "source_replay_used","filter_changed","temporal_force_slew_assumed","raw_aw_component_used_for_x_latent_residual",
            "large_scalar_gain_multiplied_by_full_residual_norm","first_attitude_PSD_cross_axis_remainder_included",
            "sample1_S_covariance_update_included","sample1_S_attitude_injection_included",
            "complete_sample1_branch_closed_here","q8_word_promoted_here","whole_word_promoted_here","N_H_words_set_here",
        ):
            self.assertIs(self.d[k],False)


if __name__=="__main__": unittest.main()
