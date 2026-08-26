from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_structured_full_gain_v6 as V6


class StructuredFullGainV6Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=V6.build(source_pieces=4,source_cell_index=0,p_pieces=4,tangent_pieces=4,axial_pieces=4)

    def test_validate(self):
        self.assertEqual(V6.validate(self.d),[])

    def test_direct_residual_semantics(self):
        for k in (
            "first_residual_components_are_primary_subdivision_coordinates",
            "first_gain_interval_never_inverted_to_reconstruct_residual",
            "first_residual_norm_ball_enforced_directly",
            "first_gravity_tangent_directional_bound_used",
            "first_gravity_axial_directional_bound_used",
            "first_latent_and_bias_source_component_caps_used",
            "first_attitude_and_aw_corrections_derived_forward_from_same_residual_cell",
            "first_posterior_aw_error_residual_identity_used",
            "sample1_nonaxial_force_included",
            "complete_3x3_accelerometer_innovation_used",
            "complete_3x3_attitude_gain_used",
            "full_propagated_aw_error_norm_used",
        ):
            self.assertIs(self.d[k],True)
        self.assertLessEqual(self.d["first_axial_residual_source_cap_mps2"],self.d["first_residual_norm_upper_mps2"]+1e-12)
        self.assertGreater(self.d["evaluated_joint_cells"],0)

    def test_fail_closed_numeric_result(self):
        st=self.d["P5_SAMPLE1_DIRECT_FIRST_RESIDUAL_FULL_GAIN_V6"]
        self.assertIn(st,("PASS","NOT_ESTABLISHED"))
        self.assertTrue(math.isfinite(self.d["max_post_prediction_aw_error_norm_upper_mps2"]))
        self.assertTrue(math.isfinite(self.d["max_correction_norm_upper_rad"]))
        self.assertEqual(self.d["first_unclosed_joint_cell"] is None,st=="PASS")

    def test_no_promotion(self):
        for k in (
            "source_replay_used","filter_changed","temporal_force_slew_assumed",
            "first_attitude_PSD_cross_axis_remainder_included","sample1_S_covariance_update_included",
            "sample1_S_attitude_injection_included","complete_sample1_branch_closed_here",
            "q8_word_promoted_here","whole_word_promoted_here","N_H_words_set_here",
        ):
            self.assertIs(self.d[k],False)


if __name__=="__main__": unittest.main()
