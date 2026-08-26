from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_structured_full_gain_v11 as V11


class StructuredFullGainV11Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # This deliberately coarse fixture is a semantic/fail-closed test, not
        # the authoritative V10/V11 numerical certificate. V10 requires the
        # focused 24^3 grid to close; the workflow emitter below supplies that.
        cls.d=V11.build(source_pieces=4,source_cell_index=0,p_pieces=4,tangent_pieces=4,axial_pieces=4)

    def test_validate_fail_closed_coarse_prerequisite(self):
        failures=V11.validate(self.d)
        self.assertEqual(failures,["V10 prerequisite did not pass"])
        self.assertEqual(self.d["P5_SAMPLE1_PSD_S_PERTURBATION_V11"],"NOT_ESTABLISHED")

    def test_omitted_covariance_terms_are_now_explicit(self):
        for k in (
            "source_generated_not_trajectory_fit","V10_canonical_core_retained",
            "V10_exact_combined_perpendicular_residual_retained",
            "first_attitude_PSD_diagonal_remainder_already_in_V10",
            "first_attitude_PSD_cross_axis_remainder_included",
            "second_prediction_attitude_process_remainder_included",
            "sample1_S_covariance_update_included","sample1_S_attitude_injection_included",
            "sample1_S_aw_mean_correction_included",
            "sample1_S_solver_identity_branch_contained_as_zero_perturbation",
        ):
            self.assertIs(self.d[k],True)

    def test_no_dependency_regression_or_promotion(self):
        for k in (
            "source_replay_used","filter_changed","broad_sample1_3x3_interval_inverse_reintroduced",
            "temporal_force_slew_assumed","complete_sample1_branch_closed_here",
            "signed_cayley_q8_composed_here","q8_word_promoted_here",
            "whole_word_promoted_here","N_H_words_set_here",
        ):
            self.assertIs(self.d[k],False)

    def test_perturbation_is_nonnegative_and_additive(self):
        self.assertTrue(self.d["rows"])
        for r in self.d["rows"]:
            if "V11_correction_norm_upper_rad" not in r:
                continue
            self.assertGreaterEqual(r["total_residual_perturbation_upper_mps2"],0.0)
            self.assertGreaterEqual(r["total_reduced_covariance_perturbation_upper"],0.0)
            self.assertGreaterEqual(r["sample1_gain_operator_perturbation_upper"],0.0)
            self.assertGreaterEqual(
                r["V11_correction_norm_upper_rad"]+1e-12,
                r["V10_directional_correction_upper_rad"],
            )

    def test_fail_closed_numeric_result(self):
        st=self.d["P5_SAMPLE1_PSD_S_PERTURBATION_V11"]
        self.assertIn(st,("PASS","NOT_ESTABLISHED"))
        for k in (
            "max_total_residual_perturbation_upper_mps2",
            "max_total_reduced_covariance_perturbation_upper",
            "max_sample1_gain_operator_perturbation_upper",
            "max_V11_correction_norm_upper_rad",
        ):
            self.assertTrue(math.isfinite(float(self.d[k])))
        self.assertEqual(self.d["first_unclosed_joint_cell"] is None,st=="PASS")
        self.assertEqual(self.d["unclosed_joint_cells"]==0,st=="PASS")
        if st=="PASS":
            self.assertLess(self.d["max_V11_correction_norm_upper_rad"],9.0)
            self.assertEqual(self.d["next_obligation"],"SIGNED_CAYLEY_COMPOSE_SAMPLE1_S_AND_ACCELERATOR_INSIDE_Q8")

    def test_S_branch_bounds_are_finite(self):
        s=self.d["sample1_S_perturbation_bounds"]
        for k in (
            "sample1_S_attitude_correction_upper_rad","sample1_S_aw_correction_upper_mps2",
            "sample1_S_covariance_decrement_upper","sample1_S_reset_covariance_perturbation_upper",
            "sample1_S_total_reduced_covariance_perturbation_upper",
        ):
            self.assertTrue(math.isfinite(float(s[k])))
            self.assertGreaterEqual(float(s[k]),0.0)


if __name__=="__main__": unittest.main()
