from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_structured_full_gain_v12 as V12


class StructuredFullGainV12Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Semantic/fail-closed fixture only; the workflow's 24^3 emitter is
        # the authoritative V12 certificate.
        cls.d=V12.build(source_pieces=4,source_cell_index=0,p_pieces=4,tangent_pieces=4,axial_pieces=4)

    def test_interval_ldlt_lambda_lower_on_diagonal_spd(self):
        I=V12.Interval.outward_bounds
        z=V12.FULL.I(0.0)
        S=[[I(2.0,2.0),z,z],[z,I(3.0,3.0),z],[z,z,I(5.0,5.0)]]
        lam=V12._nominal_lambda_min_lower(S)
        self.assertGreater(lam,1.999999999999)
        self.assertLessEqual(lam,2.0)

    def test_validate_fail_closed_on_the_coarse_grid(self):
        # V10 now closes on this coarse grid, so no prerequisite failure is
        # raised here; V12's own joint cells still do not close on 4^3, and the
        # certificate has to stay NOT_ESTABLISHED with a recorded witness
        # instead of inheriting the parent's PASS.
        failures=V12.validate(self.d)
        self.assertEqual(failures,[])
        self.assertEqual(self.d["P5_SAMPLE1_PSD_S_ACTUAL_INNOVATION_V12"],"NOT_ESTABLISHED")
        self.assertGreater(self.d["unclosed_joint_cells"],0)
        self.assertIsNotNone(self.d["first_unclosed_joint_cell"])

    def test_actual_nominal_innovation_replaces_noise_only_floor(self):
        self.assertIs(self.d["actual_nominal_innovation_lambda_certified_by_interval_LDLT"],True)
        self.assertIs(self.d["measurement_noise_only_resolvent_floor_used"],False)
        self.assertGreaterEqual(
            float(self.d["minimum_nominal_innovation_lambda_lower"])+1e-15,
            0.0,
        )
        for r in self.d["rows"]:
            if "nominal_innovation_lambda_lower" not in r or "measurement_noise_floor_lower" not in r:
                continue
            self.assertGreaterEqual(
                float(r["nominal_innovation_lambda_lower"])+1e-12,
                float(r["measurement_noise_floor_lower"]),
            )

    def test_V11_perturbations_retained(self):
        for k in (
            "V10_canonical_core_retained","V11_PSD_and_S_perturbation_magnitudes_retained",
            "first_attitude_PSD_cross_axis_remainder_included",
            "second_prediction_attitude_process_remainder_included",
            "sample1_S_covariance_update_included","sample1_S_attitude_injection_included",
            "sample1_S_aw_mean_correction_included",
            "sample1_S_solver_identity_branch_contained_as_zero_perturbation",
        ):
            self.assertIs(self.d[k],True)

    def test_no_promotion_or_limit_change(self):
        for k in (
            "source_replay_used","filter_changed","broad_sample1_3x3_interval_inverse_reintroduced",
            "temporal_force_slew_assumed","complete_sample1_branch_closed_here",
            "signed_cayley_q8_composed_here","q8_word_promoted_here",
            "whole_word_promoted_here","N_H_words_set_here","deployed_correction_limit_increased",
        ):
            self.assertIs(self.d[k],False)
        self.assertEqual(float(self.d["deployed_correction_limit_rad"]),6.0)

    def test_perturbed_result_is_fail_closed_and_additive(self):
        st=self.d["P5_SAMPLE1_PSD_S_ACTUAL_INNOVATION_V12"]
        self.assertIn(st,("PASS","NOT_ESTABLISHED"))
        for k in (
            "minimum_nominal_innovation_lambda_lower","max_total_residual_perturbation_upper_mps2",
            "max_total_reduced_covariance_perturbation_upper","max_sample1_gain_operator_perturbation_upper",
            "max_V12_correction_norm_upper_rad",
        ):
            self.assertTrue(math.isfinite(float(self.d[k])))
        for r in self.d["rows"]:
            if "V12_correction_norm_upper_rad" in r:
                self.assertGreaterEqual(
                    float(r["V12_correction_norm_upper_rad"])+1e-12,
                    float(r["V10_directional_correction_upper_rad"]),
                )
        self.assertEqual(self.d["first_unclosed_joint_cell"] is None,st=="PASS")
        self.assertEqual(self.d["unclosed_joint_cells"]==0,st=="PASS")
        if st=="PASS":
            self.assertLess(float(self.d["max_V12_correction_norm_upper_rad"]),9.0)
            self.assertEqual(self.d["next_obligation"],"SIGNED_CAYLEY_COMPOSE_SAMPLE1_S_AND_ACCELERATOR_INSIDE_Q8")


if __name__=="__main__": unittest.main()
