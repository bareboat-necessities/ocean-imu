import math
import pathlib
import sys
import unittest

TOOLS = pathlib.Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import ou3_p5_sample1_directional_innovation_row_lift_v34 as V34


class Sample1DirectionalInnovationRowLiftV34Tests(unittest.TestCase):
    def test_row_delta_s_bound_uses_row_structure(self):
        got = V34._innovation_row_perturbation_upper(
            h_row_norm=1.0, h_norm=2.0, hp_row_norm=1.5,
            p_norm=3.0, dP=0.01, dH=0.02)
        generic = V34._sum_up(
            2.0 * 2.0 * 0.01,
            2.0 * 2.0 * 3.0 * 0.02,
            3.0 * 0.02 * 0.02,
            2.0 * 2.0 * 0.01 * 0.02,
            0.01 * 0.02 * 0.02,
        )
        self.assertTrue(math.isfinite(got))
        self.assertGreaterEqual(got, 0.0)
        self.assertLess(got, generic)

    def test_directional_caps_intersect_existing_parents(self):
        parent = {
            "gain_perturbation_ball_upper_rad": 0.5,
            "yz_correction_perturbation_norm_upper_rad": 0.6,
            "total_correction_perturbation_norm_upper_rad": 0.7,
            "x_correction_perturbation_abs_upper_rad": 0.2,
        }

        def parent_fn(*, base, vr, row_detail):
            return dict(parent)

        out = V34._directional_caps(
            base={
                "Ktheta_perpendicular_block_upper": 0.1,
                "sample1_full_residual_norm_upper_mps2": 1.0,
            },
            vr={"total_residual_perturbation_upper_mps2": 0.1},
            x_detail={},
            ds_detail={
                "nominal_theta_y_gain_row_norm_upper": 0.05,
                "nominal_theta_z_gain_row_norm_upper": 0.04,
                "theta_y_gain_perturbation_intersected_upper": 0.01,
                "theta_z_gain_perturbation_intersected_upper": 0.02,
            },
            parent_fn=parent_fn,
        )
        self.assertLessEqual(
            out["rowwise_yz_gain_perturbation_norm_upper_rad"],
            parent["gain_perturbation_ball_upper_rad"])
        self.assertLessEqual(
            out["yz_correction_perturbation_norm_upper_rad"],
            parent["yz_correction_perturbation_norm_upper_rad"])
        self.assertLessEqual(
            out["total_correction_perturbation_norm_upper_rad"],
            parent["total_correction_perturbation_norm_upper_rad"])

    def test_validation_keeps_promotion_guards(self):
        d = {
            "schema": V34.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_DIRECTIONAL_INNOVATION_ROW_LIFT_V34",
            "source_generated_not_trajectory_fit": True,
            "V31_current_subbox_lift_parent_retained": True,
            "V32_exact_theta_x_construction_retained": True,
            "V33_invalid_theta_yz_DeltaC_route_retired": True,
            "V12D_full_DeltaC_parent_retained_for_theta_yz": True,
            "exact_sparse_theta_yz_gain_rows_used": True,
            "first_measurement_DeltaS_row_bounded_from_exact_nominal_HP": True,
            "first_measurement_DeltaS_row_intersected_with_V12D_parent": True,
            "theta_yz_gain_row_bounds_intersected_with_V12D_parent": True,
            "theta_yz_correction_bounds_intersected_with_V29_parent": True,
            "all_candidate_current_subboxes_accounted": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "deployed_correction_limit_rad": 6.0,
            "q_target": V34.Q_TARGET,
            "directional_innovation_detail": {
                "V12D_full_innovation_perturbation_upper": 0.2,
                "first_measurement_row_DeltaS_intersected_upper": 0.1,
                "V12D_full_attitude_gain_perturbation_operator_upper": 0.3,
                "theta_y_gain_perturbation_intersected_upper": 0.1,
                "theta_z_gain_perturbation_intersected_upper": 0.1,
                "theta_y_gain_row_exact_scalar_x_channel": True,
                "theta_z_gain_row_exact_scalar_x_channel": True,
                "V12D_full_DeltaC_parent_retained_for_theta_yz": True,
            },
            "P5_SAMPLE1_DIRECTIONAL_INNOVATION_ROW_LIFT_V34": "PASS",
            "failures": [],
        }
        self.assertEqual(V34.validate(d), [])
        d["N_H_words_set_here"] = True
        self.assertIn("N_H_words_set_here is not false", V34.validate(d))


if __name__ == "__main__":
    unittest.main()
