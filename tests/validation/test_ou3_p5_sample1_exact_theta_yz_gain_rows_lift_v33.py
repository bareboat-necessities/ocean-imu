import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import Interval
import ou3_p5_sample1_exact_theta_yz_gain_rows_lift_v33 as V33


class Sample1ExactThetaYZGainRowsLiftV33Tests(unittest.TestCase):
    def test_row_norm_upper(self):
        r = [Interval(-1.0, 1.0), Interval(2.0, 2.0), Interval(0.0, 0.0)]
        self.assertGreaterEqual(V33._row_norm_upper(r), 5.0 ** 0.5)

    def test_rowwise_yz_caps_intersect_parent(self):
        base = {
            "Ktheta_perpendicular_block_upper": 1.0,
            "sample1_full_residual_norm_upper_mps2": 2.0,
        }
        vr = {"total_residual_perturbation_upper_mps2": 0.1}
        yz_detail = {
            "nominal_theta_y_gain_row_norm_upper": 0.3,
            "nominal_theta_z_gain_row_norm_upper": 0.4,
            "theta_y_gain_perturbation_operator_upper": 0.01,
            "theta_z_gain_perturbation_operator_upper": 0.02,
        }
        def parent_fn(*, base, vr, row_detail):
            return {
                "gain_perturbation_ball_upper_rad": 0.2,
                "x_correction_perturbation_abs_upper_rad": 0.05,
                "yz_correction_perturbation_norm_upper_rad": 0.3,
                "total_correction_perturbation_norm_upper_rad": 0.31,
            }
        got = V33._yz_refined_caps(
            base=base, vr=vr, x_detail={},
            yz_detail=yz_detail, parent_fn=parent_fn)
        self.assertLessEqual(got["yz_correction_perturbation_norm_upper_rad"], 0.3)
        self.assertLessEqual(got["total_correction_perturbation_norm_upper_rad"], 0.31)

    def test_validation_keeps_promotion_guards(self):
        d = {
            "schema": V33.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_EXACT_THETA_YZ_GAIN_ROWS_LIFT_V33",
            "source_generated_not_trajectory_fit": True,
            "V31_current_subbox_lift_parent_retained": True,
            "V32_exact_theta_x_construction_retained": True,
            "exact_sparse_theta_gain_row_structure_used": True,
            "theta_yz_nominal_covariance_rows_bounded_separately": True,
            "theta_yz_gain_perturbation_rows_bounded_separately": True,
            "theta_yz_rowwise_bounds_intersect_full_operator_parent": True,
            "all_candidate_current_subboxes_accounted": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "theta_x_exact_gain_detail": {
                "theta_x_gain_perturbation_operator_upper": 0.01,
                "V12D_full_attitude_gain_perturbation_operator_upper": 0.03,
            },
            "theta_yz_exact_gain_detail": {
                "theta_y_gain_perturbation_operator_upper": 0.01,
                "theta_z_gain_perturbation_operator_upper": 0.02,
                "V12D_full_attitude_gain_perturbation_operator_upper": 0.03,
                "theta_y_gain_row_exact_scalar_x_channel": True,
                "theta_z_gain_row_exact_scalar_x_channel": True,
            },
            "deployed_correction_limit_rad": 6.0,
            "q_target": V33.Q_TARGET,
            "P5_SAMPLE1_EXACT_THETA_YZ_GAIN_ROWS_LIFT_V33": "PASS",
            "failures": [],
        }
        self.assertEqual(V33.validate(d), [])
        d["whole_word_promoted_here"] = True
        self.assertIn("whole_word_promoted_here is not false", V33.validate(d))


if __name__ == "__main__":
    unittest.main()
