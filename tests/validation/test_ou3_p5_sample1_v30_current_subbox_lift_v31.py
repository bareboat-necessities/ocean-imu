import math
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p5_sample1_v30_current_subbox_lift_v31 as V31


class Sample1V30CurrentSubboxLiftV31Tests(unittest.TestCase):
    def test_refined_caps_only_shrink_x_gain_direction(self):
        base = {
            "Ktheta_perpendicular_block_upper": 1.1,
            "Ktheta_parallel_block_upper": 0.1,
            "sample1_full_residual_norm_upper_mps2": 2.0,
        }
        vr = {
            "total_residual_perturbation_upper_mps2": 0.02,
            "sample1_attitude_gain_operator_perturbation_upper": 0.01,
        }
        row_detail = {"theta_x_gain_perturbation_operator_upper": 0.001}
        parent = V31.V29._directional_perturbation_caps(
            k_perp=1.1, k_parallel=0.1, drho=0.02, dk=0.01, rho=2.0)
        got = V31._refined_caps(base=base, vr=vr, row_detail=row_detail)
        self.assertLess(got["x_correction_perturbation_abs_upper_rad"],
                        parent["x_correction_perturbation_abs_upper_rad"])
        self.assertEqual(got["yz_correction_perturbation_norm_upper_rad"],
                         parent["yz_correction_perturbation_norm_upper_rad"])
        self.assertEqual(got["total_correction_perturbation_norm_upper_rad"],
                         parent["total_correction_perturbation_norm_upper_rad"])

    def test_four_piece_current_partition_has_64_boxes(self):
        I = V31.Interval
        parent = [I(-1.0, 1.0), I(-2.0, 2.0), I(-3.0, 3.0)]
        boxes = V31.V23._current_subboxes(parent, 4)
        self.assertEqual(len(boxes), 64)

    def test_validation_keeps_promotion_guards(self):
        d = {
            "schema": V31.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_V30_CURRENT_SUBBOX_LIFT_V31",
            "source_generated_not_trajectory_fit": True,
            "V22_exact_current_residual_parent_retained": True,
            "V28_split_gravity_signed_source_enclosure_retained": True,
            "V29_yz_and_radial_perturbation_parents_retained": True,
            "V30_theta_x_row_resolvent_retained": True,
            "V23_current_partition_and_q_ball_projection_retained": True,
            "current_dependent_and_source_directional_correction_enclosures_intersected": True,
            "V16_axis_cone_V15_geodesic_V18_yz_support_retained": True,
            "all_candidate_current_subboxes_accounted": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "current_component_pieces": 4,
            "candidate_current_subboxes": 64,
            "closed_current_subboxes": 64,
            "open_current_subboxes": 0,
            "theta_x_gain_perturbation_detail": {
                "theta_x_gain_perturbation_operator_upper": 0.001,
                "V12D_full_attitude_gain_perturbation_operator_upper": 0.01,
            },
            "directional_perturbation_detail": {
                "x_correction_perturbation_abs_upper_rad": 0.03,
                "yz_correction_perturbation_norm_upper_rad": 0.2,
                "total_correction_perturbation_norm_upper_rad": 0.2,
            },
            "previous_isotropic_V12D_correction_perturbation_upper_rad": 0.21,
            "deployed_correction_limit_rad": 6.0,
            "q_target": V31.Q_TARGET,
            "P5_SAMPLE1_V30_CURRENT_SUBBOX_LIFT_V31": "PASS",
            "failures": [],
        }
        self.assertEqual(V31.validate(d), [])
        d["N_H_words_set_here"] = True
        self.assertIn("N_H_words_set_here is not false", V31.validate(d))


if __name__ == "__main__":
    unittest.main()
