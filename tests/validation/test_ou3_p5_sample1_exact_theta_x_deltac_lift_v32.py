import math
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p5_sample1_exact_theta_x_deltac_lift_v32 as V32


class Sample1ExactThetaXDeltaCLiftV32Tests(unittest.TestCase):
    def test_exact_theta_x_deltac_uses_row_norm(self):
        dP = 0.01
        dH = 0.02
        htheta = 10.0
        a = 0.001
        full_ptheta = 0.2
        exact = V32._theta_x_deltac_upper(
            dP=dP, dH=dH, htheta_norm=htheta, a_row_norm=a)
        generic = V32._sum_up(
            V32.FULL.up(dP * htheta),
            V32.FULL.up(full_ptheta * dH),
            V32.FULL.up(dP * dH),
            dP)
        self.assertLess(exact, generic)

    def test_exact_theta_x_deltac_rejects_negative_inputs(self):
        with self.assertRaises(ValueError):
            V32._theta_x_deltac_upper(
                dP=-1.0, dH=0.1, htheta_norm=1.0, a_row_norm=0.1)

    def test_validation_keeps_promotion_guards(self):
        d = {
            "schema": V32.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_EXACT_THETA_X_DELTAC_LIFT_V32",
            "source_generated_not_trajectory_fit": True,
            "V31_current_subbox_lift_parent_retained": True,
            "V30_theta_x_row_resolvent_retained": True,
            "V12D_full_DeltaC_parent_retained": True,
            "theta_x_nominal_covariance_row_exact_a_0_0_used": True,
            "theta_x_Ptheta_DeltaH_uses_exact_a_row_norm": True,
            "V22_exact_current_residual_parent_retained": True,
            "V28_split_gravity_signed_source_enclosure_retained": True,
            "V29_yz_and_radial_perturbation_parents_retained": True,
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
            "theta_x_exact_DeltaC_gain_detail": {
                "theta_x_DeltaC_operator_upper": 0.001,
                "V12D_full_DeltaC_operator_upper": 0.01,
                "theta_x_gain_perturbation_operator_upper": 0.002,
                "V12D_full_attitude_gain_perturbation_operator_upper": 0.02,
                "theta_x_covariance_row_exactly_a_0_0": True,
            },
            "deployed_correction_limit_rad": 6.0,
            "q_target": V32.Q_TARGET,
            "P5_SAMPLE1_EXACT_THETA_X_DELTAC_LIFT_V32": "PASS",
            "failures": [],
        }
        self.assertEqual(V32.validate(d), [])
        d["q8_word_promoted_here"] = True
        self.assertIn("q8_word_promoted_here is not false", V32.validate(d))


if __name__ == "__main__":
    unittest.main()
