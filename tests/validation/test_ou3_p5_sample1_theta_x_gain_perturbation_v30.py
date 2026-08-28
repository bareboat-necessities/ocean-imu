from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_theta_x_gain_perturbation_v30 as V30


class Sample1ThetaXGainPerturbationV30Tests(unittest.TestCase):
    def test_row_resolvent_uses_parallel_nominal_gain(self):
        vr = {
            "sample1_attitude_cross_covariance_perturbation_upper": 0.1,
            "sample1_innovation_perturbation_upper": 0.2,
            "actual_innovation_inverse_operator_upper": 0.5,
            "sample1_attitude_gain_operator_perturbation_upper": 0.16,
        }
        base = {"Ktheta_parallel_block_upper": 0.1}
        d = V30._theta_x_gain_perturbation_upper(vr, base)
        self.assertAlmostEqual(d["theta_x_gain_perturbation_operator_upper"], 0.06, places=12)
        self.assertLess(d["theta_x_gain_perturbation_operator_upper"],
                        d["V12D_full_attitude_gain_perturbation_operator_upper"])

    def test_row_resolvent_refuses_parent_violation(self):
        vr = {
            "sample1_attitude_cross_covariance_perturbation_upper": 1.0,
            "sample1_innovation_perturbation_upper": 1.0,
            "actual_innovation_inverse_operator_upper": 1.0,
            "sample1_attitude_gain_operator_perturbation_upper": 0.1,
        }
        base = {"Ktheta_parallel_block_upper": 1.0}
        with self.assertRaises(RuntimeError):
            V30._theta_x_gain_perturbation_upper(vr, base)

    def test_validation_keeps_promotion_guards(self):
        d = {
            "schema": V30.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_THETA_X_GAIN_PERTURBATION_V30",
            "source_generated_not_trajectory_fit": True,
            "V29_directional_V12D_parent_retained": True,
            "V12D_resolvent_reused_rowwise": True,
            "theta_x_DeltaC_row_bounded_by_full_DeltaC_operator_norm": True,
            "theta_x_nominal_K_row_uses_parallel_block_bound": True,
            "theta_yz_gain_perturbation_parent_unchanged": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "theta_x_gain_perturbation_detail": {
                "theta_x_gain_perturbation_operator_upper": 0.1,
                "V12D_full_attitude_gain_perturbation_operator_upper": 0.2,
            },
            "deployed_correction_limit_rad": 6.0,
            "q_target": 8.0,
            "P5_SAMPLE1_THETA_X_GAIN_PERTURBATION_V30": "PASS",
            "failures": [],
        }
        self.assertEqual(V30.validate(d), [])


if __name__ == "__main__":
    unittest.main()
