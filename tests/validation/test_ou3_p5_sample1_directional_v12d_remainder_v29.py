from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_directional_v12d_remainder_v29 as V29


class Sample1DirectionalV12DRemainderV29Tests(unittest.TestCase):
    def test_directional_caps_keep_parallel_x_smaller(self):
        d = V29._directional_perturbation_caps(
            k_perp=1.0, k_parallel=0.1, drho=0.2, dk=0.01, rho=2.0)
        self.assertLess(d["x_correction_perturbation_abs_upper_rad"],
                        d["yz_correction_perturbation_norm_upper_rad"])
        self.assertLessEqual(d["yz_correction_perturbation_norm_upper_rad"],
                             d["total_correction_perturbation_norm_upper_rad"] + 1e-15)

    def test_yz_projection_rejects_disjoint_box(self):
        y = Interval(2.0, 3.0); z = Interval(2.0, 3.0)
        self.assertIsNone(V29._clip_yz_to_radius(y, z, 1.0))

    def test_yz_projection_tightens(self):
        y = Interval(-2.0, 2.0); z = Interval(0.9, 1.1)
        out = V29._clip_yz_to_radius(y, z, 1.2)
        self.assertIsNotNone(out)
        yy, zz = out
        self.assertLess(yy.hi, 2.0)
        self.assertGreaterEqual(zz.lo, 0.9)

    def test_validation_keeps_promotion_guards(self):
        d = {
            "schema": V29.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_DIRECTIONAL_V12D_REMAINDER_V29",
            "source_generated_not_trajectory_fit": True,
            "V28_split_gravity_parent_retained": True,
            "V12D_exact_deltaK_deltaR_decomposition_retained": True,
            "V10_one_plus_two_orthogonal_gain_blocks_used": True,
            "V12D_residual_perturbation_mapped_directionally": True,
            "V12D_gain_perturbation_retained_as_radial_ball": True,
            "V23_first_open_current_subbox_retained": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "directional_perturbation_detail": {
                "total_correction_perturbation_norm_upper_rad": 0.2,
                "x_correction_perturbation_abs_upper_rad": 0.1,
                "yz_correction_perturbation_norm_upper_rad": 0.2,
            },
            "previous_isotropic_V12D_correction_perturbation_upper_rad": 0.2,
            "source_constraints_incompatible": False,
            "joint_directional_correction_box_rad": [[0,1],[0,1],[0,1]],
            "deployed_correction_limit_rad": 6.0,
            "q_target": 8.0,
            "P5_SAMPLE1_DIRECTIONAL_V12D_REMAINDER_V29": "PASS",
            "failures": [],
        }
        self.assertEqual(V29.validate(d), [])


if __name__ == "__main__":
    unittest.main()
