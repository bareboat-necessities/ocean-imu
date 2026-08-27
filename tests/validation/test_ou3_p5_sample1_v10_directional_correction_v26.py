from __future__ import annotations
import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_v10_directional_correction_v26 as V26


class Sample1V10DirectionalCorrectionV26Tests(unittest.TestCase):
    def test_directional_caps_preserve_one_plus_two_assignment(self):
        row = {
            "sample1_full_residual_norm_upper_mps2": 10.0,
            "sample1_combined_source_x_residual_upper_mps2": 2.0,
            "Ktheta_perpendicular_block_upper": 3.0,
            "Ktheta_parallel_block_upper": 0.5,
            "combined_directional_correction_norm_upper_rad": 7.0,
        }
        d = V26._directional_caps(row, 0.25)
        self.assertGreaterEqual(d["dx_abs_upper_rad"], 5.25)
        self.assertLess(d["dx_abs_upper_rad"], 5.25 + 1e-12)
        self.assertGreaterEqual(d["dyz_norm_upper_rad"], 6.25)
        self.assertLess(d["dyz_norm_upper_rad"], 6.25 + 1e-12)
        self.assertGreaterEqual(d["radial_upper_rad"], 7.25)
        self.assertLess(d["radial_upper_rad"], 7.25 + 1e-12)

    def test_yz_ball_projection_rejects_disjoint_box(self):
        self.assertIsNone(V26._clip_yz_ball(Interval(2.0, 3.0), Interval(2.0, 3.0), 1.0))

    def test_yz_ball_projection_tightens_partial_box(self):
        yz = V26._clip_yz_ball(Interval(-2.0, 2.0), Interval(0.8, 2.0), 1.0)
        self.assertIsNotNone(yz)
        y, z = yz
        self.assertLess(y.hi, 0.61)
        self.assertGreater(y.hi, 0.59)
        self.assertLessEqual(z.hi, 1.0 + 8.0 * math.ulp(1.0))

    def test_validation_keeps_source_and_promotion_guards(self):
        d = {
            "schema": V26.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_V10_DIRECTIONAL_FIRST_OPEN_WITNESS_V26",
            "source_generated_not_trajectory_fit": True,
            "V23_first_open_subbox_retained": True,
            "V10_combined_perpendicular_residual_identity_revalidated": True,
            "V10_one_plus_two_directional_caps_used": True,
            "V12D_correction_perturbation_retained_as_single_ball": True,
            "V16_axis_cone_and_V18_signed_product_retained": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "deployed_correction_limit_rad": 6.0,
            "q_target": 8.0,
            "V10_directional_source_detail": {
                "dx_abs_upper_rad": 0.4,
                "dyz_norm_upper_rad": 0.7,
                "radial_upper_rad": 0.8,
            },
            "directional_constraints_incompatible": False,
            "directional_correction_box_rad": [[-0.4, 0.4], [-0.3, 0.3], [-0.3, 0.3]],
            "P5_SAMPLE1_V10_DIRECTIONAL_FIRST_OPEN_WITNESS_V26": "PASS",
            "failures": [],
        }
        self.assertEqual(V26.validate(d), [])


if __name__ == "__main__":
    unittest.main()
