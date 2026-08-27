from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_signed_post_first_aw_components_v27 as V27


class Sample1SignedPostFirstAwComponentsV27Tests(unittest.TestCase):
    def test_v10_witness_row_is_exact(self):
        row = {"p_cell": 0, "tangent_residual_cell": 0, "axial_residual_cell": 19}
        self.assertIs(V27._v10_witness_row({"rows": [row]}), row)

    def test_nominal_signed_gain_mapping(self):
        parent = {
            "gain_detail": {
                "perpendicular_gain_components": [[2.0, 2.0], [3.0, 3.0]],
                "parallel_gain_components": [[4.0, 4.0], [5.0, 5.0]],
            }
        }
        r = [Interval(1.0, 1.0), Interval(2.0, 2.0), Interval(3.0, 3.0)]
        d = V27._nominal_correction(r, parent)
        self.assertLessEqual(d[0].lo, 23.0); self.assertGreaterEqual(d[0].hi, 23.0)
        self.assertLessEqual(d[1].lo, 2.0); self.assertGreaterEqual(d[1].hi, 2.0)
        self.assertLessEqual(d[2].lo, 3.0); self.assertGreaterEqual(d[2].hi, 3.0)

    def test_validation_keeps_promotion_guards(self):
        d = {
            "schema": V27.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_SIGNED_POST_FIRST_AW_COMPONENTS_V27",
            "source_generated_not_trajectory_fit": True,
            "V23_first_open_subbox_retained": True,
            "V10_exact_first_update_OU_cancellation_used": True,
            "signed_tangent_axial_first_residual_cell_retained": True,
            "V21_signed_one_plus_two_gain_components_used": True,
            "V10_transport_series_and_bias_remainders_retained": True,
            "V12D_correction_perturbation_retained_as_single_ball": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "deployed_correction_limit_rad": 6.0,
            "q_target": 8.0,
            "sample1_signed_residual_box_mps2": [[-1,1],[-1,1],[-1,1]],
            "nominal_signed_correction_box_rad": [[-1,1],[-1,1],[-1,1]],
            "source_constraints_incompatible": False,
            "joint_signed_correction_box_rad": [[-1,1],[-1,1],[-1,1]],
            "P5_SAMPLE1_SIGNED_POST_FIRST_AW_COMPONENTS_V27": "PASS",
            "failures": [],
        }
        self.assertEqual(V27.validate(d), [])


if __name__ == "__main__":
    unittest.main()
