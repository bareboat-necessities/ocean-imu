from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_v10_gain_radial_ball_v25 as V25


class Sample1V10GainRadialBallV25Tests(unittest.TestCase):
    def test_certified_v10_gain_norm_selects_exact_witness_row(self):
        core = {"rows": [
            {
                "p_cell": 0, "tangent_residual_cell": 0, "axial_residual_cell": 18,
                "Ktheta_perpendicular_block_upper": 9.0,
                "Ktheta_parallel_block_upper": 9.0,
            },
            {
                "p_cell": 0, "tangent_residual_cell": 0, "axial_residual_cell": 19,
                "Ktheta_perpendicular_block_upper": 0.25,
                "Ktheta_parallel_block_upper": 0.4,
            },
        ]}
        k, detail = V25._certified_v10_gain_norm(core)
        self.assertEqual(k, 0.4)
        self.assertEqual(detail["axial_residual_cell"], 19)

    def test_validation_keeps_v10_and_promotion_guards(self):
        d = {
            "schema": V25.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_V10_GAIN_RADIAL_BALL_SUBDIVISION_V25",
            "source_generated_not_trajectory_fit": True,
            "V24_radial_nuisance_ball_parent_retained": True,
            "V10_structured_gain_row_revalidated": True,
            "V10_orthogonal_block_operator_norm_used": True,
            "V24_component_gain_norm_reconstruction_superseded": True,
            "all_candidate_current_subboxes_accounted": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "deployed_correction_limit_rad": 6.0,
            "q_target": 8.0,
            "V10_structured_gain_operator_norm_upper": 0.4,
            "gain_operator_norm_upper": 0.4,
            "open_current_subboxes": 2,
            "focused_first_witness_signed_subcell_closed_by_V10_gain_ball": False,
            "P5_SAMPLE1_V10_GAIN_RADIAL_BALL_SUBDIVISION_V25": "PASS",
            "failures": [],
        }
        self.assertEqual(V25.validate(d), [])


if __name__ == "__main__":
    unittest.main()
