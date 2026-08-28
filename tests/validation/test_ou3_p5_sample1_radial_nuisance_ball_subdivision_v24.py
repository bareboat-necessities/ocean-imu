from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_radial_nuisance_ball_subdivision_v24 as V24


class Sample1RadialNuisanceBallSubdivisionV24Tests(unittest.TestCase):
    def test_structured_gain_operator_norm_uses_larger_orthogonal_block(self):
        parent = {
            "gain_detail": {
                "perpendicular_gain_components": [[3.0, 3.0], [4.0, 4.0]],
                "parallel_gain_components": [[0.0, 0.0], [2.0, 2.0]],
            }
        }
        k = V24._gain_operator_norm(parent)
        self.assertGreaterEqual(k, 5.0)
        self.assertLess(k, 5.0 + 1e-12)

    def test_nominal_exact_correction_is_zero_at_zero_current_rotation(self):
        parent = {
            "sample1_force_components_yz_mps2": [[-1.0, -1.0], [9.0, 9.0]],
            "gain_detail": {
                "perpendicular_gain_components": [[0.1, 0.1], [0.2, 0.2]],
                "parallel_gain_components": [[-0.3, -0.3], [0.4, 0.4]],
            },
        }
        c = [Interval(0.0, 0.0) for _ in range(3)]
        d, y = V24._nominal_exact_correction(c, parent)
        # Outward interval arithmetic may leave a few subnormal ulps around
        # mathematical zero.  The rigorous contract is that zero is enclosed
        # and the outward dust remains negligible, not bit-exact endpoints.
        for x in y + d:
            self.assertLessEqual(x.lo, 0.0)
            self.assertGreaterEqual(x.hi, 0.0)
            self.assertLessEqual(max(abs(x.lo), abs(x.hi)), 1e-320)

    def test_validation_keeps_radial_ball_and_promotion_guards(self):
        d = {
            "schema": V24.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_RADIAL_NUISANCE_BALL_SUBDIVISION_V24",
            "source_generated_not_trajectory_fit": True,
            "V23_current_exact_residual_subdivision_parent_retained": True,
            "structured_gain_operator_norm_used_for_physical_nuisance_ball": True,
            "V12D_correction_perturbation_retained_as_radial_ball": True,
            "componentwise_correction_box_retained_for_sign_and_feasibility": True,
            "independent_radial_ball_intersected_before_q8_test": True,
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
            "candidate_current_subboxes": 64,
            "q_ball_rejected_current_subboxes": 0,
            "component_incompatible_current_subboxes": 0,
            "radial_incompatible_current_subboxes": 0,
            "compatible_current_subboxes": 64,
            "closed_current_subboxes": 60,
            "open_current_subboxes": 4,
            "focused_first_witness_signed_subcell_closed_by_radial_ball": False,
            "P5_SAMPLE1_RADIAL_NUISANCE_BALL_SUBDIVISION_V24": "PASS",
            "failures": [],
        }
        self.assertEqual(V24.validate(d), [])


if __name__ == "__main__":
    unittest.main()
