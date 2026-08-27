from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_base_row_direction_v20 as V20


class Sample1BaseRowDirectionV20Tests(unittest.TestCase):
    def test_exact_gravity_residual_matches_shipping_small_angle_orientation(self):
        z = Interval.point(0.0)
        c = [Interval.point(0.2), z, z]
        r = V20._gravity_residual_from_cayley(c, 9.80665)
        self.assertLessEqual(r[0].lo, 0.0)
        self.assertGreaterEqual(r[0].hi, 0.0)
        self.assertLess(r[1].hi, 0.0)
        self.assertLessEqual(r[2].hi, 0.0)

    def test_base_row_direction_rejects_opposite_tangent_branch_with_small_nuisance(self):
        z = Interval.point(0.0)
        good = V20._gravity_residual_from_cayley(
            [Interval.point(0.2), z, z], 9.80665)
        rt = Interval.outward_bounds(-good[1].hi, -good[1].lo)
        rz = good[2]
        target = [z, Interval.outward_bounds(-rt.hi, -rt.lo), rz]
        self.assertEqual(V20._vector_box_distance_lower(good, target), 0.0)

        opposite = V20._gravity_residual_from_cayley(
            [Interval.point(-0.2), z, z], 9.80665)
        self.assertGreater(V20._vector_box_distance_lower(opposite, target), 0.1)

    def test_direction_subdivision_retains_ball_compatible_cover(self):
        boxes, rejected_tangent, rejected_full, grid = V20._candidate_cayley_boxes(
            2.0, 0.5, tangent_direction_pieces=4, yaw_direction_pieces=4)
        self.assertGreater(len(boxes), 0)
        self.assertEqual(grid, 64)
        self.assertGreaterEqual(rejected_tangent, 0)
        self.assertGreaterEqual(rejected_full, 0)
        self.assertLessEqual(len(boxes), grid)

    def test_validation_forbids_q8_or_word_promotion(self):
        d = {
            "schema": V20.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_SOURCE_CORRELATED_BASE_ROW_DIRECTION_V20",
            "source_generated_not_trajectory_fit": True,
            "V10_direct_first_residual_coordinate_prerequisite_retained": True,
            "V10_canonical_SO2_gravity_gauge_retained": True,
            "positive_x_correction_uses_negative_y_tangent_residual_branch": True,
            "exact_cayley_gravity_residual_used": True,
            "first_aw_and_bias_nuisance_combined_before_direction_rejection": True,
            "P1_full_cayley_ball_retained": True,
            "P1_gravity_tangent_cayley_ball_retained": True,
            "base_row_direction_subdivision_is_source_complete": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "deployed_correction_limit_rad": 6.0,
            "evaluated_base_direction_rows": 1,
            "directionally_refined_base_direction_rows": 1,
            "source_incompatible_base_direction_rows": 0,
            "total_direction_compatibility_checks": 4,
            "total_direction_incompatible_subboxes": 1,
            "minimum_survival_fraction": 0.75,
            "maximum_survival_fraction": 0.75,
            "first_v19_q8_witness_base_direction_row": {"directionally_refined": True},
            "P5_SAMPLE1_SOURCE_CORRELATED_BASE_ROW_DIRECTION_V20": "PASS",
            "failures": [],
        }
        self.assertEqual(V20.validate(d), [])


if __name__ == "__main__":
    unittest.main()
