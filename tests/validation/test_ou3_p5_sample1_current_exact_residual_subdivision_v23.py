from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_current_exact_residual_subdivision_v23 as V23


class Sample1CurrentExactResidualSubdivisionV23Tests(unittest.TestCase):
    def test_q_ball_projection_rejects_disjoint_box(self):
        box = [
            Interval.outward_bounds(0.8, 1.0),
            Interval.outward_bounds(0.8, 1.0),
            Interval.outward_bounds(0.8, 1.0),
        ]
        self.assertIsNone(V23._clip_box_to_q_ball(box, 1.0))

    def test_q_ball_projection_tightens_partially_intersecting_box(self):
        box = [
            Interval.outward_bounds(-2.0, 2.0),
            Interval.outward_bounds(0.6, 0.7),
            Interval.outward_bounds(0.6, 0.7),
        ]
        clipped = V23._clip_box_to_q_ball(box, 1.0)
        self.assertIsNotNone(clipped)
        self.assertGreater(clipped[0].lo, box[0].lo)
        self.assertLess(clipped[0].hi, box[0].hi)
        self.assertLess(clipped[0].abs_upper(), 0.54)
        for c, p in zip(clipped, box):
            self.assertGreaterEqual(c.lo, p.lo)
            self.assertLessEqual(c.hi, p.hi)
        self.assertLessEqual(V23.V14.CAYLEY2._norm_lower(clipped), 1.0)

    def test_current_component_partition_is_cartesian(self):
        parent = [
            Interval.outward_bounds(-1.0, 1.0),
            Interval.outward_bounds(-2.0, 2.0),
            Interval.outward_bounds(-3.0, 3.0),
        ]
        boxes = V23._current_subboxes(parent, 2)
        self.assertEqual(len(boxes), 8)
        with self.assertRaises(ValueError):
            V23._current_subboxes(parent, 1)

    def test_validation_keeps_estimator_and_promotion_guards(self):
        d = {
            "schema": V23.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_CURRENT_EXACT_RESIDUAL_SUBDIVISION_V23",
            "source_generated_not_trajectory_fit": True,
            "V22_exact_nonlinear_residual_parent_retained": True,
            "authoritative_current_component_box_partitioned": True,
            "current_q_ball_projected_per_subbox": True,
            "current_subbox_and_exact_residual_kept_joint_through_q8_test": True,
            "V13E_signed_correction_subcell_intersected_per_current_subbox": True,
            "V16_axis_cone_V15_geodesic_V18_yz_support_retained": True,
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
            "candidate_current_subboxes": 8,
            "q_ball_rejected_current_subboxes": 0,
            "evaluated_q_compatible_current_subboxes": 8,
            "source_correction_rejected_current_subboxes": 0,
            "radial_constraint_rejected_current_subboxes": 0,
            "compatible_current_subboxes": 8,
            "closed_current_subboxes": 6,
            "open_current_subboxes": 2,
            "focused_first_witness_signed_subcell_closed_by_subdivision": False,
            "P5_SAMPLE1_CURRENT_EXACT_RESIDUAL_SUBDIVISION_V23": "PASS",
            "failures": [],
        }
        self.assertEqual(V23.validate(d), [])


if __name__ == "__main__":
    unittest.main()
