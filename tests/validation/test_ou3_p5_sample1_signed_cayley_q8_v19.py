from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_signed_cayley_q8_v19 as G


class Sample1JointYzDirectionSubdivisionV19Tests(unittest.TestCase):
    def test_two_piece_partition_covers_parent(self):
        parent = Interval(-3.25, 7.5)
        parts = G._parts2(parent)
        self.assertEqual(len(parts), 2)
        self.assertLessEqual(parts[0].lo, parent.lo)
        self.assertGreaterEqual(parts[-1].hi, parent.hi)
        self.assertLessEqual(parts[0].hi, parts[1].lo)

    def test_joint_subdivision_can_strictly_tighten_global_parent(self):
        vd = [
            Interval.point(0.0),
            Interval(-2.0, 1.0),
            Interval(-1.0, 2.0),
        ]
        chart = {
            "cy": Interval(-2.0, -1.5),
            "cz": Interval(-2.0, -1.5),
            "cyz_norm_upper": 2.0,
        }
        parent = Interval(-5.7, 5.7)
        joint, pairs, empty = G._joint_yz_dot_subdivision(vd, chart, parent)
        self.assertEqual(pairs, 16)
        self.assertGreaterEqual(empty, 0)
        self.assertGreater(joint.lo, parent.lo)
        self.assertLess(joint.hi, parent.hi)
        self.assertLessEqual(joint.lo, joint.hi)

    def test_validate_keeps_shipping_limit_and_promotion_guards(self):
        d = {
            "schema": G.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_JOINT_YZ_DIRECTION_SUBDIVISION_Q8_V19",
            "source_generated_not_trajectory_fit": True,
            "V18B_signed_full_angle_parent_retained": True,
            "joint_current_correction_yz_subdivision_used": True,
            "correction_quaternion_vector_unit_ball_retained": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_rad": 6.0,
            "deployed_correction_limit_increased": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "joint_yz_subdivision_calls": 10,
            "joint_yz_subdivision_attempted_cells": 4,
            "joint_yz_subdivision_skipped_parent_closed_cells": 5,
            "joint_yz_subdivision_skipped_high_q_cells": 1,
            "joint_yz_subdivision_pair_evaluations": 64,
            "joint_yz_subdivision_refined_cells": 3,
            "joint_yz_subdivision_newly_closed_cells": 1,
            "P5_SAMPLE1_JOINT_YZ_DIRECTION_SUBDIVISION_Q8_V19": "NOT_ESTABLISHED",
            "first_unclosed_q8_cell": {"post_sample1_cayley_norm_upper": 9.0},
            "failures": [],
        }
        self.assertEqual(G.validate(d), [])


if __name__ == "__main__":
    unittest.main()
