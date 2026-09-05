#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import ou3_sea3_hard_window_behavior as BEHAVIOR


class HardWindowBehaviorTest(unittest.TestCase):
    def test_sampled_complete_sea3_behavior_is_compact(self) -> None:
        d = BEHAVIOR.build()
        self.assertEqual(BEHAVIOR.validate(d), [])
        self.assertEqual(d["behavior_set_symbol"], "B^601_SEA3")
        self.assertTrue(d["sampled_behavior_set_compact"])
        self.assertTrue(d["membership_requires_common_SEA3_witness"])
        self.assertTrue(d["continuum_phase_coordinate_set_closed"])
        self.assertTrue(d["phase_continuous_propagation_closed"])
        self.assertTrue(d["machine_readable_R_lambda_closed"])
        self.assertFalse(d["validated_membership_or_separation_oracle_closed"])
        self.assertFalse(d["validated_correlated_outer_enclosure_closed"])
        self.assertFalse(d["P3_promoted"])

    def test_norm_caps_do_not_become_a_source_generator(self) -> None:
        d = BEHAVIOR.build()
        self.assertFalse(d["normal_live_caps_are_membership_sufficient"])
        self.assertFalse(d["arbitrary_bounded_sequence_is_member"])
        self.assertFalse(d["independent_sample_boxes_define_behavior_set"])
        self.assertFalse(d["independent_axis_boxes_define_behavior_set"])
        self.assertFalse(d["finite_frequency_grid_used"])
        self.assertFalse(d["seeded_simulator_used"])
        self.assertFalse(d["gaussian_good_event_used"])
        self.assertFalse(d["spectral_moments_alone_used_as_membership"])


if __name__ == "__main__":
    unittest.main()
