#!/usr/bin/env python3
from __future__ import annotations

import sys
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
TOOLS = REPO / "tools" / "stability"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_hard_spectral_driver as DRIVER


class HardSpectralDriverTests(unittest.TestCase):
    def test_complete_sea3_continuum_driver_closes_only_semantic_source_definition(self):
        d = DRIVER.build()
        self.assertEqual(DRIVER.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(d["hard_spectral_driver_set_closed"])
        self.assertTrue(d["joint_source_output_map_semantics_closed"])
        self.assertTrue(d["conditional_source_definition_closed"])
        self.assertFalse(d["validated_numerical_gram_oracle_closed"])
        self.assertFalse(d["physical_left_inclusion_closed"])
        self.assertFalse(d["P3_changed"])
        self.assertFalse(d["P3_promoted_from_this_module"])
        self.assertFalse(d["P4_promoted_from_this_module"])

    def test_one_common_continuum_coordinate_preserves_time_and_axis_correlation(self):
        d = DRIVER.build()
        source = d["source_definition"]
        op = d["sampled_operator"]
        self.assertTrue(source["same_coordinate_drives_translation_and_rotation"])
        self.assertTrue(source["same_coordinate_drives_all_601_samples"])
        self.assertTrue(source["lambda_history_remains_coupled"])
        self.assertTrue(source["response_witness_remains_coupled"])
        self.assertTrue(op["cross_time_blocks_preserved"])
        self.assertTrue(op["cross_axis_blocks_preserved"])
        self.assertTrue(op["common_witness_required_for_every_block"])
        self.assertTrue(op["not_cartesian_product_of_sample_caps"])

    def test_no_forbidden_surrogate_is_reintroduced(self):
        d = DRIVER.build()
        c = d["continuum_semantics"]
        for key in (
            "finite_frequency_grid_used",
            "finite_direction_grid_used",
            "finite_seeded_harmonic_used",
            "trajectory_replay_used",
            "gaussian_good_event_used",
            "spectral_moment_only_membership_used",
            "independent_axis_boxes_used",
            "independent_sample_boxes_used",
        ):
            self.assertFalse(c[key], key)
        self.assertFalse(
            d["bounded_operator_basis"]["pointwise_Normal_Live_caps_generate_source"]
        )


if __name__ == "__main__":
    unittest.main()
