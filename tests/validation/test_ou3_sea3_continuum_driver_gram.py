#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(TOOLS / "stability"))

import ou3_sea3_continuum_driver_gram as DRIVER


class CompleteSea3ContinuumDriverGramTest(unittest.TestCase):
    def test_contract_closes_driver_definition_but_not_complete_family_enclosure(self) -> None:
        d = DRIVER.build()
        self.assertEqual(DRIVER.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertEqual(d["hard_realization_set_symbol"], "X^s_SEA3(lambda_{0:N_W})")
        self.assertEqual(d["sample_count"], 601)
        self.assertEqual(d["sampled_source_core_dimension"], 3606)
        self.assertTrue(d["hard_spectral_driver_set_closed"])
        self.assertTrue(d["exact_fixed_history_correlated_oracle_formula_closed"])
        self.assertTrue(d["same_driver_field_entire_window"])
        self.assertTrue(d["same_driver_field_translation_and_rotation"])
        self.assertFalse(d["validated_complete_family_gram_enclosure_closed"])
        self.assertFalse(d["provider_word_materialized_here"])
        self.assertFalse(d["P4_promoted_here"])
        self.assertFalse(d["P5_promoted_here"])

    def test_forbidden_shortcuts_are_absent(self) -> None:
        d = DRIVER.build()
        for key in (
            "finite_frequency_grid_used",
            "finite_direction_grid_used",
            "finite_harmonic_source_used",
            "trajectory_replay_used",
            "gaussian_good_event_used",
            "independent_sample_boxes_used",
            "independent_translation_rotation_sources_used",
        ):
            self.assertFalse(d[key], key)

    def test_fixed_history_oracle_is_one_correlated_gram_not_axis_boxes(self) -> None:
        d = DRIVER.build()
        self.assertEqual(d["fixed_history_operator"], "y=K_{lambda,G} a")
        self.assertEqual(d["fixed_history_gram"], "Q_{lambda,G}=K_{lambda,G} K_{lambda,G}^*")
        self.assertIn("Q^dagger", d["fixed_history_membership"])
        self.assertIn("one common driver field", d["complete_family"])


if __name__ == "__main__":
    unittest.main()
