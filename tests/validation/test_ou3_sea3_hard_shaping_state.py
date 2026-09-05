#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import ou3_sea3_hard_shaping_state as SHAPING


class HardShapingStateContractTest(unittest.TestCase):
    def test_contract_is_valid_but_not_promoted(self) -> None:
        d = SHAPING.build()
        self.assertEqual(SHAPING.validate(d), [])
        self.assertTrue(d["SEA3_parameter_domain_compact"])
        self.assertTrue(d["compactness_is_not_an_open_obligation"])
        self.assertTrue(d["theorem_has_deterministic_shaping_contract"])
        self.assertTrue(d["theorem_separates_probabilistic_random_sea_corollary"])
        self.assertTrue(d["complete_source_rejects_gaussian_word_generator"])
        self.assertFalse(d["hard_shaping_state_or_excitation_bound_closed"])
        self.assertFalse(d["complete_SEA3_family_materialized_here"])
        self.assertFalse(d["P3_promoted"])

    def test_forbidden_surrogates_cannot_close_xs(self) -> None:
        d = SHAPING.build()
        for key in (
            "power_spectrum_alone_is_hard_pathwise_bound",
            "spectral_moments_alone_may_close_xs",
            "gaussian_good_event_may_close_xs",
            "replay_may_close_xs",
            "seeded_128_frequency_generator_may_close_xs",
            "finite_RAO_grid_may_close_xs",
            "arbitrary_bounded_input_box_may_close_xs",
        ):
            self.assertFalse(d[key], key)

    def test_every_executable_ingredient_remains_code_owned_false(self) -> None:
        d = SHAPING.build()
        self.assertEqual(
            d["executable_ingredients"],
            {
                "hard_compact_phase_or_driver_set_closed": False,
                "phase_continuous_propagation_closed": False,
                "complete_SEA3_left_inclusion_closed": False,
                "joint_source_output_map_closed": False,
            },
        )


if __name__ == "__main__":
    unittest.main()
