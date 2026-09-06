#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))
sys.path.insert(0, str(TOOLS / "stability"))

import ou3_sea3_hard_shaping_state as SHAPING


class HardShapingStateContractTest(unittest.TestCase):
    def test_contract_is_valid_but_not_promoted(self) -> None:
        d = SHAPING.build()
        self.assertEqual(SHAPING.validate(d), [])
        self.assertTrue(d["SEA3_parameter_domain_compact"])
        self.assertTrue(d["compactness_is_not_an_open_obligation"])
        self.assertTrue(d["theorem_has_deterministic_shaping_contract"])
        self.assertTrue(d["theorem_has_explicit_hard_realization_set"])
        self.assertTrue(d["theorem_rejects_statistical_or_seeded_surrogates"])
        self.assertTrue(d["theorem_separates_probabilistic_random_sea_corollary"])
        self.assertTrue(d["complete_source_rejects_gaussian_word_generator"])
        self.assertEqual(d["hard_realization_set_symbol"], "X^s_SEA3(lambda_{0:N_W})")

        phase = d["continuum_phase_certificate"]
        self.assertTrue(phase["continuum_phase_coordinate_set_closed"])
        self.assertTrue(phase["phase_continuous_propagation_closed"])
        self.assertFalse(phase["finite_frequency_grid_used"])
        self.assertFalse(phase["finite_direction_grid_used"])
        self.assertFalse(phase["phase_reset_on_lambda_transition_allowed"])

        driver = d["hard_driver_certificate"]
        self.assertTrue(driver["same_driver_field_entire_window"])
        self.assertTrue(driver["same_driver_field_translation_and_rotation"])
        self.assertTrue(driver["hard_spectral_driver_set_closed"])
        self.assertTrue(driver["exact_fixed_history_correlated_oracle_formula_closed"])
        self.assertFalse(driver["validated_complete_family_gram_enclosure_closed"])

        behavior = d["sampled_behavior_target"]
        self.assertEqual(behavior["symbol"], "B^601_SEA3")
        self.assertTrue(behavior["compact"])
        self.assertTrue(behavior["membership_requires_common_SEA3_witness"])
        self.assertFalse(behavior["normal_live_caps_are_membership_sufficient"])
        self.assertFalse(behavior["independent_sample_boxes_define_behavior_set"])
        self.assertFalse(behavior["validated_membership_or_separation_oracle_closed"])
        self.assertFalse(behavior["validated_correlated_outer_enclosure_closed"])

        self.assertTrue(d["hard_shaping_state_or_excitation_bound_closed"])
        self.assertFalse(d["conditional_SEA3_source_output_closed"])
        self.assertFalse(d["deployment_left_inclusion_closed"])
        self.assertFalse(d["complete_SEA3_family_materialized_here"])
        self.assertFalse(d["P3_promoted"])

    def test_conditional_execution_still_needs_complete_family_output_enclosure(self) -> None:
        d = SHAPING.build()
        self.assertFalse(d["conditional_execution_requires_deployment_left_inclusion"])
        self.assertFalse(SHAPING.COMPLETE_SEA3_LEFT_INCLUSION_CLOSED)
        self.assertFalse(SHAPING.DEPLOYMENT_LEFT_INCLUSION_CLOSED)
        self.assertTrue(SHAPING.HARD_SPECTRAL_DRIVER_SET_CLOSED)
        self.assertTrue(SHAPING.HARD_SHAPING_STATE_OR_EXCITATION_BOUND_CLOSED)
        self.assertFalse(SHAPING.JOINT_SOURCE_OUTPUT_MAP_CLOSED)
        self.assertFalse(SHAPING.CONDITIONAL_SEA3_SOURCE_OUTPUT_CLOSED)

    def test_forbidden_surrogates_remain_forbidden(self) -> None:
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

    def test_driver_is_closed_inside_complete_sea3_but_output_enclosure_is_open(self) -> None:
        d = SHAPING.build()
        self.assertEqual(
            d["executable_ingredients"],
            {
                "continuum_phase_coordinate_set_closed": True,
                "phase_continuous_propagation_closed": True,
                "hard_spectral_driver_set_closed": True,
                "complete_SEA3_left_inclusion_closed": False,
                "joint_source_output_map_closed": False,
            },
        )
        self.assertTrue(d["hard_shaping_state_or_excitation_bound_closed"])
        self.assertFalse(d["conditional_SEA3_source_output_closed"])
        self.assertFalse(d["deployment_left_inclusion_closed"])


if __name__ == "__main__":
    unittest.main()
