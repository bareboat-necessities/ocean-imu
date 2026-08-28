#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_v40_full_source_cell0_q8_lift_v41 as V41


class Sample1V40FullSourceCell0Q8LiftV41Tests(unittest.TestCase):
    def test_build_installs_v40_helper_and_restores_it(self):
        original_hook = V41.V12D._first_psd_perturbation_tangent
        original_refined = V41.V40._first_psd_perturbation_exact_joseph_components
        original_build = V41.V18B.build
        original_validate = V41.V18B.validate
        seen = {}

        def fake_refined(**kwargs):
            seen["kwargs"] = kwargs
            return {
                "first_offaxis_attitude_correction_upper_rad": 0.1,
                "first_posterior_covariance_perturbation_upper": 0.2,
                "sample1_reduced_covariance_PSD_perturbation_upper": 0.3,
                "first_PSD_Joseph_tangent_column_norm_upper": 0.4,
            }

        def fake_build(*args, **kwargs):
            seen["patched"] = V41.V12D._first_psd_perturbation_tangent is not original_hook
            self.assertEqual(
                V41.V12D._first_psd_perturbation_tangent(dummy=1)[
                    "first_offaxis_attitude_correction_upper_rad"],
                0.1,
            )
            return {
                "source_generated_not_trajectory_fit": True,
                "source_replay_used": False,
                "filter_changed": False,
                "P5_SAMPLE1_SIGNED_FULL_ANGLE_CURRENT_YZ_Q8_V18B": "NOT_ESTABLISHED",
                "evaluated_signed_cayley_cells": 10,
                "unclosed_q8_cells": 2,
                "first_unclosed_q8_cell": {"p_cell": 0},
            }

        V41.V40._first_psd_perturbation_exact_joseph_components = fake_refined
        V41.V18B.build = fake_build
        V41.V18B.validate = lambda _d: []
        try:
            d = V41.build(V41.DEFAULT_DOMAIN)
        finally:
            V41.V40._first_psd_perturbation_exact_joseph_components = original_refined
            V41.V18B.build = original_build
            V41.V18B.validate = original_validate

        self.assertTrue(seen.get("patched"))
        self.assertEqual(seen.get("kwargs"), {"dummy": 1})
        self.assertIs(V41.V12D._first_psd_perturbation_tangent, original_hook)
        self.assertEqual(d["V40_exact_Joseph_first_PSD_helper_calls"], 1)
        self.assertEqual(d["P5_SAMPLE1_V40_FULL_SOURCE_CELL0_Q8_LIFT_V41"], "NOT_ESTABLISHED")
        self.assertFalse(d["q8_word_promoted_here"])
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])

    def test_validation_keeps_scope_and_promotion_guards(self):
        d = {
            "schema": V41.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_V40_FULL_SOURCE_CELL0_Q8_LIFT_V41",
            "failures": [],
            "source_generated_not_trajectory_fit": True,
            "V18B_full_signed_angle_current_yz_parent_retained": True,
            "V40_exact_Joseph_first_PSD_installed_globally": True,
            "refined_V18_current_chart_recomputed_from_V40_bound": True,
            "V18_support_intersections_with_parent_retained": True,
            "temporary_V12D_helper_restored_after_build": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "V40_exact_Joseph_first_PSD_helper_calls": 1,
            "deployed_correction_limit_rad": 6.0,
            "q_target": V41.Q_TARGET,
            "first_unclosed_q8_cell": {"p_cell": 0},
            "P5_SAMPLE1_V40_FULL_SOURCE_CELL0_Q8_LIFT_V41": "NOT_ESTABLISHED",
        }
        self.assertEqual(V41.validate(d), [])
        d["whole_word_promoted_here"] = True
        self.assertIn("whole_word_promoted_here is not false", V41.validate(d))


if __name__ == "__main__":
    unittest.main()
