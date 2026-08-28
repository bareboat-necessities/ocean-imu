#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_v36_psd_cone_current_subbox_lift_v37 as V37


class Sample1V36PSDConeCurrentSubboxLiftV37Tests(unittest.TestCase):
    def test_build_installs_and_restores_v36_psd_cone(self):
        original_build = V37.V34.build
        original_validate = V37.V34.validate
        original_cone = V37.V36._first_psd_perturbation_psd_cone
        original_helper = V37.V12D._first_psd_perturbation_tangent
        seen = {}

        def fake_cone(**kwargs):
            seen["cone_kwargs"] = dict(kwargs)
            return {"fake_psd_cone": True}

        def fake_build(*args, **kwargs):
            seen["helper_is_patched"] = (
                V37.V12D._first_psd_perturbation_tangent is not original_helper
            )
            got = V37.V12D._first_psd_perturbation_tangent(dummy=1)
            self.assertEqual(got, {"fake_psd_cone": True})
            return {
                "P5_SAMPLE1_DIRECTIONAL_INNOVATION_ROW_LIFT_V34": "PASS",
                "open_current_subboxes": 1,
            }

        V37.V36._first_psd_perturbation_psd_cone = fake_cone
        V37.V34.build = fake_build
        V37.V34.validate = lambda _d: []
        try:
            d = V37.build(V37.DEFAULT_DOMAIN)
        finally:
            V37.V34.build = original_build
            V37.V34.validate = original_validate
            V37.V36._first_psd_perturbation_psd_cone = original_cone

        self.assertTrue(seen.get("helper_is_patched"))
        self.assertEqual(seen.get("cone_kwargs"), {"dummy": 1})
        self.assertIs(V37.V12D._first_psd_perturbation_tangent, original_helper)
        self.assertEqual(d["V36_PSD_cone_helper_calls"], 1)
        self.assertEqual(
            d["next_obligation"],
            "REFINE_FIRST_REMAINING_V37_SUBBOX_WITH_EXACT_FIRST_PSD_RESET_COMPONENT_MATRIX",
        )

    def test_validation_keeps_promotion_guards(self):
        d = {
            "schema": V37.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_V36_PSD_CONE_CURRENT_SUBBOX_LIFT_V37",
            "failures": [],
            "source_generated_not_trajectory_fit": True,
            "V34_directional_innovation_construction_retained": True,
            "V36_PSD_offdiagonal_cone_installed": True,
            "V36_changes_only_first_PSD_offdiagonal_operator_bound": True,
            "all_candidate_current_subboxes_accounted": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "V36_PSD_cone_helper_calls": 1,
            "deployed_correction_limit_rad": 6.0,
            "q_target": V37.Q_TARGET,
            "P5_SAMPLE1_V36_PSD_CONE_CURRENT_SUBBOX_LIFT_V37": "PASS",
        }
        self.assertEqual(V37.validate(d), [])
        d["N_H_words_set_here"] = True
        self.assertIn("N_H_words_set_here is not false", V37.validate(d))


if __name__ == "__main__":
    unittest.main()
