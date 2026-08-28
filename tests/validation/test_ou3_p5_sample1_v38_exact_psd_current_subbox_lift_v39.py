#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_v38_exact_psd_current_subbox_lift_v39 as V39


class Sample1V38ExactPSDCurrentSubboxLiftV39Tests(unittest.TestCase):
    def test_build_isolates_refined_psd_from_authoritative_current_chart(self):
        original_v12_build = V39.V12D.build
        original_v12_validate = V39.V12D.validate
        original_v34_build = V39.V34.build
        original_v34_validate = V39.V34.validate
        original_exact = V39.V38._first_psd_perturbation_exact_correction
        original_helper = V39.V12D._first_psd_perturbation_tangent
        original_chart = V39.V21._current_component_chart
        original_witness = V39.V30._witness_row
        original_matches = V39.V21B._matches_reference
        seen = {}

        baseline_row = {"first_offaxis_attitude_correction_upper_rad": 0.25}

        V39.V12D.build = lambda *a, **k: {"rows": [baseline_row]}
        V39.V12D.validate = lambda _d: []
        V39.V30._witness_row = lambda _d: baseline_row
        V39.V21B._matches_reference = lambda q: abs(float(q) - 0.6415230535178351) < 1e-12

        def fake_exact(**kwargs):
            seen["exact_kwargs"] = dict(kwargs)
            return {"fake_exact": True}

        def fake_chart(*, first, base, vr, dom, src, sample1_s_angle):
            seen["chart_vr"] = vr
            return {"q1": 0.6415230535178351}

        def fake_v34_build(*args, **kwargs):
            seen["helper_patched"] = V39.V12D._first_psd_perturbation_tangent is not original_helper
            seen["chart_patched"] = V39.V21._current_component_chart is not original_chart
            self.assertEqual(
                V39.V12D._first_psd_perturbation_tangent(dummy=3),
                {"fake_exact": True},
            )
            chart = V39.V21._current_component_chart(
                first={}, base={}, vr={"refined": True}, dom={}, src={},
                sample1_s_angle=0.0)
            self.assertEqual(chart["q1"], 0.6415230535178351)
            return {
                "P5_SAMPLE1_DIRECTIONAL_INNOVATION_ROW_LIFT_V34": "PASS",
                "sample1_current_cayley_norm_upper": 0.6415230535178351,
                "open_current_subboxes": 1,
            }

        V39.V38._first_psd_perturbation_exact_correction = fake_exact
        V39.V21._current_component_chart = fake_chart
        V39.V34.build = fake_v34_build
        V39.V34.validate = lambda _d: []
        try:
            d = V39.build(V39.DEFAULT_DOMAIN)
        finally:
            V39.V12D.build = original_v12_build
            V39.V12D.validate = original_v12_validate
            V39.V34.build = original_v34_build
            V39.V34.validate = original_v34_validate
            V39.V38._first_psd_perturbation_exact_correction = original_exact
            V39.V21._current_component_chart = original_chart
            V39.V30._witness_row = original_witness
            V39.V21B._matches_reference = original_matches

        self.assertTrue(seen.get("helper_patched"))
        self.assertTrue(seen.get("chart_patched"))
        self.assertEqual(seen.get("chart_vr"), baseline_row)
        self.assertEqual(seen.get("exact_kwargs"), {"dummy": 3})
        self.assertIs(V39.V12D._first_psd_perturbation_tangent, original_helper)
        self.assertIs(V39.V21._current_component_chart, original_chart)
        self.assertEqual(d["V38_exact_first_PSD_helper_calls"], 1)
        self.assertEqual(d["authoritative_V18B_current_chart_freeze_calls"], 1)
        self.assertEqual(
            d["next_obligation"],
            "REFINE_FIRST_REMAINING_V39_SUBBOX_WITH_EXACT_FIRST_PSD_JOSEPH_COMPONENT_MATRIX",
        )

    def test_validation_keeps_authoritative_q_and_promotion_guards(self):
        q = V39.V21B.V18B_FIRST_WITNESS_CURRENT_Q
        d = {
            "schema": V39.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_V38_ISOLATED_EXACT_PSD_CURRENT_SUBBOX_LIFT_V39",
            "failures": [],
            "source_generated_not_trajectory_fit": True,
            "V34_directional_innovation_construction_retained": True,
            "V38_exact_first_PSD_correction_installed": True,
            "V38_exact_canonical_tangent_geometry_retained": True,
            "V36_full_Joseph_gain_operator_parent_retained": True,
            "authoritative_V18B_current_chart_frozen_to_baseline_V12D_witness": True,
            "current_q_matches_authoritative_V18B_reference": True,
            "refined_PSD_used_only_outside_authoritative_current_chart": True,
            "temporary_helpers_restored_after_build": True,
            "all_candidate_current_subboxes_accounted": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "V38_exact_first_PSD_helper_calls": 1,
            "authoritative_V18B_current_chart_freeze_calls": 1,
            "deployed_correction_limit_rad": 6.0,
            "q_target": V39.Q_TARGET,
            "sample1_current_cayley_norm_upper": q,
            "P5_SAMPLE1_V38_ISOLATED_EXACT_PSD_CURRENT_SUBBOX_LIFT_V39": "PASS",
        }
        self.assertEqual(V39.validate(d), [])
        d["N_H_words_set_here"] = True
        self.assertIn("N_H_words_set_here is not false", V39.validate(d))


if __name__ == "__main__":
    unittest.main()
