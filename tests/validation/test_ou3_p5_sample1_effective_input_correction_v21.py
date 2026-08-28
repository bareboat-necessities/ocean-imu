from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_effective_input_correction_v21 as V21


class Sample1EffectiveInputCorrectionV21Tests(unittest.TestCase):
    def test_box_intersection_tightens_and_detects_empty(self):
        a = [Interval.outward_bounds(-2.0, 2.0), Interval.outward_bounds(-1.0, 1.0)]
        b = [Interval.outward_bounds(0.5, 3.0), Interval.outward_bounds(-0.25, 0.25)]
        c = V21._intersect_boxes(a, b)
        self.assertIsNotNone(c)
        self.assertEqual(c[0].lo, max(a[0].lo, b[0].lo))
        self.assertEqual(c[0].hi, min(a[0].hi, b[0].hi))
        self.assertEqual(c[1].lo, max(a[1].lo, b[1].lo))
        self.assertEqual(c[1].hi, min(a[1].hi, b[1].hi))
        self.assertIsNone(V21._intersect_boxes(
            [Interval.outward_bounds(-2.0, -1.0)],
            [Interval.outward_bounds(1.0, 2.0)]))

    def test_q_target_and_witness_indices_are_fixed(self):
        self.assertEqual(V21.Q_TARGET, 8.0)
        self.assertEqual(V21.WITNESS_P_CELL, 0)
        self.assertEqual(V21.WITNESS_TANGENT_CELL, 0)
        self.assertEqual(V21.WITNESS_AXIAL_CELL, 19)
        self.assertEqual(V21.WITNESS_RX_CELL, 0)
        self.assertEqual(V21.WITNESS_PARALLEL_CELL, 0)

    def test_validation_keeps_estimator_and_promotion_guards(self):
        d = {
            "schema": V21.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_EFFECTIVE_INPUT_CORRECTION_WITNESS_V21",
            "source_generated_not_trajectory_fit": True,
            "V12D_PSD_S_perturbation_retained": True,
            "V10_one_plus_two_gain_retained": True,
            "exact_accelerometer_effective_input_lemma_used": True,
            "current_cayley_and_sample1_correction_jointly_mapped": True,
            "V13E_signed_subcell_intersected_not_replaced": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "deployed_correction_limit_rad": 6.0,
            "q_target": 8.0,
            "sample1_current_cayley_norm_upper": 0.64,
            "nominal_effective_residual_nuisance_norm_upper_mps2": 1.0,
            "V12D_correction_perturbation_norm_upper_rad": 0.1,
            "baseline_correction_radial_upper_rad": 2.08,
            "refined_correction_radial_upper_rad": 2.0,
            "source_subcell_incompatible": False,
            "joint_correction_box_rad": [[-1.0, -0.5], [-0.2, 0.2], [-0.2, 0.2]],
            "P5_SAMPLE1_EFFECTIVE_INPUT_CORRECTION_WITNESS_V21": "PASS",
            "failures": [],
        }
        self.assertEqual(V21.validate(d), [])


if __name__ == "__main__":
    unittest.main()
