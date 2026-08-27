from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_effective_input_correction_v21b as V21B


class Sample1EffectiveInputCorrectionV21BTests(unittest.TestCase):
    def test_v14d_binding_is_scoped_and_restored(self):
        original = V21B.V21.V14._normalized_shipping_quaternion
        seen = V21B._with_v14d_quaternion(
            lambda: V21B.V21.V14._normalized_shipping_quaternion)
        self.assertIs(seen, V21B.V14D.radial_sinc_normalized_shipping_quaternion)
        self.assertIs(V21B.V21.V14._normalized_shipping_quaternion, original)

    def test_authoritative_v18b_reference_is_fixed(self):
        self.assertTrue(V21B._matches_reference(V21B.V18B_FIRST_WITNESS_CURRENT_Q))
        self.assertFalse(V21B._matches_reference(0.6593778441001633))
        self.assertEqual(V21B.Q_TARGET, 8.0)

    def test_validation_keeps_parent_and_promotion_guards(self):
        d = {
            "schema": V21B.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_V14D_BOUND_EFFECTIVE_INPUT_WITNESS_V21B",
            "source_generated_not_trajectory_fit": True,
            "V21_effective_input_parent_retained": True,
            "V14D_radial_sinc_quaternion_installed_for_V21_audit": True,
            "V14D_quaternion_restored_after_audit": True,
            "current_q_matches_authoritative_V18B_reference": True,
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
            "sample1_current_cayley_norm_upper": V21B.V18B_FIRST_WITNESS_CURRENT_Q,
            "P5_SAMPLE1_V14D_BOUND_EFFECTIVE_INPUT_WITNESS_V21B": "PASS",
            "failures": [],
        }
        self.assertEqual(V21B.validate(d), [])


if __name__ == "__main__":
    unittest.main()
