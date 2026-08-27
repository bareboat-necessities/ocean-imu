from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_split_gravity_signed_components_v28 as V28


class Sample1SplitGravitySignedComponentsV28Tests(unittest.TestCase):
    def test_axial_decay_is_tighter_for_positive_cosine(self):
        d = V28._gravity_component_decay_bounds(
            cosine_lower=0.99, alpha_lower=0.9, gravity=9.80665)
        self.assertGreater(d["tangent_gravity_decay_remainder_upper_mps2"], 0.0)
        self.assertGreaterEqual(d["tangent_gravity_decay_remainder_upper_mps2"],
                                d["axial_gravity_decay_remainder_upper_mps2"])

    def test_axial_decay_is_only_outward_dust_at_exact_alignment(self):
        d = V28._gravity_component_decay_bounds(
            cosine_lower=1.0, alpha_lower=0.9, gravity=9.80665)
        self.assertGreaterEqual(d["gravity_axial_residual_upper_mps2"], 0.0)
        self.assertLess(d["gravity_axial_residual_upper_mps2"], 1e-320)
        self.assertGreaterEqual(d["axial_gravity_decay_remainder_upper_mps2"], 0.0)
        self.assertLess(d["axial_gravity_decay_remainder_upper_mps2"], 1e-320)

    def test_validation_keeps_split_and_promotion_guards(self):
        d = {
            "schema": V28.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_SPLIT_GRAVITY_SIGNED_COMPONENTS_V28",
            "source_generated_not_trajectory_fit": True,
            "V27_signed_post_first_aw_parent_retained": True,
            "first_accel_source_certificate_revalidated": True,
            "post_prediction_gravity_cosine_source_bound_used": True,
            "tangent_and_axial_gravity_decay_split": True,
            "shared_tangent_gravity_decay_for_axial_channel_retired": True,
            "V23_first_open_subbox_retained": True,
            "V10_exact_first_update_OU_cancellation_used": True,
            "signed_tangent_axial_first_residual_cell_retained": True,
            "V21_signed_one_plus_two_gain_components_used": True,
            "V12D_correction_perturbation_retained_as_single_ball": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "deployed_correction_limit_rad": 6.0,
            "q_target": 8.0,
            "split_gravity_detail": {
                "post_prediction_true_gravity_cosine_lower": 0.9,
                "tangent_gravity_decay_remainder_upper_mps2": 0.2,
                "axial_gravity_decay_remainder_upper_mps2": 0.1,
            },
            "P5_SAMPLE1_SPLIT_GRAVITY_SIGNED_COMPONENTS_V28": "PASS",
            "failures": [],
        }
        self.assertEqual(V28.validate(d), [])


if __name__ == "__main__":
    unittest.main()
