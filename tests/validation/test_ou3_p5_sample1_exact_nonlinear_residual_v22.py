from __future__ import annotations
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_exact_nonlinear_residual_v22 as V22


class Sample1ExactNonlinearResidualV22Tests(unittest.TestCase):
    def test_exact_rotation_residual_is_zero_at_zero_cayley(self):
        z = Interval.point(0.0)
        f = [z, Interval.point(-1.0), Interval.point(9.0)]
        r = V22.exact_rotation_residual([z, z, z], f)
        for x in r:
            self.assertLessEqual(x.lo, 0.0)
            self.assertGreaterEqual(x.hi, 0.0)

    def test_exact_rotation_residual_matches_shipping_small_angle_sign(self):
        z = Interval.point(0.0)
        c = [Interval.point(0.1), z, z]
        f = [z, Interval.point(0.0), Interval.point(9.0)]
        r = V22.exact_rotation_residual(c, f)
        # H_theta c = c x f has negative y for positive c_x and positive f_z.
        self.assertLess(r[1].hi, 0.0)
        self.assertLessEqual(r[2].hi, 0.0)

    def test_validation_keeps_exact_algebra_and_promotion_guards(self):
        d = {
            "schema": V22.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_EXACT_NONLINEAR_RESIDUAL_WITNESS_V22",
            "source_generated_not_trajectory_fit": True,
            "V21B_authoritative_current_chart_retained": True,
            "exact_Cayley_rotation_residual_used": True,
            "linear_attitude_plus_defect_double_charge_retired": True,
            "linear_aw_plus_latent_cross_double_charge_retired": True,
            "rotated_physical_aw_norm_preserved_by_orthogonality": True,
            "accelerometer_bias_norm_retained": True,
            "V12D_PSD_S_correction_perturbation_retained": True,
            "V16_axis_cone_and_V18_yz_support_used_for_product_check": True,
            "source_replay_used": False,
            "filter_changed": False,
            "deployed_correction_limit_increased": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "deployed_correction_limit_rad": 6.0,
            "q_target": 8.0,
            "combined_rotated_aw_plus_bias_nuisance_norm_upper_mps2": 7.0,
            "V21B_previous_nuisance_norm_upper_mps2": 14.0,
            "source_subcell_incompatible_under_exact_residual": False,
            "baseline_correction_radial_upper_rad": 2.08,
            "refined_correction_radial_upper_rad": 2.0,
            "P5_SAMPLE1_EXACT_NONLINEAR_RESIDUAL_WITNESS_V22": "PASS",
            "failures": [],
        }
        self.assertEqual(V22.validate(d), [])


if __name__ == "__main__":
    unittest.main()
