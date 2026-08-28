#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_first_psd_exact_correction_geometry_v38 as V38


class Sample1FirstPSDExactCorrectionGeometryV38Tests(unittest.TestCase):
    def test_exact_canonical_correction_is_nonworsening_on_synthetic_cell(self):
        kwargs = dict(
            t=Interval.outward_bounds(0.01, 0.0101),
            Y=Interval.outward_bounds(0.02, 0.0201),
            p=Interval.outward_bounds(0.001, 0.0011),
            r=Interval.outward_bounds(0.0005, 0.0005),
            g=9.80665,
            eps=1.0e-6,
            rho0=1.0,
            dhi=0.1,
            rt=Interval.outward_bounds(0.4, 0.5),
            rz=Interval.outward_bounds(-0.2, 0.2),
            alpha_hi=0.99,
            aw_pre=0.1,
        )
        parent = V38.V36._first_psd_perturbation_psd_cone(**kwargs)
        refined = V38._first_psd_perturbation_exact_correction(**kwargs)
        for key in (
            "first_offaxis_attitude_correction_upper_rad",
            "first_aw_x_correction_upper_mps2",
            "PSD_induced_sample1_residual_perturbation_upper_mps2",
            "reset_gauge_transform_perturbation_upper",
            "sample1_reduced_covariance_PSD_perturbation_upper",
        ):
            self.assertLessEqual(float(refined[key]), V38.FULL.up(float(parent[key])))
        self.assertTrue(refined["first_PSD_exact_canonical_residual_direction_used"])
        self.assertTrue(refined["first_PSD_axial_residual_has_zero_attitude_effect"])
        self.assertTrue(refined["V36_full_gain_operator_retained_for_Joseph_covariance"])

    def test_validation_keeps_promotion_guards(self):
        d = {
            "schema": V38.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_FIRST_PSD_EXACT_CORRECTION_GEOMETRY_V38",
            "failures": [],
            "source_generated_not_trajectory_fit": True,
            "V36_PSD_cone_parent_revalidated": True,
            "canonical_first_residual_direction_used": True,
            "exact_2x2_tangent_inverse_geometry_used": True,
            "axial_residual_zero_PSD_attitude_effect_used": True,
            "V36_Joseph_gain_operator_parent_retained": True,
            "source_replay_used": False,
            "filter_changed": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "P5_SAMPLE1_FIRST_PSD_EXACT_CORRECTION_GEOMETRY_V38": "PASS",
        }
        self.assertEqual(V38.validate(d), [])
        d["whole_word_promoted_here"] = True
        self.assertIn("whole_word_promoted_here is not false", V38.validate(d))


if __name__ == "__main__":
    unittest.main()
