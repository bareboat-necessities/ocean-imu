#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_p5_sample1_first_psd_exact_joseph_components_v40 as V40


class Sample1FirstPSDExactJosephComponentsV40Tests(unittest.TestCase):
    def test_attitude_supported_transport_is_nonworsening_on_synthetic_cell(self):
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
        parent = V40.V38._first_psd_perturbation_exact_correction(**kwargs)
        refined = V40._first_psd_perturbation_exact_joseph_components(**kwargs)

        for key in (
            "first_offaxis_attitude_correction_upper_rad",
            "first_aw_x_correction_upper_mps2",
            "PSD_induced_sample1_residual_perturbation_upper_mps2",
            "first_posterior_covariance_perturbation_upper",
            "reset_gauge_transform_perturbation_upper",
            "sample1_reduced_covariance_PSD_perturbation_upper",
        ):
            self.assertLessEqual(float(refined[key]), V40.FULL.up(float(parent[key])))

        self.assertLess(
            float(refined["first_posterior_covariance_perturbation_upper"]),
            float(parent["first_posterior_covariance_perturbation_upper"]),
        )
        self.assertTrue(refined["first_PSD_Joseph_attitude_support_used"])
        self.assertTrue(refined["first_PSD_Joseph_nominal_columns_orthogonal"])
        self.assertTrue(refined["first_PSD_Joseph_zero_diagonal_component_matrix_used"])
        self.assertTrue(refined["V38_exact_mean_correction_geometry_retained"])

    def test_validation_keeps_promotion_guards(self):
        d = {
            "schema": V40.SCHEMA,
            "qualification": "OU3_P5_SAMPLE1_FIRST_PSD_EXACT_JOSEPH_COMPONENTS_V40",
            "failures": [],
            "source_generated_not_trajectory_fit": True,
            "V38_exact_correction_parent_revalidated": True,
            "attitude_supported_Joseph_columns_used": True,
            "zero_diagonal_PSD_component_matrix_used": True,
            "deltaK_cross_terms_retained": True,
            "source_replay_used": False,
            "filter_changed": False,
            "q8_composed_here": False,
            "q8_word_promoted_here": False,
            "whole_word_promoted_here": False,
            "N_H_words_set_here": False,
            "P5_SAMPLE1_FIRST_PSD_EXACT_JOSEPH_COMPONENTS_V40": "PASS",
        }
        self.assertEqual(V40.validate(d), [])
        d["q8_word_promoted_here"] = True
        self.assertIn("q8_word_promoted_here is not false", V40.validate(d))


if __name__ == "__main__":
    unittest.main()
