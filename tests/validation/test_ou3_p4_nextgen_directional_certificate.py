from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_nextgen_directional_certificate as P4D
import ou3_validate_enclosure as ENC


class Ou3P4NextgenDirectionalCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = P4D.build()

    def test_directional_exact_budget_certificate_passes(self):
        d = self.d
        self.assertEqual(P4D.validate(d), [])
        self.assertEqual(d["P4_NEXTGEN_DIRECTIONAL_WORD_CERTIFICATE"], "PASS")
        self.assertTrue(d["nextgen_directional_operator_refinement"])
        self.assertTrue(d["nextgen_exact_endpoint_budget_refinement"])

    def test_each_stage_is_monotone_and_total_widening_exceeds_stage1(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            self.assertTrue(m["directional_defect_sum_monotone"])
            self.assertLessEqual(
                m["transported_word_defect_B_upper"],
                m["transported_word_defect_B_upper_previous"],
            )
            self.assertGreaterEqual(
                m["certified_level_W"], m["certified_level_W_previous_nextgen"]
            )
            self.assertGreaterEqual(m["secondgen_W_widening_factor_lower"], 1.0)
            self.assertGreaterEqual(
                m["total_W_widening_factor_vs_legacy_lower"],
                m["certified_level_W_widening_factor_lower"],
            )

    def test_measurement_specific_H_bounds_refine_global_bound(self):
        for mode in ("H", "A"):
            h = self.d["modes"][mode]["directional_measurement_operator_norm_upper"]
            self.assertLessEqual(h["S_zero"], h["global_previous"])
            self.assertLessEqual(h["accelerometer"], h["global_previous"])
            self.assertLessEqual(h["magnetometer"], h["global_previous"])
            self.assertLess(h["S_zero"], h["magnetometer"])
            self.assertLess(h["accelerometer"], h["magnetometer"])

    def test_exact_endpoint_and_prefix_budgets_close(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            self.assertTrue(m["legacy_delta_over_8_budget_replaced"])
            self.assertGreater(m["exact_endpoint_sqrt_gap_lower"], 0.0)
            self.assertLessEqual(m["prefix_bootstrap_B_sqrt_W_upper"], 1.0)
            self.assertLess(m["prefix_canonical_error_norm_upper"], m["cayley_norm_limit"])
            self.assertLess(m["accepted_correction_norm_prefix_upper"], 1.0e-2)

    def test_rows_remain_compatible_with_enclosure_validator(self):
        for mode in ("H", "A"):
            out = ENC.validate_mode(
                mode,
                self.d["modes"][mode],
                {"required_path_metric": "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC"},
            )
            self.assertTrue(out["linear_pass"], out["failures"])
            self.assertTrue(out["nonlinear_pass"], out["failures"])

    def test_A_projection_remains_interior(self):
        p = self.d["modes"]["A"]["active_bias_projection"]
        self.assertFalse(p["projection_surface_reached_in_certified_funnel"])
        self.assertLess(p["certified_error_norm_prefix_upper"], p["interior_margin_lower_mps2"])

    def test_no_replay_or_sampling_route(self):
        text = (ROOT / "tools" / "ou3_p4_nextgen_directional_certificate.py").read_text(
            encoding="utf-8"
        ).casefold()
        self.assertNotIn("monte carlo", text)
        self.assertNotIn("sampled trajectory", text)
        self.assertNotIn("numpy", text)


if __name__ == "__main__":
    unittest.main()
