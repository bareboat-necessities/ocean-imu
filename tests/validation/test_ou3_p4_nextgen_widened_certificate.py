from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_nextgen_widened_certificate as P4W
import ou3_validate_enclosure as ENC


class Ou3P4NextgenWidenedCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = P4W.build()

    def test_nextgen_H_A_certificate_passes_and_strictly_widens(self):
        d = self.d
        self.assertEqual(P4W.validate(d), [])
        self.assertEqual(d["P4_NEXTGEN_WIDENED_WORD_CERTIFICATE"], "PASS")
        self.assertTrue(d["nextgen_refinement_source_only"])
        self.assertTrue(d["legacy_P4_retained_as_baseline"])
        self.assertFalse(d["source_replay_used"])
        for mode in ("H", "A"):
            m = d["modes"][mode]
            self.assertTrue(m["nextgen_operation_specific_defect_transport"])
            self.assertTrue(m["operation_specific_sum_no_larger_than_legacy_max_budget"])
            self.assertLess(m["transported_word_defect_B_upper"],
                            m["transported_word_defect_B_upper_legacy"])
            self.assertGreater(m["certified_level_W"], m["certified_level_W_legacy"])
            self.assertGreater(m["certified_level_W_widening_factor_lower"], 1.0)
            self.assertGreater(m["certified_level_sqrt_W_widening_factor_lower"], 1.0)
            self.assertLess(m["prefix_canonical_error_norm_upper"], m["cayley_norm_limit"])
            self.assertLess(m["accepted_correction_norm_prefix_upper"], 1.0e-2)

    def test_operation_class_sum_is_the_bound_used_for_transport(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            c = m["operation_specific_quadratic_defect_constants_upper"]
            self.assertEqual(set(c), {
                "prediction", "S_zero_accepted", "accelerometer_accepted", "magnetometer_accepted"
            })
            self.assertGreater(c["S_zero_accepted"], 0.0)
            self.assertGreater(c["accelerometer_accepted"], 0.0)
            self.assertGreater(c["magnetometer_accepted"], 0.0)
            self.assertGreaterEqual(m["operation_specific_defect_sum_per_sample_upper"], sum(c.values()))
            self.assertLess(m["operation_specific_defect_sum_per_sample_upper"],
                            m["legacy_four_operation_max_defect_per_sample_upper"])

    def test_widened_rows_still_satisfy_schema4_enclosure_validator(self):
        for mode in ("H", "A"):
            out = ENC.validate_mode(
                mode,
                self.d["modes"][mode],
                {"required_path_metric": "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC"},
            )
            self.assertTrue(out["linear_pass"], out["failures"])
            self.assertTrue(out["nonlinear_pass"], out["failures"])
            self.assertGreater(out["mu_W_lower"], 0.0)

    def test_A_widened_funnel_remains_inside_bias_projection_interior(self):
        p = self.d["modes"]["A"]["active_bias_projection"]
        self.assertFalse(p["projection_surface_reached_in_certified_funnel"])
        self.assertLess(p["certified_error_norm_prefix_upper"], p["interior_margin_lower_mps2"])

    def test_refinement_does_not_use_replay_sampling_or_numpy(self):
        text = (ROOT / "tools" / "ou3_p4_nextgen_widened_certificate.py").read_text(
            encoding="utf-8"
        ).casefold()
        self.assertNotIn("numpy", text)
        self.assertNotIn("monte carlo", text)
        self.assertNotIn("sampled trajectory", text)
        self.assertIn("operation-specific", text)


if __name__ == "__main__":
    unittest.main()
