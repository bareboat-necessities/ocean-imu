from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_nextgen_gain_certificate as P4G
import ou3_validate_enclosure as ENC


class Ou3P4NextgenGainCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = P4G.build()

    def test_gain_refined_certificate_passes(self):
        self.assertEqual(P4G.validate(self.d), [])
        self.assertEqual(self.d["P4_NEXTGEN_GAIN_WORD_CERTIFICATE"], "PASS")
        self.assertTrue(self.d["nextgen_measurement_specific_gain_refinement"])

    def test_class_local_gain_bounds_are_monotone(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            self.assertTrue(m["measurement_specific_gain_bounds_monotone"])
            g = m["measurement_specific_gain_norm_upper"]
            for value in g.values():
                self.assertLessEqual(value, m["global_gain_norm_upper_previous"])
            self.assertTrue(m["gain_refined_defect_sum_monotone"])
            self.assertLessEqual(
                m["transported_word_defect_B_upper"],
                m["transported_word_defect_B_upper_previous_gain_stage"],
            )

    def test_W_never_regresses_and_total_factor_is_at_least_one(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            self.assertGreaterEqual(m["certified_level_W"], m["certified_level_W_previous_gain_stage"])
            self.assertGreaterEqual(m["gain_stage_W_widening_factor_lower"], 1.0)
            self.assertGreaterEqual(m["total_W_widening_factor_vs_legacy_lower"], 1.0)

    def test_schema4_enclosure_compatibility_and_safety(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            out = ENC.validate_mode(
                mode, m, {"required_path_metric": "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC"}
            )
            self.assertTrue(out["linear_pass"], out["failures"])
            self.assertTrue(out["nonlinear_pass"], out["failures"])
            self.assertLess(m["accepted_correction_norm_prefix_upper"], 1.0e-2)
        p = self.d["modes"]["A"]["active_bias_projection"]
        self.assertFalse(p["projection_surface_reached_in_certified_funnel"])

    def test_no_replay_sampling_or_numpy(self):
        text = (ROOT / "tools" / "ou3_p4_nextgen_gain_certificate.py").read_text(
            encoding="utf-8"
        ).casefold()
        self.assertNotIn("monte carlo", text)
        self.assertNotIn("sampled trajectory", text)
        self.assertNotIn("numpy", text)


if __name__ == "__main__":
    unittest.main()
