import importlib.util
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location(
    "ou3_vector_uco_certificate", ROOT / "tools" / "ou3_vector_uco_certificate.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class VectorUcoCertificateTests(unittest.TestCase):
    def test_conditional_vector_uco_is_strict_and_explicit(self):
        d = mod.build()
        self.assertEqual(mod.validate(d), [])
        self.assertTrue(d["pass"])
        self.assertTrue(d["conditional_source_complete"])
        self.assertFalse(d["unconditional_full_heading_claim"])
        self.assertTrue(d["persistent_excitation_is_theorem_hypothesis"])
        self.assertGreater(d["vector_pair"]["mu_theta_lower"], 0.0)
        self.assertGreater(
            d["gyro_bias_two_packet"]["alpha_6_information_lower"], 0.0
        )

    def test_rate_gap_condition_is_strict(self):
        d = mod.build()
        g = d["gyro_bias_two_packet"]
        self.assertLess(g["omega_times_gap_upper"], 2.0)
        self.assertGreater(g["gamma_bracket_lower"], 0.0)
        self.assertGreater(g["Gamma_g_sigma_min_lower_s"], 0.0)

    def test_measurement_bounds_come_from_configured_source(self):
        d = mod.build()
        c = d["configured_measurement_bounds"]
        self.assertAlmostEqual(c["mag_odr_hz"], 25.0)
        self.assertGreater(c["acc_measurement_variance_upper"], 0.0)
        self.assertGreater(c["mag_measurement_variance_upper"], 0.0)
        # Never use a magnetic norm floor stronger than the source acquisition guard.
        self.assertLessEqual(
            d["operating_envelope"]["magnetic_vector_norm_lower_uT"],
            c["mag_init_norm_guard_uT"],
        )

    def test_no_replay_or_observed_extrema_are_used(self):
        text = (ROOT / "tools" / "ou3_vector_uco_certificate.py").read_text()
        for forbidden in (
            "ou3_exact_replay", "path_metrics", "neighborhood_radius_search",
            "np.quantile", "observed_min", "replay_min",
        ):
            self.assertNotIn(forbidden, text)

    def test_collinearity_is_not_silently_certified(self):
        self.assertGreater(mod.PE["vector_sine_separation_lower"], 0.0)
        d = mod.build()
        self.assertIn("collinearity", d["theorem_scope"])


if __name__ == "__main__":
    unittest.main()
