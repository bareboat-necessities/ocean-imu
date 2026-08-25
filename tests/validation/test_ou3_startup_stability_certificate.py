import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_implementation_proof_manifest as MANIFEST
import ou3_startup_stability_certificate as STARTUP


class Ou3ImplementationManifestTests(unittest.TestCase):
    def test_manifest_binds_deployed_startup_and_hybrid_semantics(self):
        d = MANIFEST.build()
        self.assertEqual(MANIFEST.validate(d), [])
        self.assertEqual(d["state_coordinates"]["H_dimension"], 18)
        self.assertEqual(d["state_coordinates"]["A_dimension"], 21)
        self.assertAlmostEqual(d["startup"]["two_kp"], 0.2, places=7)
        self.assertAlmostEqual(d["startup"]["two_ki"], 0.02, places=8)
        self.assertEqual(d["startup"]["gravity_align_max_sin"], 0.075)
        self.assertEqual(d["startup"]["proxy_startup_timeout_sec"], 150.0)
        self.assertTrue(d["startup"]["timeout_cannot_handoff_antipodal_branch"])
        self.assertTrue(d["startup"]["go_live_bias_learning_held"])
        self.assertIn("held_to_active", d["hybrid_events"])
        self.assertIn("magnetic_regauge_refinement", d["hybrid_events"])
        self.assertIn("tilt_relock", d["hybrid_events"])


class Ou3StartupStabilityCertificateTests(unittest.TestCase):
    def test_declared_operating_domain_is_not_trajectory_fitted(self):
        import json
        d = json.loads(STARTUP.DEFAULT_DOMAIN.read_text())
        self.assertFalse(d["trajectory_fit"])
        self.assertEqual(d["normal_live"]["vector_pe_recurrence_window_s"], 1.0)

    def test_startup_reset_and_handoff_have_numeric_margins(self):
        d = STARTUP.build()
        self.assertEqual(STARTUP.validate(d), [])
        self.assertTrue(d["startup_certificate_pass"])
        self.assertGreater(d["source_global_reset"]["post_reset_true_gravity_cosine_lower"], 0.0)
        self.assertGreater(d["mahony"]["chart_invariance_margin_lower"], 0.0)
        self.assertGreater(d["normal_handoff"]["true_gravity_cosine_lower"], 0.0)
        self.assertGreater(d["timeout_handoff"]["combined_true_gravity_cosine_lower"], 0.0)
        self.assertEqual(d["go_live"]["first_mode"], "H")
        self.assertTrue(d["go_live"]["bias_learning_held"])

    def test_sqrt_enclosure_is_verified_by_exact_square_comparison(self):
        from fractions import Fraction
        for x in (0.25, 0.5, 2.0, 9.80665):
            r = STARTUP.sqrt_interval_point(x)
            q = Fraction.from_float(x)
            self.assertLessEqual(Fraction.from_float(r.lo) ** 2, q)
            self.assertGreaterEqual(Fraction.from_float(r.hi) ** 2, q)


if __name__ == "__main__":
    unittest.main()
