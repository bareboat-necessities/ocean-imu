"""Deterministic startup magnetic capture-bound regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_startup_capture_bound as X

class Tests(unittest.TestCase):
    def test_current_envelope_exact_open_tilt_margin_is_four_over_75(self):
        self.assertEqual(X.default_tilt_margin_half_sin(),F(4,75))
        at=X.derive(F(4,75))
        self.assertEqual(at.earth_rotation_error_max,8)
        self.assertEqual(at.mean_perturbation_max,15)
        self.assertFalse(at.capture_nonzero)
        with self.assertRaisesRegex(ValueError,'erase horizontal north'): X.require_capture(at)

    def test_declared_0p02_rad_startup_tilt_supplies_exact_rational_capture_bound(self):
        b=X.declared_startup_tilt_capture()
        self.assertEqual(X.DECLARED_STARTUP_TILT_RAD_MAX,F(1,50))
        self.assertEqual(b.half_tilt_sin_max,F(1,100))
        self.assertEqual(b.earth_rotation_error_max,F(3,2))
        self.assertEqual(b.mean_perturbation_max,F(17,2))
        self.assertEqual(b.sin_yaw_error_max,F(17,30))
        self.assertTrue(b.capture_nonzero); self.assertTrue(X.require_capture(b))

    def test_perfect_tilt_still_retains_full_deterministic_hardiron_and_noise(self):
        b=X.derive(0)
        self.assertEqual(b.mean_perturbation_max,7)
        self.assertEqual(b.sin_yaw_error_max,F(7,15))
        self.assertTrue(b.capture_nonzero)
        self.assertNotEqual(b.mean_perturbation_max,F(7,250))

    def test_no_fake_capture_at_large_tilt(self):
        b=X.derive(F(3,50))
        self.assertEqual(b.earth_rotation_error_max,9)
        self.assertEqual(b.mean_perturbation_max,16)
        self.assertFalse(b.capture_nonzero)
        with self.assertRaisesRegex(ValueError,'erase horizontal north'): X.require_capture(b)

    def test_readiness_attaches_mahony_tilt_but_leaves_binary32_and_SO3_comparison_open(self):
        r=X.readiness()
        self.assertTrue(r['deterministic_average_does_not_claim_sqrtN_improvement'])
        self.assertTrue(r['current_envelope_half_tilt_sin_margin_is_4_over_75'])
        self.assertTrue(r['declared_startup_tilt_0p02_rad_attached_via_sin_x_le_x'])
        self.assertTrue(r['declared_tilt_plus_magnetic_envelope_yields_E_8p5_uT'])
        self.assertTrue(r['declared_tilt_plus_magnetic_envelope_yields_sin_yaw_le_17_over_30'])
        self.assertFalse(r['atan2_binary32_yaw_bound_attached'])
        self.assertFalse(r['full_SO3_45deg_startup_entrance_closed'])
        self.assertFalse(r['startup_capture_closed']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
