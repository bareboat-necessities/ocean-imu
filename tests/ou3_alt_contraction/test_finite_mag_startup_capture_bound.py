"""Deterministic startup magnetic capture-bound regressions."""
from fractions import Fraction as F
from pathlib import Path
import sys, unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))

from tools.stability.ou3_alt_contraction import finite_mag_startup_capture_bound as X


class Tests(unittest.TestCase):
    def test_current_envelope_exact_tilt_margin_is_one_over_75_half_angle_sine(self):
        self.assertEqual(X.default_tilt_margin_half_sin(),F(1,75))
        at=X.derive(F(1,75))
        self.assertEqual(at.earth_rotation_error_max,2)
        self.assertEqual(at.mean_perturbation_max,15)
        self.assertFalse(at.capture_nonzero)
        with self.assertRaisesRegex(ValueError,'erase horizontal north'):
            X.require_capture(at)

    def test_strictly_inside_margin_gives_source_uniform_yaw_sine_bound(self):
        b=X.derive(F(1,100))
        self.assertEqual(b.earth_rotation_error_max,F(3,2))
        self.assertEqual(b.mean_perturbation_max,F(29,2))
        self.assertEqual(b.sin_yaw_error_max,F(29,30))
        self.assertTrue(b.capture_nonzero); self.assertTrue(X.require_capture(b))

    def test_perfect_tilt_still_retains_full_deterministic_hardiron_and_noise(self):
        b=X.derive(0)
        self.assertEqual(b.mean_perturbation_max,13)
        self.assertEqual(b.sin_yaw_error_max,F(13,15))
        self.assertTrue(b.capture_nonzero)
        # There is intentionally no 1/sqrt(250) factor from tuner sample count.
        self.assertNotEqual(b.mean_perturbation_max,F(13,250))

    def test_no_fake_capture_when_tilt_bound_is_too_large(self):
        b=X.derive(F(1,50))
        self.assertEqual(b.earth_rotation_error_max,3)
        self.assertEqual(b.mean_perturbation_max,16)
        self.assertFalse(b.capture_nonzero)
        with self.assertRaisesRegex(ValueError,'erase horizontal north'):
            X.require_capture(b)

    def test_readiness_leaves_mahony_and_binary32_open(self):
        r=X.readiness()
        self.assertTrue(r['deterministic_average_does_not_claim_sqrtN_improvement'])
        self.assertTrue(r['current_envelope_requires_half_tilt_sin_lt_1_over_75'])
        self.assertFalse(r['Mahony_startup_tilt_bound_attached'])
        self.assertFalse(r['atan2_binary32_yaw_bound_attached'])
        self.assertFalse(r['startup_capture_closed']); self.assertFalse(r['ALT_STARTUP_PASS'])

if __name__=='__main__': unittest.main()
