import unittest
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_theorem.planar_frontend_domain import certificate
from tools.stability.ou3_theorem.mahony_raw_normalization import certificate as norm_certificate


class FrontendDomainTests(unittest.TestCase):
    def test_rational_source_charges_close_subordinate_domain(self):
        c = certificate()
        self.assertGreater(F(c['strict_radial_margin']), 0)
        self.assertLess(F(c['source_derived_one_step_charges']['pitch_charge']), F(46, 10**9))
        self.assertLess(F(c['source_derived_one_step_charges']['scaled_integral_charge']), F(2, 10**9))
        self.assertFalse(c['inherited_nominal_mekf_domain_verified'])
        self.assertFalse(c['target_toolchain_libm_qualification_verified'])
        self.assertFalse(c['complete_word_storage_contraction_verified'])
        self.assertFalse(c['all_time_magnetic_service_verified'])
        ref = c['reference_geometry_real']
        self.assertTrue(ref['source_magnetic_rounding_charge_retained'])
        self.assertEqual(F(ref['norm_upper']), 75*(1+F(1, 2**24)))

    def test_accelerometer_normalization_is_in_proved_exponent_domain(self):
        c = norm_certificate()
        self.assertEqual(c['extended_input_squared_norm_domain'], ['1/4', '128'])
        lo, hi = map(F, c['pre_component_scalar_squared_norm_bounds'])
        self.assertGreater(lo, F(995, 1000))
        self.assertLess(hi, F(1001, 1000))

    def test_source_operation_order_and_guard_initialization(self):
        root = Path(__file__).resolve().parents[2]
        text = (root/'src/ahrs/Mahony_AHRS.h').read_text()
        imu = text[text.index('void update('):text.index('void updateMag(')]
        self.assertLess(imu.index('integralFBy +='), imu.index('gy += integralFBy'))
        self.assertLess(imu.index('gy += integralFBy'), imu.index('gy += twoKp * halfey'))
        self.assertIn('q0 += (-qb * gx - qc * gy - q3 * gz)', imu)
        self.assertIn('q2 += ( qa * gy - qb * gz + q3 * gx)', imu)
        guard = (root/'src/tuner/AccelVibrationGuard.h').read_text()
        self.assertIn('detect_stages_[0] = acc;', guard)
        self.assertIn('high_passed -= detect_stages_[i];', guard)
        self.assertIn('if (weight_ <= 0.0f) return acc;', guard)


if __name__ == '__main__':
    unittest.main()
