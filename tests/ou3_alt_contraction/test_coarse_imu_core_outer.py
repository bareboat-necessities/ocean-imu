import unittest
from tools.stability.ou3_alt_contraction import coarse_imu_core_outer as C


class CoarseImuCoreOuterTest(unittest.TestCase):
    def test_mag_is_not_hidden_in_finite_sample_family(self):
        d=C.build();self.assertEqual(C.validate(d),[])
        print('ALT_COARSE_IMU_CORE',C.summary(d))
        self.assertTrue(d['async_magnetometer_excluded_from_finite_IMU_core_family'])
        self.assertFalse(d['async_magnetometer_count_upper_bound_assumed'])
        self.assertTrue(d['async_magnetometer_star_theorem_consumed'])
        self.assertFalse(d['async_magnetometer_nonexpansive_same_M_closed'])
        self.assertFalse(d['complete_async_sample_family_closed'])
        self.assertFalse(d['ALT_LIVE_PASS'])


if __name__=='__main__':unittest.main()
