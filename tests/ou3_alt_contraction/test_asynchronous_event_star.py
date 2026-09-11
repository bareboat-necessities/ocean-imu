import unittest
from tools.stability.ou3_alt_contraction import asynchronous_event_star as A


class AsynchronousEventStarTest(unittest.TestCase):
    def test_finite_nonexpansive_star_needs_no_count_cap(self):
        self.assertEqual(A.finite_star_rho(1.0,0),1.0)
        self.assertEqual(A.finite_star_rho(1.0,1000000),1.0)
        self.assertLess(A.finite_star_rho(0.9,10),1.0)

    def test_shipping_kernel_requires_star_not_zero_one_shortcut(self):
        d=A.build();self.assertEqual(A.validate(d),[])
        self.assertTrue(d['typed_kernel_executes_every_async_event'])
        self.assertFalse(d['magnetometer_event_count_upper_bound_assumed'])
        self.assertFalse(d['zero_or_one_mag_per_IMU_is_universal_cover'])
        self.assertTrue(d['finite_star_composition_theorem_closed'])
        self.assertFalse(d['single_source_uniform_magnetometer_nonexpansive_certificate_closed'])
        self.assertFalse(d['arbitrary_finite_async_magnetometer_sequence_closed'])
        self.assertFalse(d['ALT_LIVE_PASS'])


if __name__=='__main__':unittest.main()
