import unittest
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_startup_direction_sampling as GAP


class DirectionSamplingGapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = GAP.build()

    def test_continuous_physical_witness_and_sampled_separation(self):
        result = self.result
        self.assertEqual(GAP.validate(result), [])
        exact = {k: F(v) for k, v in result['exact'].items()}
        self.assertLess(exact['continuous_mean_norm_upper'], F(1, 10))
        self.assertLess(exact['continuous_primitive_norm_upper_s'], F(3, 2))
        self.assertGreater(exact['sampled_mean_norm2_lower'], F(1, 100))
        # The separation is well above the outward quadrature uncertainty.
        self.assertGreater(exact['sampled_mean_norm2_excess_lower'], F(3, 100000))
        self.assertLess(exact['exact_real_guard_detector_abs_upper'], F(3, 100))

    def test_certified_pi_interval_and_periodic_ancestry(self):
        a, b = GAP.pi_interval()
        self.assertGreater(a, F('3.141592653589793238462643383279'))
        self.assertLess(b, F('3.141592653589793238462643383280'))
        self.assertEqual(GAP.DT*GAP.COUNT, GAP.T)
        self.assertEqual(GAP.RESIDUAL_CYCLES_PER_PERIOD, GAP.COUNT)

    def test_capture_cannot_be_promoted_from_the_obstruction(self):
        for key in ('shipping_startup_nonreachability_proved',
                    'universal_startup_capture_closed',
                    'storage_search_allowed'):
            modified = dict(self.result, **{key: True})
            self.assertTrue(GAP.validate(modified))


if __name__ == '__main__':
    unittest.main()
