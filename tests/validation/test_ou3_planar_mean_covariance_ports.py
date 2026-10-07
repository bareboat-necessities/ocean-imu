from fractions import Fraction as F
import unittest
from tools.stability.ou3_theorem.planar_mean_covariance_ports import certificate


class MeanCovariancePortTests(unittest.TestCase):
    def test_missing_cell_bounds_fail_closed(self):
        c = certificate()
        self.assertIsNone(c['acc_H_attitude_coefficient_squared'])
        self.assertIsNone(c['innovation_inverse_bound'])
        self.assertFalse(c['cell_bounds_supplied'])
        self.assertFalse(c['physical_acceleration_substituted_for_nominal_aw'])

    def test_conditional_rational_cell_bounds(self):
        c = certificate(g=F(10), nominal_aw_bound=F(2),
                        reference_norm_bound=F(76), measurement_noise_floor=F(1, 25))
        self.assertEqual(c['acc_H_attitude_coefficient_squared'], '145')
        self.assertEqual(c['innovation_inverse_bound'], '25')
        self.assertFalse(c['cell_bounds_certified'])
        self.assertFalse(c['joint_cell_forward_invariant'])
        with self.assertRaises(ValueError):
            certificate(measurement_noise_floor=F(0))


if __name__ == '__main__':
    unittest.main()
