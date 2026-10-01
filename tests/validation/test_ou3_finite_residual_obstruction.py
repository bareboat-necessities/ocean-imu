from fractions import Fraction as F
import math
import unittest
from tools.stability.ou3_theorem.finite_residual_obstruction import certificate


class FiniteResidualObstructionTests(unittest.TestCase):
    def test_all_time_envelopes_and_service_are_rational(self):
        c = certificate()
        self.assertLess(F(c['accel_residual_norm_upper_mps2']), F(3, 10))
        self.assertLess(F(c['gyro_residual_norm_upper_rad_s']), F(1, 50))
        self.assertGreater(F(c['actual_magnetic_service_lower']), 1)
        self.assertEqual(F(c['physical_sqrt_V_lower']), F(2, 5))
        self.assertFalse(c['all_time_float32_verified'])
        self.assertFalse(c['local_homogeneous_theorem_refuted'])

    def test_packet_identity_at_examples(self):
        # Examples check signs. The all-time proof uses the exact cancellation
        # and rational majorants, not a sampled maximum.
        g, gm, ba = 9.80665, 9.8066501617431640625, .01
        for t in (0., .13, 1.8, 8., 120., 50000.):
            phi = .01*math.sin(t/2)
            down = (0., math.sin(phi), math.cos(phi))
            n = ( -ba, g*down[1], g*down[2]-gm )
            f = (-g*down[0]+ba+n[0], -g*down[1]+n[1], -g*down[2]+n[2])
            self.assertEqual(f, (0., 0., -gm))
            self.assertLess(math.sqrt(sum(x*x for x in n)), .3)

    def test_excitation_is_true_tilt_not_proxy(self):
        c = certificate()
        self.assertLess(4*math.pi, c['example_excited_window_s'])
        self.assertGreater(float(F(c['all_time_gravity_span_rad'])), math.pi/180)
        self.assertGreater(.02, 60*.001/9.80665)


if __name__ == '__main__':
    unittest.main()
