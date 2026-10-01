import math
import unittest

from tools.stability.ou3_theorem.stationary_covariance import quiet_covariance_bounds


class OuterEntryObstructionTests(unittest.TestCase):
    def test_constant_rest_attitude_ba_gauge_excludes_point_inner_ball(self):
        g = 9.80665
        alpha = 1.0e-3
        b_a = 2.0 * g * math.sin(alpha / 2.0)
        self.assertLess(b_a, 0.22516660498395405)

        # The source-audited quiet comparison has P_ba,ba <= (1/40)^2 I.
        self.assertTrue(quiet_covariance_bounds()["full_21_upper_diagonal"])
        sqrt_v_ba_lower = b_a / (1.0 / 40.0)
        self.assertGreater(sqrt_v_ba_lower, 0.392)
        self.assertGreater(sqrt_v_ba_lower, 0.15)

    def test_constant_rotation_about_field_is_packet_indistinguishable(self):
        # Choose B along x and rotate about x. Q B=B; compensating BA makes
        # -g Q e_z+b_a=-g e_z exactly.
        g = 9.80665
        a = 1.0e-3
        c, s = math.cos(a), math.sin(a)
        qez = (0.0, -s, c)
        ez = (0.0, 0.0, 1.0)
        ba = tuple(g * (qez[i] - ez[i]) for i in range(3))
        f = tuple(-g * qez[i] + ba[i] for i in range(3))
        self.assertAlmostEqual(f[0], 0.0, places=14)
        self.assertAlmostEqual(f[1], 0.0, places=14)
        self.assertAlmostEqual(f[2], -g, places=14)


if __name__ == "__main__":
    unittest.main()
