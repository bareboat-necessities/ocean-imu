from pathlib import Path
import math
import sys
import unittest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import ou3_p4_p3_metric_attachment as M
import ou3_p4_signed_joseph_feasibility as J


class Ou3P4P3MetricAttachmentTests(unittest.TestCase):
    def test_shipping_H_A_state_order_is_exact(self):
        self.assertEqual(len(M.STATE_ORDER["H"]), 18)
        self.assertEqual(len(M.STATE_ORDER["A"]), 21)
        self.assertEqual(M.STATE_GROUPS["H"], ["theta", "b_g", "v", "p", "S", "a_w"])
        self.assertEqual(M.STATE_GROUPS["A"], ["theta", "b_g", "v", "p", "S", "a_w", "b_a"])
        self.assertEqual(M.STATE_ORDER["H"][6:18], [
            "v_x", "v_y", "v_z", "p_x", "p_y", "p_z",
            "S_x", "S_y", "S_z", "a_wx", "a_wy", "a_wz",
        ])
        self.assertEqual(M.STATE_ORDER["A"][18:21], ["b_ax", "b_ay", "b_az"])

    def test_scaled_translation_floor_maps_back_to_physical_units(self):
        rho = 2.0
        h = 0.5
        g = M._translation_physical_lower_groups(rho, h)
        self.assertGreater(g["v"], 0.0)
        self.assertGreater(g["p"], 0.0)
        self.assertGreater(g["S"], 0.0)
        self.assertGreater(g["a_w"], 0.0)
        self.assertLessEqual(g["v"], rho * h**2)
        self.assertLessEqual(g["p"], rho * h**4)
        self.assertLessEqual(g["S"], rho * h**6)
        self.assertLessEqual(g["a_w"], rho)

    def test_finite_H_metric_pays_exact_precision_join_factor_once(self):
        cond = {
            "attitude_conditional_posterior_lower": 4.0,
            "gyro_bias_conditional_posterior_lower": 9.0,
        }
        d = M._finite_mode_metric("H", rho_z=2.0, h=0.5, cond=cond)
        self.assertEqual(d["dimension"], 18)
        self.assertEqual(d["each_group_repeated_coordinates"], 3)
        lower = d["covariance_lower_group_diagonal"]
        info = d["information_metric_upper_group_diagonal"]
        self.assertLessEqual(lower["theta"], 2.0)
        self.assertLessEqual(lower["b_g"], 4.5)
        self.assertGreater(lower["theta"], 0.0)
        self.assertGreater(lower["b_g"], 0.0)
        for group in M.STATE_GROUPS["H"]:
            self.assertTrue(math.isfinite(lower[group]) and lower[group] > 0.0)
            self.assertTrue(math.isfinite(info[group]) and info[group] > 0.0)
            self.assertGreaterEqual(info[group], 1.0 / lower[group])

    def test_active_accel_bias_mode_remains_21_dimensional(self):
        cond = {
            "attitude_conditional_posterior_lower": 4.0,
            "gyro_bias_conditional_posterior_lower": 9.0,
            "accel_bias_conditional_posterior_lower": 16.0,
        }
        d = M._finite_mode_metric("A", rho_z=2.0, h=0.5, cond=cond)
        self.assertEqual(d["dimension"], 21)
        self.assertIn("b_a", d["covariance_lower_group_diagonal"])
        self.assertIn("b_a", d["information_metric_upper_group_diagonal"])
        self.assertLessEqual(d["covariance_lower_group_diagonal"]["b_a"], 8.0)

    def test_frozen_metric_uses_certified_relative_floor_not_covariance_inverse(self):
        cond = {
            "attitude_conditional_posterior_lower": 4.0,
            "gyro_bias_conditional_posterior_lower": 9.0,
        }
        delta = 0.25
        upper = [8.0, 12.0, 20.0, 5.0]
        d = M._frozen_mode_metric("H", delta, upper, cond)
        lower = d["covariance_lower_group_diagonal"]
        # P_T >= delta*diag(U_T), then the same P3 precision join pays 1/2.
        self.assertLessEqual(lower["v"], 0.5 * delta * upper[0])
        self.assertLessEqual(lower["p"], 0.5 * delta * upper[1])
        self.assertLessEqual(lower["S"], 0.5 * delta * upper[2])
        self.assertLessEqual(lower["a_w"], 0.5 * delta * upper[3])
        self.assertTrue(all(x > 0.0 for x in lower.values()))

    def test_signed_joseph_accelerometer_bound_retains_psd_cross_covariance(self):
        env = {
            "translation_covariance_upper_groups": {"a_w": 9.0},
            "H_bias_covariance_upper": {
                "theta_covariance_upper": 4.0,
                "accel_bias_covariance_upper": None,
            },
            "A_bias_covariance_upper": {
                "theta_covariance_upper": 4.0,
                "accel_bias_covariance_upper": 16.0,
            },
        }
        rv = {"acc_upper": 1.0, "mag_upper": 1.0}
        h = J._innovation_bounds(env, "H", fmax=2.0, mmax=3.0, rv=rv)
        a = J._innovation_bounds(env, "A", fmax=2.0, mmax=3.0, rv=rv)
        # Arbitrary PSD cross blocks are retained through
        # (|f|sqrt(Utheta)+sqrt(Uaw)[+sqrt(Uba)])^2.
        self.assertGreaterEqual(h["accelerometer_state_innovation_covariance_upper"], 49.0)
        self.assertGreaterEqual(a["accelerometer_state_innovation_covariance_upper"], 121.0)
        self.assertGreaterEqual(h["magnetometer_state_innovation_covariance_upper"], 36.0)
        self.assertGreater(a["accelerometer_S_lambda_max_upper"], h["accelerometer_S_lambda_max_upper"])

    def test_signed_joseph_local_scalar_criterion_has_correct_sign(self):
        q = 0.8
        eta_over_y2 = J.up(J.mul_up(q, q) / 4.0)
        self.assertGreater(J.down(0.25 - eta_over_y2), 0.0)
        self.assertLess(J.down(0.10 - eta_over_y2), 0.0)

    def test_metric_module_has_no_theorem_promotion_shortcut(self):
        self.assertEqual(M.JOIN_FACTOR, 0.5)
        self.assertEqual(M.PHASES, tuple(range(26)))
        self.assertEqual(M.BASE.MIN_USEFUL_DELTA, 1.0e-18)


if __name__ == "__main__":
    unittest.main()
