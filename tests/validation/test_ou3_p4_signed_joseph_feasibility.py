from pathlib import Path
import math
import sys
import unittest

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import ou3_p4_signed_joseph_feasibility as J


class Ou3P4SignedJosephFeasibilityTests(unittest.TestCase):
    def test_accelerometer_innovation_bound_retains_psd_cross_covariance(self):
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
        # H: (2*sqrt(4)+sqrt(9))^2 = 49.  A additionally pays sqrt(16),
        # giving 121.  Outward rounding may only increase these ceilings.
        self.assertGreaterEqual(h["accelerometer_state_innovation_covariance_upper"], 49.0)
        self.assertGreaterEqual(a["accelerometer_state_innovation_covariance_upper"], 121.0)
        self.assertGreaterEqual(h["magnetometer_state_innovation_covariance_upper"], 36.0)
        self.assertGreaterEqual(a["magnetometer_state_innovation_covariance_upper"], 36.0)
        self.assertGreater(a["accelerometer_S_lambda_max_upper"], h["accelerometer_S_lambda_max_upper"])

    def test_measurement_variance_bounds_are_outward(self):
        vector = {
            "configured_measurement_bounds": {
                "acc_measurement_std_mps2": 2.0,
                "mag_measurement_std_uT": 3.0,
            }
        }
        r = J._measurement_variances(vector)
        self.assertLessEqual(r["acc_lower"], 4.0)
        self.assertGreaterEqual(r["acc_upper"], 4.0)
        self.assertLessEqual(r["mag_lower"], 9.0)
        self.assertGreaterEqual(r["mag_upper"], 9.0)
        self.assertTrue(all(math.isfinite(x) and x > 0.0 for x in r.values()))

    def test_signed_fraction_criterion_is_dimensionless(self):
        # At q=0.8 the exact pure-rotation eta/y squared ratio is q^2/4=0.16.
        # A local Joseph scalarization survives iff R/S exceeds that value.
        q = 0.8
        eta = J.up(J.mul_up(q, q) / 4.0)
        good = J.down(0.25 - eta)
        bad = J.down(0.10 - eta)
        self.assertGreater(good, 0.0)
        self.assertLess(bad, 0.0)

    def test_a_mode_requires_accel_bias_upper(self):
        env = {
            "translation_covariance_upper_groups": {"a_w": 1.0},
            "A_bias_covariance_upper": {
                "theta_covariance_upper": 1.0,
                "accel_bias_covariance_upper": None,
            },
        }
        with self.assertRaises((TypeError, ValueError, RuntimeError)):
            J._innovation_bounds(env, "A", 1.0, 1.0, {"acc_upper": 1.0, "mag_upper": 1.0})


if __name__ == "__main__":
    unittest.main()
