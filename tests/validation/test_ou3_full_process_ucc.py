import importlib.util
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location(
    "ou3_full_process_ucc", ROOT / "tools" / "ou3_full_process_ucc.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class FullProcessUccTests(unittest.TestCase):
    def test_complete_H_A_process_covariance_is_strict(self):
        d = mod.build()
        self.assertEqual(mod.validate(d), [])
        self.assertTrue(d["pass"])
        self.assertTrue(d["full_process_ucc_pass"])
        for mode, dim in (("H", 18), ("A", 21)):
            self.assertEqual(d["modes"][mode]["dimension"], dim)
            self.assertTrue(d["modes"][mode]["pass"])
            self.assertGreater(d["modes"][mode]["prediction_Q_lambda_min_lower"], 0.0)

    def test_attitude_bias_lower_bound_is_rate_independent_and_positive(self):
        d = mod.build()
        a = d["attitude_gyro_bias"]
        self.assertGreater(a["q_gyro_lower"], 0.0)
        self.assertGreater(a["q_gyro_bias_lower"], 0.0)
        self.assertGreater(a["Q_attitude_gyro_bias_lambda_min_lower"], 0.0)
        self.assertLess(a["cross_norm_upper"], a["gyro_bias_diagonal_lower"])

    def test_active_accelerometer_bias_uses_validated_expm1(self):
        d = mod.build()
        b = d["active_accelerometer_bias"]
        self.assertGreaterEqual(b["two_h_over_tau_interval"][0], 0.0)
        self.assertLess(b["expm1_minus_interval"][1], 0.0)
        self.assertGreater(b["qd_scale_interval_s"][0], 0.0)
        self.assertGreater(b["Q_accel_bias_lambda_min_lower"], 0.0)

    def test_translation_is_the_limiting_but_strict_source_block(self):
        d = mod.build()
        qt = d["translation"]["Q_translation_lambda_min_lower"]
        qa = d["attitude_gyro_bias"]["Q_attitude_gyro_bias_lambda_min_lower"]
        qb = d["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"]
        self.assertGreater(qt, 0.0)
        self.assertLess(qt, qa)
        self.assertLess(qt, qb)
        self.assertLessEqual(d["modes"]["H"]["prediction_Q_lambda_min_lower"], qt)
        self.assertLessEqual(d["modes"]["A"]["prediction_Q_lambda_min_lower"], qt)

    def test_no_replay_inputs(self):
        text = (ROOT / "tools" / "ou3_full_process_ucc.py").read_text()
        for forbidden in ("ou3_exact_replay", "path_metrics", "neighborhood_radius_search"):
            self.assertNotIn(forbidden, text)


if __name__ == "__main__":
    unittest.main()
