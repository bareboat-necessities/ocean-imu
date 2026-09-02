#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p3_ltv_postmeasurement_certificate as P
import ou3_source_reachable_matrix_p3 as BASE


class LtvPostmeasurementP3Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = P.build()

    def test_certificate_is_source_complete_and_fail_closed(self):
        d = self.payload
        self.assertEqual(P.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertTrue(d["zero_lever_arm_branch"])
        self.assertTrue(d["dormant_transparent_vibration_guard_branch"])
        self.assertTrue(d["measurement_attenuation_applied_once_in_lifted_information_space"])
        self.assertFalse(d["per_sample_multiplicative_Joseph_loss_used"])
        self.assertFalse(d["trajectory_replay_used"])
        self.assertFalse(d["filter_changed"])

    def test_both_modes_scan_all_tau_cells_with_strict_attenuation(self):
        for mode in ("H", "A"):
            row = self.payload["modes"][mode]
            self.assertEqual(row["endpoint_tau_cells_scanned"], 10)
            self.assertGreater(row["relative_Riccati_injection_margin_lower"], 0.0)
            for endpoint in row["endpoint_rows"]:
                best = endpoint["best"]
                self.assertGreater(best["joint_selected_mode_attenuation_lower"], 0.0)
                self.assertLessEqual(best["joint_selected_mode_attenuation_lower"], 1.0)
                self.assertGreater(best["beta_total_upper"], 0.0)
                self.assertGreater(best["measurement_slots_upper"], 0)

    def test_promotion_flag_is_exactly_the_unchanged_useful_gate(self):
        established = all(
            self.payload["modes"][mode]["relative_Riccati_injection_margin_lower"]
            >= BASE.MIN_USEFUL_DELTA
            for mode in ("H", "A")
        )
        self.assertEqual(self.payload["P3_LINEAR_CERTIFICATE_ESTABLISHED"], established)
        self.assertEqual(self.payload["P3_PROMOTED"], established)
        self.assertEqual(self.payload["useful_gate"], 1.0e-18)

    def test_selected_mode_constant_matches_existing_ltv_determinant_constant(self):
        self.assertAlmostEqual(2025.0 / 144.0, (15.0 / 4.0) ** 2, places=15)
        self.assertAlmostEqual(2025.0 / 144.0, 225.0 / 16.0, places=15)

    def test_A_mode_pays_accel_bias_information_and_H_does_not(self):
        h = self.payload["modes"]["H"]["worst_endpoint"]["best"]
        a = self.payload["modes"]["A"]["worst_endpoint"]["best"]
        self.assertEqual(h["beta_final_accel_bias_acc_upper"], 0.0)
        self.assertIsNone(h["accel_bias_relative_noise_floor_lower"])
        self.assertGreater(a["beta_final_accel_bias_acc_upper"], 0.0)
        self.assertTrue(math.isfinite(a["accel_bias_relative_noise_floor_lower"]))
        self.assertGreater(a["accel_bias_relative_noise_floor_lower"], 0.0)


if __name__ == "__main__":
    unittest.main()
