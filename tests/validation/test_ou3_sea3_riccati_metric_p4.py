from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p4_moving_metric_rebind as rebind  # noqa: E402
import ou3_sea3_riccati_metric_p4 as mod  # noqa: E402


class Sea3MovingRiccatiP4Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()
        cls.r = rebind.build()

    def test_p4_uses_same_moving_metric_and_not_800_endpoint_scan(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(rebind.validate(self.r), [])
        self.assertEqual(
            self.d["canonical_P4_architecture"],
            "NONLINEAR_WORD_IN_MOVING_SHIPPING_RICCATI_METRIC",
        )
        self.assertFalse(self.d["old_800_endpoint_signed_Joseph_scan_consumed"])
        self.assertFalse(self.d["old_terminal_source_phase_metric_attachment_consumed"])
        self.assertFalse(self.d["old_group_isotropic_P3_P4_metric_assumed"])
        self.assertGreaterEqual(self.d["outer_angle_rad"], 0.80)

    def test_exact_moving_metric_rebind_is_closed(self):
        self.assertTrue(self.d["P3_CANONICAL_PASS_consumed"])
        self.assertTrue(self.d["full_nonlinear_measurement_metric_rebind_closed"])
        self.assertFalse(self.d["exact_vector_accelerometer_congruence_rebind_pending"])
        self.assertTrue(self.d["moving_metric_coordinate_congruence_exact"])
        self.assertTrue(self.d["Joseph_nonlinear_injection_metric_closed"])
        self.assertTrue(self.r["prediction_linear_map_nonexpansive"])
        self.assertTrue(self.r["Joseph_linear_map_nonexpansive"])
        self.assertTrue(self.r["left_error_reset_exact_metric_isometry"])
        self.assertFalse(self.r["group_isotropic_metric_attachment_used"])
        self.assertFalse(self.r["endpoint_source_word_scan_used"])

    def test_only_live_blocker_is_finite_nonlinear_remainder(self):
        self.assertGreaterEqual(self.d["P3_H_delta_consumed"], 1.0e-18)
        self.assertGreaterEqual(self.d["P3_A_delta_consumed"], 1.0e-18)
        self.assertFalse(self.d["P4_CANONICAL_PASS"])
        self.assertFalse(self.d["P5_MAY_START"])
        self.assertFalse(self.d["nonlinear_remainder_dominated_on_full_sector"])
        self.assertEqual(1, len(self.d["P4_CANONICAL_FAIL_REASONS"]))
        self.assertIn(
            "nonlinear remainder",
            self.d["P4_CANONICAL_FAIL_REASONS"][0].lower(),
        )


if __name__ == "__main__":
    unittest.main()
