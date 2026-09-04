from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_sea3_riccati_metric_p4 as mod  # noqa: E402


class Sea3MovingRiccatiP4Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_p4_uses_same_moving_metric_and_not_800_endpoint_scan(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(
            self.d["canonical_P4_architecture"],
            "NONLINEAR_WORD_IN_MOVING_SHIPPING_RICCATI_METRIC",
        )
        self.assertFalse(self.d["old_800_endpoint_signed_Joseph_scan_consumed"])
        self.assertFalse(self.d["old_terminal_source_phase_metric_attachment_consumed"])
        self.assertFalse(self.d["old_group_isotropic_P3_P4_metric_assumed"])
        self.assertGreaterEqual(self.d["outer_angle_rad"], 0.80)

    def test_p4_consumes_real_p3_and_only_reports_live_blockers(self):
        self.assertIsInstance(self.d["P3_CANONICAL_PASS_consumed"], bool)
        self.assertGreater(self.d["P3_H_delta_consumed"], 0.0)
        self.assertGreater(self.d["P3_A_delta_consumed"], 0.0)
        self.assertFalse(self.d["P4_CANONICAL_PASS"])
        self.assertFalse(self.d["P5_MAY_START"])
        self.assertTrue(self.d["exact_vector_accelerometer_congruence_rebind_pending"])
        self.assertTrue(self.d["P4_CANONICAL_FAIL_REASONS"])

        reasons = " ".join(self.d["P4_CANONICAL_FAIL_REASONS"]).lower()
        if self.d["P3_CANONICAL_PASS_consumed"]:
            self.assertNotIn("p3 h/a quantitative margin", reasons)
        else:
            self.assertIn("p3 h/a quantitative margin", reasons)


if __name__ == "__main__":
    unittest.main()
