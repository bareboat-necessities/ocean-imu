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

    def test_p4_remains_blocked_until_real_p3_and_nonlinear_margin_close(self):
        self.assertFalse(self.d["P3_CANONICAL_PASS_consumed"])
        self.assertFalse(self.d["P4_CANONICAL_PASS"])
        self.assertFalse(self.d["P5_MAY_START"])
        self.assertTrue(self.d["P4_CANONICAL_FAIL_REASONS"])


if __name__ == "__main__":
    unittest.main()
