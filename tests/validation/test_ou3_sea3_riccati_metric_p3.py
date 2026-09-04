from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_sea3_riccati_metric_p3 as mod  # noqa: E402


class Sea3MovingRiccatiP3Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_architecture_is_moving_shipping_covariance(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(
            self.d["canonical_P3_architecture"],
            "MOVING_SHIPPING_RICCATI_COVARIANCE_METRIC",
        )
        self.assertTrue(self.d["SEA3_dynamic_source_consumed"])
        self.assertFalse(self.d["old_P2_800_state_graph_consumed"])
        self.assertFalse(self.d["old_P2_V1_history_frontier_consumed"])
        self.assertFalse(self.d["old_terminal_source_phase_metric_attachment_consumed"])

    def test_open_numeric_obligation_is_fail_closed(self):
        self.assertTrue(self.d["P3_FOUNDATION_PASS"])
        self.assertFalse(self.d["P3_CANONICAL_PASS"])
        self.assertFalse(self.d["P4_MAY_CONSUME_P3"])
        self.assertEqual(self.d["useful_gate"], 1e-18)
        self.assertTrue(self.d["P3_CANONICAL_FAIL_REASONS"])

    def test_exact_metric_change_needs_no_artificial_jump_penalty(self):
        self.assertFalse(self.d["metric_derivative_or_jump_penalty_required"])
        self.assertTrue(self.d["metric_change_handled_by_exact_Riccati_recursion"])


if __name__ == "__main__":
    unittest.main()
