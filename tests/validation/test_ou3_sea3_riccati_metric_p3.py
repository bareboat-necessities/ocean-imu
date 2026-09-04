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
        self.assertTrue(self.d["quantitative_Riccati_tube_consumed"])
        self.assertTrue(self.d["current_source_interval_cover_only"])
        self.assertFalse(self.d["source_history_graph_consumed"])
        self.assertFalse(self.d["predecessor_path_enumeration_consumed"])
        self.assertFalse(self.d["old_P2_800_state_graph_consumed"])
        self.assertFalse(self.d["old_P2_V1_history_frontier_consumed"])
        self.assertFalse(self.d["old_terminal_source_phase_metric_attachment_consumed"])
        self.assertEqual(self.d["cell_cover"]["history_depth"], 0)

    def test_canonical_verdict_is_exactly_quantitative_tube_verdict(self):
        self.assertTrue(self.d["P3_FOUNDATION_PASS"])
        self.assertEqual(self.d["useful_gate"], 1e-18)

        expected = True
        for mode in ("H", "A"):
            row = self.d["modes"][mode]
            mode_expected = row["relative_Riccati_injection_margin_lower"] >= self.d["useful_gate"]
            expected = expected and mode_expected
            self.assertTrue(row["riccati_covariance_upper_bound_closed"])
            self.assertTrue(row["word_injection_comparison_closed"])
            self.assertEqual(row["pass"], mode_expected)
            self.assertGreater(row["Pbar_lambda_max_trace_upper"], 0.0)

        self.assertEqual(self.d["quantitative_Riccati_tube_pass"], expected)
        self.assertEqual(self.d["P3_CANONICAL_PASS"], expected)
        self.assertEqual(self.d["P4_MAY_CONSUME_P3"], expected)
        if expected:
            self.assertEqual(self.d["P3_CANONICAL_FAIL_REASONS"], [])
        else:
            self.assertTrue(self.d["P3_CANONICAL_FAIL_REASONS"])

    def test_exact_metric_change_needs_no_artificial_jump_penalty(self):
        self.assertFalse(self.d["metric_derivative_or_jump_penalty_required"])
        self.assertTrue(self.d["metric_change_handled_by_exact_Riccati_recursion"])


if __name__ == "__main__":
    unittest.main()
