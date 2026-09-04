from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_sea3_riccati_tube as mod  # noqa: E402


class Sea3RiccatiTubeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_no_history_graph_is_consumed(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertTrue(self.d["SEA3_dynamic_source_consumed"])
        self.assertTrue(self.d["current_source_interval_cover_only"])
        self.assertFalse(self.d["source_history_graph_consumed"])
        self.assertFalse(self.d["predecessor_path_enumeration_consumed"])
        self.assertFalse(self.d["P2_800_state_partition_consumed"])
        self.assertEqual(self.d["cell_cover"]["history_depth"], 0)

    def test_margins_are_real_positive_comparisons(self):
        self.assertEqual(self.d["useful_gate"], 1e-18)
        for mode in ("H", "A"):
            row = self.d["modes"][mode]
            self.assertGreater(row["relative_Riccati_injection_margin_lower"], 0.0)
            self.assertGreater(
                row["worst_current_source_cell"]["post_measurement_scaled_Omega_lambda_min_lower"],
                0.0,
            )
            self.assertGreater(row["Pbar_lambda_max_trace_upper"], 0.0)

    def test_cross_covariances_are_paid_by_trace_bound(self):
        text = self.d["PSD_cross_covariance_handling"].lower()
        self.assertIn("trace", text)
        self.assertIn("no max-diagonal shortcut", text)


if __name__ == "__main__":
    unittest.main()
