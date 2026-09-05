from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_sea3_riccati_metric_p3 as mod  # noqa: E402


class Sea3LiftedFiniteWordP3Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_architecture_is_lifted_finite_word(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(
            self.d["canonical_P3_architecture"],
            "LIFTED_FINITE_WORD_SELECTED_PROCESS_MODES",
        )
        self.assertTrue(self.d["P3_ARCHITECTURE_READY"])
        self.assertTrue(self.d["samplewise_nonexpansion_closed"])
        self.assertEqual(self.d["strictness_location"], "RECURRENT_FINITE_WORD_ONLY")
        self.assertFalse(self.d["source_history_graph_consumed"])
        self.assertFalse(self.d["predecessor_path_enumeration_consumed"])
        self.assertFalse(self.d["old_P2_800_state_graph_consumed"])
        self.assertFalse(self.d["old_P2_V1_history_frontier_consumed"])

    def test_dead_end_numerical_routes_cannot_reenter(self):
        self.assertFalse(self.d["one_sample_strict_Riccati_margin_consumed"])
        self.assertFalse(self.d["commit_aligned_source_word_consumed"])
        self.assertFalse(self.d["per_sample_SPD_lower_required"])
        self.assertFalse(self.d["determinant_trace_scalarization_used"])
        self.assertFalse(self.d["scalar_information_beta_used"])
        self.assertEqual(self.d["useful_gate"], 1e-18)

    def test_gate_fails_closed_until_full_lifted_matrices_close(self):
        self.assertFalse(self.d["P3_QUANTITATIVE_WORD_MATRIX_CLOSED"])
        self.assertFalse(self.d["P3_CANONICAL_PASS"])
        self.assertFalse(self.d["P4_MAY_CONSUME_P3"])
        self.assertTrue(self.d["P3_CANONICAL_FAIL_REASONS"])
        for mode in ("H", "A"):
            row = self.d["modes"][mode]
            self.assertEqual(row["selected_process_mode_dimension"], row["dimension"])
            self.assertFalse(row["lifted_endpoint_map_B_closed"])
            self.assertFalse(row["lifted_measurement_information_J_upper_closed"])
            self.assertFalse(row["pass"])
            self.assertEqual(row["relative_Riccati_injection_margin_lower"], 0.0)


if __name__ == "__main__":
    unittest.main()
