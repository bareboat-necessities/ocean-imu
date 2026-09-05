from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_sea3_riccati_tube_factored as mod  # noqa: E402


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

    def test_endpoint_referenced_translation_ceiling_is_canonical(self):
        timing = self.d["covariance_memory"]
        self.assertTrue(self.d["endpoint_referenced_translation_covariance"])
        self.assertFalse(self.d["post_reconstruction_forward_propagation_used"])
        self.assertEqual(timing["translation_reference"], "word_endpoint")
        self.assertTrue(timing["endpoint_referenced_observability"])
        self.assertTrue(timing["endpoint_p_sign_similarity_applied"])
        self.assertFalse(timing["forward_propagation_after_endpoint_reconstruction"])
        self.assertTrue(timing["full_word_process_noise_dominator_retained"])

    def test_endpoint_memories_overlap_instead_of_serializing_PE(self):
        timing = self.d["covariance_memory"]
        g = timing["pseudo_gap_s_upper"]
        spacing = timing["S_observation_window_spacing_s"]
        tobs = timing["observation_window_s_upper"]
        tpe = timing["vector_PE_window_s_upper"]
        word = timing["covariance_memory_window_s_upper"]
        self.assertFalse(timing["S_observation_spacing_uses_vector_PE"])
        self.assertTrue(timing["S_and_vector_PE_memories_overlap_at_endpoint"])
        self.assertTrue(timing["covariance_memory_is_max_not_sum"])
        self.assertEqual(timing["S_observation_window_layout"], "[0,g],[2g,3g],[4g,5g]")
        self.assertGreaterEqual(spacing, 2.0 * g)
        self.assertGreaterEqual(tobs, 5.0 * g)
        self.assertLess(tobs, 5.0 * g + 1e-12)
        self.assertGreaterEqual(word, max(tobs, tpe))
        self.assertLess(word, max(tobs, tpe) + 1e-12)
        # With the deployed 150 ms cadence guard, the translation memory is
        # strictly shorter than the retained 1 s vector-PE recurrence.
        self.assertLess(tobs, tpe)

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
