from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_sea3_riccati_metric_p3 as mod  # noqa: E402


class Sea3RsInnovationP3Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_architecture_uses_rs_and_sea3_as_primary_strictness(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(
            self.d["canonical_P3_architecture"],
            "SEA3_RS_INNOVATION_DISSIPATION_WORD",
        )
        self.assertTrue(self.d["R_S_is_primary_translation_correction_mechanism"])
        self.assertTrue(self.d["pseudo_update_recurrence_is_primary_word_structure"])
        self.assertTrue(self.d["SEA3_physical_parameter_coupling_consumed"])
        self.assertTrue(self.d["same_operating_point_tau_sigma_RS_TS_required"])
        self.assertTrue(self.d["independent_source_extrema_product_forbidden"])
        self.assertEqual(self.d["strictness_location"], "RECURRENT_SEA3_MEASUREMENT_WORD")

    def test_exact_innovation_identity_is_canonical(self):
        self.assertTrue(self.d["exact_measurement_dissipation_identity_consumed"])
        self.assertTrue(self.d["batch_innovation_information_identity_consumed"])
        self.assertTrue(self.d["process_UCC_used_as_metric_lower_not_primary_strictness"])
        self.assertIn("V_minus - V_plus", self.d["innovation_identity"]["identity"])
        self.assertIn("D_W", self.d["batch_identity"]["correction_information"])

    def test_dead_end_routes_cannot_reenter(self):
        self.assertFalse(self.d["source_history_graph_consumed"])
        self.assertFalse(self.d["predecessor_path_enumeration_consumed"])
        self.assertFalse(self.d["old_P2_800_state_graph_consumed"])
        self.assertFalse(self.d["one_sample_strict_Riccati_margin_consumed"])
        self.assertFalse(self.d["commit_aligned_source_word_consumed"])
        self.assertFalse(self.d["per_sample_SPD_lower_required"])
        self.assertFalse(self.d["selected_process_mode_strictness_used"])
        self.assertFalse(self.d["determinant_trace_scalarization_used"])
        self.assertFalse(self.d["scalar_information_beta_used"])
        self.assertEqual(self.d["useful_gate"], 1e-18)

    def test_gate_fails_closed_until_dw_lw_matrix_comparison(self):
        self.assertFalse(self.d["P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED"])
        self.assertFalse(self.d["P3_UCC_METRIC_LOWER_CLOSED"])
        self.assertFalse(self.d["P3_FULL_MATRIX_COMPARISON_CLOSED"])
        self.assertFalse(self.d["P3_CANONICAL_PASS"])
        self.assertFalse(self.d["P4_MAY_CONSUME_P3"])
        self.assertTrue(self.d["P3_CANONICAL_FAIL_REASONS"])
        for mode in ("H", "A"):
            self.assertFalse(self.d["modes"][mode]["pass"])
            self.assertEqual(
                self.d["modes"][mode]["relative_Riccati_injection_margin_lower"], 0.0
            )


if __name__ == "__main__":
    unittest.main()
