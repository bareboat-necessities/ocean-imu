from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_h18_information_composition as mod  # noqa: E402


class Sea3H18InformationCompositionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_full_H18_information_clears_useful_gate(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(self.d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        c = self.d["triangular_information_composition"]
        self.assertTrue(c["full_18x18_matrix_information_lower_closed"])
        self.assertTrue(math.isfinite(c["D_H18_lambda_min_lower"]))
        self.assertGreaterEqual(c["D_H18_lambda_min_lower"], mod.USEFUL_GATE)
        self.assertTrue(self.d["H18_information_useful_gate_pass"])

    def test_corrected_tight_four_s_is_consumed(self):
        self.assertTrue(self.d["tight_four_S_measurement_covariance_structure_consumed"])
        self.assertTrue(self.d["four_S_process_cross_record_trace_bound_retained"])
        self.assertIn("TIGHT_COVARIANCE", self.d["translation_information_source"])

    def test_four_s_is_component_not_replacement_word(self):
        self.assertTrue(self.d["component_of_complete_SEA3_full_word"])
        self.assertFalse(self.d["source_family_replaced"])
        self.assertFalse(self.d["selected_four_S_events_replace_complete_word"])
        self.assertTrue(self.d["all_due_S_updates_remain_in_literal_word"])
        self.assertTrue(self.d["actual_applied_SpectralMSE_R_S_consumed"])

    def test_no_scalar_or_blockwise_promotion_shortcut(self):
        self.assertFalse(self.d["determinant_trace_scalarization_of_18x18_matrix_used"])
        self.assertFalse(self.d["blockwise_minimum_ratio_used"])
        self.assertFalse(self.d["scalar_information_beta_used"])
        self.assertFalse(self.d["P3_promoted"])


if __name__ == "__main__":
    unittest.main()
