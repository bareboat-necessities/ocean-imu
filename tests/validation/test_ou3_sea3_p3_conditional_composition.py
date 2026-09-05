from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_p3_conditional_composition as mod  # noqa: E402


class Sea3ConditionalP3CompositionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()

    def test_conditional_p3_and_physical_left_inclusion_are_separate(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(
            self.d["conditional_P3_quantifier"],
            "FOR_EVERY_ADMITTED_COMPLETE_SEA3_NORMAL_LIVE_WORD",
        )
        self.assertTrue(self.d["theorem_conditional_on_admitted_complete_SEA3_word"])
        self.assertTrue(self.d["global_physical_deployment_left_inclusion_is_separate_obligation"])
        self.assertFalse(self.d["global_physical_left_inclusion_required_before_conditional_P3_math"])
        self.assertTrue(self.d["global_physical_left_inclusion_required_before_full_deployment_theorem"])
        self.assertFalse(self.d["hard_shaping_physical_left_inclusion_currently_closed"])

    def test_complete_word_uses_four_s_information_without_source_replacement(self):
        self.assertTrue(self.d["four_S_information_component_consumed"])
        self.assertTrue(self.d["actual_applied_SpectralMSE_R_S_consumed"])
        self.assertTrue(self.d["all_due_S_updates_remain_in_complete_word"])
        self.assertFalse(self.d["four_S_selected_events_replace_complete_word"])
        self.assertFalse(self.d["source_family_replaced"])

    def test_open_numeric_obligations_are_not_falsely_promoted(self):
        self.assertTrue(self.d["H18_full_information_matrix_lower_closed"])
        self.assertGreater(self.d["H18_information_lambda_min_lower"], 0.0)
        self.assertTrue(self.d["A21_finite_bias_correlation_route_consumed"])
        self.assertFalse(self.d["A21_uses_eta9_packet_shortcut"])
        self.assertFalse(self.d["H18_prior_free_completion_closed"])
        self.assertFalse(self.d["A21_detectability_completion_closed"])
        self.assertFalse(self.d["P3_CANONICAL_PASS"])
        self.assertFalse(self.d["P4_MAY_CONSUME_P3"])

    def test_event_algebra_includes_reset_and_preserves_first_margin(self):
        self.assertTrue(self.d["exact_complete_word_event_algebra_consumed"])
        self.assertTrue(self.d["event_algebra_covers_immediate_left_error_reset"])
        self.assertTrue(self.d["event_algebra_preserves_first_established_full_matrix_margin"])


if __name__ == "__main__":
    unittest.main()
