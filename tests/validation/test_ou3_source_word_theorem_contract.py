import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location(
    "ou3_source_word_theorem_contract",
    ROOT / "tools" / "ou3_source_word_theorem_contract.py",
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class SourceWordTheoremContractTests(unittest.TestCase):
    def test_missing_pe_recurrence_blocks_source_complete_word_claim(self):
        d = mod.build()
        self.assertEqual(mod.validate(d), [])
        self.assertFalse(d["conditional_word_language"]["ready"])
        self.assertFalse(d["source_complete_relative_to_theorem_hypotheses"])
        self.assertFalse(d["continuous_word_enclosed"])
        self.assertFalse(d["nonlinear_word_enclosed"])
        self.assertEqual(d["theorem_promotion"], "NOT_ESTABLISHED")
        self.assertTrue(any("recurrence" in x for x in d["failures"]))

    def test_declared_recurrence_creates_four_S_tileable_language(self):
        d = mod.build(pe_recurrence_window_s=1.0)
        self.assertEqual(mod.validate(d), [])
        self.assertTrue(d["conditional_word_language"]["ready"])
        self.assertTrue(d["source_complete_relative_to_theorem_hypotheses"])
        tr = d["translation_recurrence"]
        self.assertEqual(tr["full_observability_route"], "FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO")
        self.assertEqual(tr["primary_route"], "FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO")
        self.assertEqual(tr["aligned_firing_count"], 4)
        self.assertEqual(tr["primary_state_order"], ["v", "p", "S", "a_w"])
        self.assertEqual(tr["three_firing_integrator_detectability_role"], "Riccati_covariance_upper_sharpening_only")
        self.assertFalse(tr["three_firing_integrator_detectability_is_promotion_fallback"])
        self.assertGreaterEqual(tr["spread_index_q_W"], 1)
        self.assertGreater(tr["spread_selected_spacing_lower_s"], 0.0)
        self.assertGreaterEqual(tr["determinant_spacing_widening_factor_vs_adjacent"], 1)
        self.assertGreaterEqual(d["conditional_word_language"]["word_horizon_lower_s"], tr["minimum_four_firing_window_s"])
        self.assertFalse(d["conditional_word_language"]["one_sample_decrease_required"])
        self.assertTrue(d["conditional_word_language"]["word_endpoint_decrease_required"])

    def test_recurrence_cannot_be_shorter_than_vector_packet_span(self):
        base = mod.build()
        gap_hi = base["vector_persistent_excitation"]["packet_gap_s"][1]
        d = mod.build(pe_recurrence_window_s=0.5 * gap_hi)
        self.assertEqual(mod.validate(d), [])
        self.assertFalse(d["conditional_word_language"]["ready"])
        self.assertTrue(any("shorter" in x for x in d["failures"]))

    def test_required_pe_event_does_not_forbid_other_source_branches(self):
        d = mod.build(pe_recurrence_window_s=1.0)
        branches = d["source_branch_language"]
        self.assertEqual(branches["accelerometer_gate"], ["accepted", "rejected"])
        self.assertEqual(branches["magnetometer_gate"], ["not_due", "accepted", "rejected"])
        self.assertTrue(branches["joint_source_reachability_required"])
        self.assertTrue(branches["cartesian_extrema_products_not_a_valid_word"])
        pe = d["vector_persistent_excitation"]
        self.assertTrue(pe["arbitrary_rejections_between_required_pe_events_allowed"])
        self.assertTrue(pe["two_consecutive_accepted_magnetic_packets_required"])
        self.assertTrue(pe["accelerometer_required_at_both_vector_times"])

    def test_hybrid_events_are_not_multiplied_into_same_mode_words(self):
        d = mod.build(pe_recurrence_window_s=1.0)
        scope = d["normal_live_scope"]
        self.assertTrue(scope["same_mode_only"])
        self.assertFalse(scope["dimension_change_multiplied_as_square_word"])
        self.assertIn("held_to_active", scope["hybrid_transitions_separate"])
        self.assertIn("tilt_reset", scope["hybrid_transitions_separate"])

    def test_information_enclosure_contract_requires_the_same_word_language(self):
        text = (ROOT / "tools" / "ou3_information_enclosure_contract.py").read_text()
        self.assertIn("tools/ou3_source_word_theorem_contract.py", text)
        self.assertIn("OU3_CONDITIONAL_SOURCE_COMPLETE_NORMAL_LIVE_WORD_LANGUAGE", text)
        self.assertIn('"pe_recurrence_window_s"', text)
        self.assertIn("four spread-selected S", text)
        self.assertIn("single favorable vector pair", text)

    def test_contract_uses_no_replay_or_observed_extrema(self):
        text = (ROOT / "tools" / "ou3_source_word_theorem_contract.py").read_text()
        for forbidden in ("ou3_exact_replay", "path_metrics", "neighborhood_radius_search", "observed_min", "replay_min", "np.quantile"):
            self.assertNotIn(forbidden, text)


if __name__ == "__main__":
    unittest.main()
