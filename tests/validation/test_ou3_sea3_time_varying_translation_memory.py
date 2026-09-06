#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_time_varying_translation_memory as CERT


class Sea3TimeVaryingTranslationMemoryTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = CERT.build()
        cls.failures = CERT.validate(cls.d)

    def test_certificate_validates(self):
        self.assertEqual(self.failures, [])
        self.assertEqual(self.d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.d["full_4x4_time_varying_translation_memory_closed"])

    def test_no_frozen_tuner_or_history_source(self):
        self.assertFalse(self.d["constant_tau_over_memory_assumed"])
        self.assertTrue(self.d["time_varying_committed_tau_sigma_RS_allowed"])
        self.assertFalse(self.d["source_history_graph_consumed"])
        self.assertFalse(self.d["predecessor_path_enumeration_consumed"])
        self.assertFalse(self.d["source_family_replaced"])
        self.assertFalse(self.d["trajectory_replay_used"])

    def test_shipping_commit_intervals_are_used(self):
        g = self.d["commit_geometry"]
        self.assertGreaterEqual(g["complete_constant_tune_intervals_lower"], 26)
        self.assertEqual(g["intervals_retained"], 26)
        self.assertGreaterEqual(g["min_constant_commit_interval_samples_conservative"], 20)
        self.assertLessEqual(g["max_commit_interval_samples_certified"], 22)
        self.assertTrue(all(self.d["shipping_commit_parity"].values()))

    def test_measurements_are_interleaved_sample_by_sample(self):
        m = self.d["measurement_lower_recursion"]
        self.assertTrue(m["optimal_posterior_lower_for_any_Joseph_gain"])
        self.assertTrue(m["S_zero_assumed_due_every_sample"])
        self.assertTrue(m["accelerometer_cross_block_factor_three_paid"])
        self.assertTrue(m["measurements_executed_sample_by_sample"])
        self.assertFalse(m["measurements_moved_to_word_endpoint"])

    def test_full_matrix_induction_and_suffix_close(self):
        ind = self.d["candidate_induction"]
        tail = self.d["terminal_suffix"]
        self.assertEqual(ind["failures"], [])
        self.assertEqual(tail["failures"], [])
        self.assertGreater(ind["certified_leaves"], 0)
        self.assertGreater(ind["worst_conditioned_LDLT_pivot_lower"], 0.0)
        self.assertTrue(tail["positive_suffix_preserves_final_lower"])
        self.assertGreater(tail["worst_conditioned_LDLT_pivot_lower"], 0.0)
        M = self.d["word_endpoint_translation_process_measurement_noise_covariance_lower"]
        self.assertEqual(len(M), 4)
        self.assertTrue(all(len(row) == 4 for row in M))

    def test_process_is_shipping_equivalent_and_not_a_p3_promotion(self):
        self.assertTrue(self.d["process"]["exact_shipping_integrated_OU_Q_consumed"])
        self.assertFalse(self.d["P3_promoted"])


if __name__ == "__main__":
    unittest.main()
