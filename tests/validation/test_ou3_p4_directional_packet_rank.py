from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_directional_packet_rank as RANK


class DirectionalPacketRankTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = RANK.build()

    def test_structural_rank_certificate_validates(self):
        self.assertEqual(RANK.validate(self.d), [])
        self.assertTrue(self.d["directional_PSD_operation_credit_required"])
        self.assertTrue(self.d["word_level_directional_accumulation_required"])
        self.assertFalse(self.d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"])
        self.assertFalse(self.d["P4_USABLE_CERTIFICATE_PROMOTED"])

    def test_same_sample_vector_packet_is_exactly_rank_five(self):
        s = self.d["measurement_structure"]
        self.assertEqual(s["accelerometer_rank_exact"], 3)
        self.assertEqual(s["magnetometer_rank_exact"], 2)
        self.assertEqual(s["stacked_vector_packet_rank_exact"], 5)
        self.assertFalse(
            self.d["instantaneous_full_state_measurement_information_positive_definite_possible"])
        self.assertFalse(
            self.d["instantaneous_positive_scalar_full_state_packet_margin_is_valid_target"])

    def test_exact_null_witness_is_recorded(self):
        w = self.d["exact_vector_packet_null_witness"]
        self.assertTrue(w["nonzero_for_nonzero_m"])
        self.assertEqual(w["delta_theta"], "alpha * m")
        self.assertIn("[f]_x", w["delta_a_w"])
        self.assertIn("= 0", w["magnetometer_residual"])
        self.assertIn("= 0", w["accelerometer_residual"])

    def test_H_and_A_nullities_force_word_accumulation(self):
        h = self.d["modes"]["H"]
        a = self.d["modes"]["A"]
        self.assertEqual(h["stacked_vector_packet_nullity_exact_on_active_block"], 1)
        self.assertEqual(a["stacked_vector_packet_nullity_exact_on_active_block"], 4)
        self.assertEqual(h["same_sample_full_state_nullity_without_S_lower"], 13)
        self.assertEqual(a["same_sample_full_state_nullity_without_S_lower"], 16)
        self.assertEqual(h["same_sample_full_state_nullity_with_due_S_lower"], 10)
        self.assertEqual(a["same_sample_full_state_nullity_with_due_S_lower"], 13)
        self.assertTrue(self.d["P4_must_transport_directional_nullspaces_through_prediction"])

    def test_no_domain_or_filter_shortcut(self):
        self.assertFalse(self.d["source_replay_used"])
        self.assertFalse(self.d["filter_changed"])
        self.assertFalse(self.d["declared_domain_changed"])
        self.assertFalse(self.d["aw_sigma_consistency_assumption_used"])


if __name__ == "__main__":
    unittest.main()
