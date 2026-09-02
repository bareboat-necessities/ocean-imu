from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_effective_word_transport as WORD


class P4EffectiveWordTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = WORD.build()

    def test_transport_semantics_are_source_complete_and_fail_closed(self):
        self.assertEqual(WORD.validate(self.d), [])
        self.assertTrue(self.d["P4_WORD_TRANSPORT_SEMANTICS_ESTABLISHED"])
        self.assertFalse(self.d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"])
        self.assertFalse(self.d["P5_FINITE_CAPTURE_ESTABLISHED"])
        self.assertTrue(self.d["joint_source_reachability_required"])
        self.assertFalse(self.d["cartesian_extrema_products_valid"])

    def test_exact_shipping_operation_order_and_sequential_resets_are_retained(self):
        self.assertEqual(self.d["shipping_operation_order"], WORD.EXPECTED_ORDER)
        self.assertEqual(
            [row["operation"] for row in self.d["operation_transport_calculus"]],
            WORD.EXPECTED_ORDER,
        )
        acc = next(
            row for row in self.d["operation_transport_calculus"]
            if row["operation"].startswith("accelerometer_correction")
        )
        mag = next(
            row for row in self.d["operation_transport_calculus"]
            if row["operation"].startswith("asynchronous_magnetometer")
        )
        self.assertTrue(acc["effective_input_isometry"])
        self.assertFalse(acc["standalone_eta_penalty_active"])
        self.assertEqual(mag["radial_residual_gain_action"], "EXACTLY_ZERO")
        self.assertFalse(mag["standalone_eta_penalty_active"])

    def test_vector_packet_rank_requires_directional_word_accumulation(self):
        rank = self.d["measurement_directional_rank"]
        self.assertEqual(rank["accelerometer_rank_exact"], 3)
        self.assertEqual(rank["magnetometer_rank_exact"], 2)
        self.assertEqual(rank["stacked_vector_packet_rank_exact"], 5)
        self.assertEqual(rank["H_vector_packet_nullity_exact"], 1)
        self.assertEqual(rank["A_vector_packet_nullity_exact"], 4)
        self.assertFalse(
            rank["instantaneous_positive_scalar_full_state_packet_margin_possible"]
        )
        self.assertTrue(rank["directional_PSD_word_accumulation_required"])
        self.assertFalse(self.d["per_packet_scalarization_allowed"])

    def test_P3_to_P4_metric_shortcuts_cannot_return(self):
        self.assertTrue(self.d["S_timing_consumed_by_linear_P3"])
        self.assertFalse(self.d["standalone_vector_eta_penalty_active"])
        self.assertFalse(
            self.d["condition_number_conversion_inserted_between_P3_and_P4"]
        )
        self.assertFalse(self.d["reset_condition_number_multiplier_used"])
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            self.assertEqual(m["same_information_metric"], "M_i=s_m Sigma_i^-1")
            self.assertTrue(m["full_attitude_linear_cross_terms_retained"])
            self.assertGreater(m["P3_relative_Riccati_injection_margin_lower"], 0.0)
            self.assertEqual(m["P3_prefix_information_gain_upper"], 1.0)
            self.assertFalse(m["word_directional_forms_numerically_accumulated_here"])
            self.assertIsNone(m["strict_generalized_word_margin_lower"])
            self.assertIsNone(m["rho_full_nonlinear_word_upper"])
            self.assertFalse(m["P4_PROMOTED"])


if __name__ == "__main__":
    unittest.main()
