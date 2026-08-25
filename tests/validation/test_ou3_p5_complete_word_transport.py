import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_complete_word_transport as W


class Ou3P5CompleteWordTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = W.build()

    def test_every_shipping_operation_class_has_exact_transport_semantics(self):
        d = self.d
        self.assertEqual(W.validate(d), [])
        self.assertEqual(d["P5_COMPLETE_WORD_TRANSPORT_ALGEBRA_CERTIFICATE"], "PASS")
        self.assertTrue(d["all_source_operation_classes_bound_to_transport_calculus"])
        self.assertEqual(
            [x["operation"] for x in d["operation_transport_calculus"]],
            d["source_operation_order"],
        )

    def test_joseph_reset_and_S_cross_gain_are_not_replaced_by_norm_shortcuts(self):
        d = self.d
        self.assertIn("eta^T R^-1 eta", d["accepted_measurement_identity"])
        self.assertIn("Ge^-1 rho", d["reset_metric_identity"])
        self.assertFalse(d["reset_condition_number_multiplier_used"])
        self.assertTrue(d["S_zero_nonlinear_measurement_eta_exact_zero"])
        self.assertTrue(d["full_S_to_attitude_gain_retained"])
        self.assertTrue(d["sequential_immediate_quaternion_resets_retained"])

    def test_gauged_word_consumes_widened_first_S_and_names_later_budget(self):
        g = self.d["gauged_H"]
        self.assertEqual(g["finite_angle_information_geometry"], "PASS")
        self.assertEqual(g["exact_correction_transport_algebra"], "PASS")
        self.assertEqual(g["first_due_S_exact_prefix"], "PASS_WIDENED_CHART")
        self.assertFalse(g["diagnostic_q_lt_1_is_promotion_gate"])
        self.assertGreater(g["widened_prefix_antipodal_margin_lower"], 0.0)
        self.assertGreater(g["widened_prefix_vector_information_vs_goLive_metric_lower"], 0.0)
        self.assertEqual(
            g["first_unclosed_numerical_obligation"],
            "COMPLETE_WORD_ETA_RESET_INFORMATION_BUDGET_NOT_CERTIFIED",
        )
        self.assertEqual(g["complete_word_numerical_status"], "NOT_ESTABLISHED")

    def test_gravity_quotient_charges_axial_bias_as_finite_input(self):
        q = self.d["gravity_quotient_H"]
        self.assertEqual(q["reduced_detectability"], "PASS")
        self.assertFalse(q["strict_contraction_of_axial_bias_requested"])
        self.assertGreater(q["one_word_axial_bias_full_attitude_coordinate_input_norm_upper_rad"], 0.0)
        self.assertLess(q["one_word_axial_bias_full_attitude_coordinate_input_norm_upper_rad"], 0.02)
        self.assertIn("gravity_tilt_2d", q["strict_coordinates"])
        self.assertNotIn("gyro_bias_parallel", q["strict_coordinates"])

    def test_numerical_word_certificates_remain_fail_closed(self):
        d = self.d
        self.assertEqual(d["P5_GAUGED_COMPLETE_WORD_NUMERICAL_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertEqual(d["P5_GRAVITY_QUOTIENT_COMPLETE_WORD_NUMERICAL_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])


if __name__ == "__main__":
    unittest.main()
