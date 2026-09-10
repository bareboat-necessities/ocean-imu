import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_translation_observation_kernel as K


class TranslationObservationKernelTests(unittest.TestCase):
    def test_exact_kernel_collapses_only_with_both_primitive_bounds(self):
        d = K.build()
        self.assertEqual(K.validate(d), [])
        self.assertEqual(d["translation_zero_output_kernel_dimension_before_indefinite_primitives_per_axis"], 2)
        self.assertEqual(d["translation_zero_output_kernel_dimension_after_bounded_position_per_axis"], 1)
        self.assertEqual(d["translation_zero_output_kernel_dimension_after_centered_S_recurrence_per_axis"], 0)
        self.assertTrue(d["all_three_axes_kernel_collapses"])
        self.assertTrue(d["only_finiteness_of_D_S_used"])
        self.assertIsNone(d["D_S_numeric_value_used"])

    def test_exact_primitive_polynomials_are_fixed(self):
        d = K.build()["difference_basis"]
        self.assertEqual(d["velocity_constant_c_v"]["delta_v"], ["1"])
        self.assertEqual(d["velocity_constant_c_v"]["delta_p"], ["0", "1"])
        self.assertEqual(d["velocity_constant_c_v"]["delta_S_L"], ["0", "0", "1/2"])
        self.assertEqual(d["position_constant_c_p"]["delta_v"], ["0"])
        self.assertEqual(d["position_constant_c_p"]["delta_p"], ["1"])
        self.assertEqual(d["position_constant_c_p"]["delta_S_L"], ["0", "1"])

    def test_no_legacy_entry_or_coordinate_shortcut(self):
        d = K.build()
        self.assertFalse(d["absolute_S_origin_needed_for_kernel_collapse"])
        self.assertFalse(d["legacy_300_m_s_entry_ball_used"])
        self.assertFalse(d["position_reanchoring_used"])
        self.assertFalse(d["wordwise_S_rezero_used"])
        self.assertFalse(d["covariance_membership_used"])
        self.assertFalse(d["shipping_filter_changed"])

    def test_false_p4_promotion_is_rejected(self):
        d = K.build()
        d["P4_PASS"] = True
        self.assertIn("P4_PASS not false", K.validate(d))

    def test_mutated_quadratic_primitive_is_rejected(self):
        d = K.build()
        d["difference_basis"]["velocity_constant_c_v"]["delta_S_L"] = ["0", "0", "1"]
        self.assertIn("c_v S polynomial changed", K.validate(d))


if __name__ == "__main__":
    unittest.main()
