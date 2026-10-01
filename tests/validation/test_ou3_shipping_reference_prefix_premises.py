"""Exact arithmetic checks for the shipping reference/prefix premise audit."""
import math
import unittest


class ShippingReferencePrefixPremiseTests(unittest.TestCase):
    def test_continuous_safety_bias_gate_does_not_preserve_physical_horizontal_floor(self):
        # Declared physical max plus hard-iron and measurement residual envelopes.
        measured_scale_max = 75.0 + 5.0 + 2.0
        accepted_bias_max = 0.35 * measured_scale_max
        self.assertGreater(accepted_bias_max, 15.0)

    def test_six_degree_capture_does_not_imply_field_exclusion_storage_radius(self):
        g = 9.80665
        physical_charge = 0.26
        kappa_aw = 4.06
        r_max = (g * math.sin(math.radians(10.0)) - physical_charge) / kappa_aw
        handoff_storage_radius_at_six_deg = math.radians(6.0) / 0.035
        self.assertLess(r_max, 0.356)
        self.assertGreater(handoff_storage_radius_at_six_deg, 2.99)
        self.assertGreater(handoff_storage_radius_at_six_deg, 8.0 * r_max)

    def test_component_aw_tube_needed_is_far_weaker_than_generic_quarter_radius(self):
        g = 9.80665
        physical_charge = 0.26
        eps_aw_max = g * math.sin(math.radians(10.0)) - physical_charge
        generic_quarter_radius_aw_charge = 4.06 * 0.25
        self.assertGreater(eps_aw_max, 1.4429)
        self.assertLess(generic_quarter_radius_aw_charge, 1.016)
        self.assertGreater(eps_aw_max, generic_quarter_radius_aw_charge)

    def test_reference_lipschitz_constant_from_convex_rotation_average(self):
        # ||Abar||<=1 makes the horizontal norm and z component each 1-Lipschitz.
        # The canonical (h,0,z) reference is therefore sqrt(2)-Lipschitz.
        self.assertLess(math.sqrt(2.0), 1.415)


if __name__ == "__main__":
    unittest.main()
