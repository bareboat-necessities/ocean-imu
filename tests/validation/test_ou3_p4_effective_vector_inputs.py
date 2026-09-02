from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_effective_vector_inputs as EFFECTIVE
import ou3_p4_signed_joseph_sector_inputs as COMPOSE


class P4EffectiveVectorInputsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.effective = EFFECTIVE.build()

    def test_effective_input_certificate_covers_declared_08rad_sector(self):
        self.assertEqual(EFFECTIVE.validate(self.effective), [])
        self.assertTrue(self.effective["pass"])
        self.assertLess(self.effective["cayley_radius_upper"], 1.0)
        self.assertLess(
            self.effective["mag_effective_coordinate_tangent_defect_factor_upper"],
            1.0,
        )

    def test_shipping_range_and_nullspace_identities_are_active(self):
        self.assertTrue(self.effective["mag_radial_residual_gain_null_exact"])
        self.assertTrue(self.effective["mag_effective_coordinate_identity_exact"])
        self.assertTrue(self.effective["acc_eta_in_aw_measurement_range_exact"])
        self.assertTrue(self.effective["acc_effective_aw_input_isometry_exact"])
        self.assertEqual(
            self.effective["accelerometer_bias_standalone_nonlinear_penalty"],
            0.0,
        )

    def test_active_route_forbids_old_metric_penalties(self):
        self.assertTrue(
            self.effective["standalone_vector_eta_penalty_retired_from_active_word_route"]
        )
        self.assertFalse(
            self.effective["condition_number_conversion_inserted_between_P3_and_P4"]
        )
        self.assertFalse(self.effective["P4_USABLE_CERTIFICATE_PROMOTED"])

    def test_composition_contract_retains_same_information_metric(self):
        d = COMPOSE.build()
        self.assertEqual(COMPOSE.validate(d), [])
        self.assertTrue(d["P4_COMPOSITION_PREREQUISITES_ESTABLISHED"])
        self.assertTrue(d["same_information_metric_retained"])
        self.assertFalse(d["standalone_vector_eta_penalty_active"])
        self.assertFalse(
            d["condition_number_conversion_inserted_between_P3_and_P4"]
        )
        self.assertTrue(d["full_source_correlated_word_form_remaining"])
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"])
        self.assertFalse(d["P5_FINITE_CAPTURE_ESTABLISHED"])
        for mode in ("H", "A"):
            m = d["modes"][mode]
            self.assertGreater(m["P3_relative_Riccati_injection_margin_lower"], 0.0)
            self.assertEqual(m["P3_prefix_information_gain_upper"], 1.0)
            self.assertFalse(m["standalone_vector_eta_penalty_active"])
            self.assertFalse(
                m["condition_number_conversion_inserted_between_P3_and_P4"]
            )
            self.assertIsNone(m["signed_word_generalized_margin_lower"])
            self.assertIsNone(m["rho_full_nonlinear_word_upper"])


if __name__ == "__main__":
    unittest.main()
