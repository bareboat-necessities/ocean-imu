from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "stability"))
import ou3_brmm_physical_wave_condition as PHYS
import ou3_brmm_physical_wave_source as CERT


class PhysicalWaveConditionTests(unittest.TestCase):
    def test_physics_is_primary_and_certificate_methods_are_not_definition(self):
        d = PHYS.physical_condition()
        self.assertTrue(d["physics_is_primary"])
        self.assertTrue(d["certificate_methods_are_sufficient_not_definitional"])
        self.assertEqual(d["wave_coordinate"], "p_wave=x_CoG-x_local_equilibrium")
        self.assertFalse(d["D_S_may_be_chosen_from_P4_working_radius"])
        self.assertFalse(d["zero_mean_alone_is_sufficient"])
        self.assertFalse(d["power_spectrum_alone_is_sufficient"])

    def test_supported_representation_proves_condition_but_does_not_define_it(self):
        c = CERT.spectral_certificate((CERT.SpectralBand(F(1, 10), 2, 3),))
        q = PHYS.qualify_with_certificate(c)
        self.assertTrue(q["bounded_primitive_proved"])
        self.assertFalse(q["certificate_is_definition_of_wave_motion"])
        self.assertEqual(q["derived_D_S_m_s"], "60")

    def test_pr515_constant_position_witness_is_physically_inadmissible(self):
        c = CERT.spectral_certificate((CERT.SpectralBand(1, 2, 3),))
        self.assertFalse(PHYS.constant_history_admitted((F(1, 8), 0, 0), c))
        self.assertFalse(PHYS.constant_history_admitted((F(-1, 8), 0, 0), c))
        self.assertTrue(PHYS.constant_history_admitted((0, 0, 0), c))
        d = PHYS.build()
        self.assertFalse(d["constant_nonzero_position_zero_velocity_history_admitted"])
        self.assertEqual(d["historical_old_finite_window_classification"], "B")
        self.assertEqual(d["physical_source_specification_omission_classification"], "E")
        self.assertFalse(d["shipping_filter_instability_claimed"])

    def test_complete_family_numeric_bound_remains_fail_closed(self):
        d = PHYS.build()
        self.assertFalse(d["physical_D_S_numeric_qualification_closed_for_complete_family"])
        self.assertFalse(d["P4_PASS"])
        self.assertFalse(d["P5_MAY_START"])
        self.assertEqual(d["P3_delta"], 1e-18)
        self.assertEqual(PHYS.validate(d), [])


if __name__ == "__main__":
    unittest.main()
