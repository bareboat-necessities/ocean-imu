import unittest
from fractions import Fraction as F
import mpmath as mp
from tools.stability.ou3_alt_contraction.driven_storage_diagnostic import (
    necessary_supply, source_certificate,
)

class DrivenStorageTests(unittest.TestCase):
    def test_nonrelaxing_bias_member_is_separate(self):
        cert = source_certificate('0.5', 2)
        self.assertTrue(cert['bias2_phi_one'])
        self.assertIsNone(cert['bias_root_tau_s'])
        self.assertLessEqual(F(cert['bias_derivative_bound']), F(2,10000))
        self.assertEqual(F(cert['bias_driver_bound']), F(cert['bias_derivative_bound'])/200)

    def test_all_time_primitive_and_supply(self):
        for family in range(3):
            for omega in ('0.5','1'):
                cert = source_certificate(omega, family)
                w = F(omega)
                # Endpoint-difference primitive bound, not a wordwise reset.
                self.assertEqual(F(cert['source_qualification']['derived_D_S_m_s']), F(1,2)/w)
                B = F(cert['bias_magnitude_bound'])
                expected = (1+w*w+w**4)/16 + 1/(4*w*w) + B*B
                self.assertEqual(F(cert['physical_supply_squared_bound']), expected)

    def test_forced_growth_is_not_unforced_contraction(self):
        self.assertEqual(necessary_supply(mp.mpf(2),mp.mpf(5),mp.mpf('.5'),mp.mpf(4)),1)
        self.assertEqual(necessary_supply(mp.mpf(2),mp.mpf('.5'),mp.mpf('.5'),mp.mpf(4)),0)
        with self.assertRaises(ValueError):
            necessary_supply(1,2,1,1)

    def test_no_unqualified_members(self):
        for omega,family in [('0',0),('0.6',0),('1',3)]:
            with self.assertRaises(ValueError):
                source_certificate(omega,family)

if __name__ == '__main__':
    unittest.main()
