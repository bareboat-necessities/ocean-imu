from fractions import Fraction as F
import unittest
from tools.stability.ou3_alt_contraction.bias_supply_storage import (
    coefficients, bias_coordinate_squared_bound,
)

class BiasSupplyStorageTests(unittest.TestCase):
    def test_input_cross_terms_cannot_be_discarded(self):
        # Exact SPD counterexample to V(u,s)>=V(u,0); not a shipping install.
        off = F(-99,100)
        self.assertGreater(1-off*off,0)
        full = 2+2*off
        self.assertEqual(full,F(1,50))
        self.assertLess(full,1)
        # Rational squared form of reverse triangle: |E-V-B|²<=4VB.
        self.assertLessEqual((1-full-1)**2,4*full)

    def test_margin_is_retained_with_bias_supply(self):
        result = coefficients(F(989,1000),F(1,1000),2,3,F(1,100))
        self.assertEqual(result['rho'],F(989901,1000000))
        self.assertLess(result['rho'],1)
        self.assertEqual(result['C'],101*F(249,50)**2)
        self.assertFalse(result['native_premises_proved'])

    def test_cross_term_gain_is_not_free(self):
        without = coefficients(F(9,10),0,0,1,F(1,100))
        with_bias = coefficients(F(9,10),0,2,1,F(1,100))
        self.assertEqual(without['rho'],with_bias['rho'])
        self.assertGreater(with_bias['C'],without['C'])

    def test_bias_bound_keeps_practical_floor(self):
        self.assertEqual(bias_coordinate_squared_bound(F(2,5),0),F(4,25))
        self.assertEqual(bias_coordinate_squared_bound(F(2,5),F(1,10)),F(13,50))
        for invalid in (0.1,True):
            with self.assertRaises(TypeError):
                coefficients(invalid,0,0,0,1)
        with self.assertRaises(ValueError):
            coefficients(1,0,0,0,0)

if __name__ == '__main__':
    unittest.main()
