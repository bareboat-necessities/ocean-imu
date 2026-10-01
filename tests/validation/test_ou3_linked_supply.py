import unittest
from fractions import Fraction as F

from tools.stability.ou3_theorem.linked_supply import (
    compose_defects, linked_supply, completed_square_check, retained_entry_budget)


class LinkedSupplyTests(unittest.TestCase):
    def test_scalar_sharp_constant_and_equality(self):
        cert = completed_square_check([[1]], [[1]], [[F(1, 2)]], [[1]], F(1, 4), [[1]])
        self.assertEqual(cert['chi'], F(3, 2))
        self.assertEqual(cert['center'], [[F(1)]])
        self.assertEqual(cert['terminal_storage'], F(9, 4))
        self.assertFalse(cert['source_uniform_verified'])

    def test_correlated_metrics_identity(self):
        j0, jn = [[3, 1], [1, 2]], [[2, F(1, 3)], [F(1, 3), 1]]
        m, b = [[F(1, 3), F(1, 8)], [0, F(1, 2)]], [[F(2, 7)], [F(-3, 11)]]
        for e in ([[0], [0]], [[3], [-5]], [[F(1, 29)], [F(-9, 17)]]):
            self.assertTrue(completed_square_check(j0, jn, m, b, F(1, 5), e)['identity_verified'])

    def test_forcing_direction_not_only_norm_matters(self):
        j, m = [[1, 0], [0, 1]], [[F(99, 100), 0], [0, F(1, 2)]]
        weak = linked_supply(j, j, m, [[1], [0]], F(1, 100))['chi']
        strong = linked_supply(j, j, m, [[0], [1]], F(1, 100))['chi']
        self.assertEqual(weak, F(100))
        self.assertEqual(strong, F(99, 74))

    def test_signed_transport_cancellation(self):
        m, b = compose_defects([{'A': [[1]], 'd': [[1]]},
                                {'A': [[2]], 'd': [[-2]]}], 1)
        self.assertEqual(m, [[F(2)]])
        self.assertEqual(b, [[F(0)]])

    def test_singular_or_indefinite_relative_loss_is_rejected(self):
        for m in ([[1]], [[2]]):
            with self.assertRaises((ValueError, ArithmeticError, ZeroDivisionError)):
                linked_supply([[1]], [[1]], m, [[0]], F(1, 2))

    def test_nonzero_actual_innovation_can_have_zero_homogeneous_loss(self):
        # H=(1,0), P=I, R=1, homogeneous error e=(0,1).
        # H e=0, but choose actual measured innovation r=2: K r=(1,0).
        h, e, r = (F(1), F(0)), (F(0), F(1)), F(2)
        homogeneous_loss = (h[0]*e[0]+h[1]*e[1])**2/F(2)
        self.assertEqual(homogeneous_loss, 0)
        self.assertNotEqual(F(1, 2)*r, 0)

    def test_covariance_lower_bound_does_not_make_storage_sublevel_compact(self):
        for n in (1, 10, 1000):
            p, e = F(n*n), F(n)
            self.assertGreaterEqual(p, 1)
            self.assertEqual(e*e/p, 1)

    def test_boundary_equality_is_not_finite_entry(self):
        result = retained_entry_budget(F(1, 2), F(1, 2), 1)
        self.assertTrue(result['root_self_map_budget'])
        self.assertFalse(result['strict_entry_budget'])
        v = F(2)
        for _ in range(40):
            v = (v+1)/2
            self.assertGreater(v, 1)

    def test_strict_entry_budget(self):
        result = retained_entry_budget(F(1, 2), F(1, 4), 1)
        self.assertTrue(result['strict_entry_budget'])
        self.assertFalse(result['prefix_retention_verified'])

    def test_invalid_arguments(self):
        with self.assertRaises(ValueError):
            linked_supply([[1]], [[1]], [[F(1, 2)]], [[1]], 0)
        with self.assertRaises(ValueError):
            linked_supply([[1]], [[1]], [[F(1, 2)]], [[1, 2]], F(1, 4))
        with self.assertRaises(ValueError):
            compose_defects([], 0)


if __name__ == '__main__':
    unittest.main()
