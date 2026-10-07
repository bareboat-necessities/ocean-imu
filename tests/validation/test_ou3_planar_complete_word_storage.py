import unittest
from fractions import Fraction as F

from tools.stability.ou3_theorem.matrix_certificates import add, identity, matmul, transpose
from tools.stability.ou3_theorem.planar_complete_word_storage import (
    certificate, covariance_residual_charge, optimal_correction_storage, storage_identity,
    storage_inequality, telescoped_gap,
)


def energy(x, J):
    return matmul(transpose(x), matmul(J, x))[0][0]


class CompleteWordStorageTests(unittest.TestCase):
    def test_noncommuting_word_telescopes_and_retains_shared_gauge(self):
        factors = [[[F(1, 2), F(1, 5)], [F(0), F(1, 3)]],
                   [[F(2, 3), F(0)], [F(-1, 7), F(1, 2)]]]
        Bs = [[[F(1)], [F(-2)]], [[F(-1)], [F(1)]]]
        metrics = [identity(2), [[F(2), F(1, 3)], [F(1, 3), F(1)]], identity(2)]
        result = storage_identity(factors, Bs, metrics[0], metrics[-1])
        self.assertEqual(result['Q'], telescoped_gap(factors, metrics))
        x, w = [[F(2)], [F(-3)]], [[F(4, 5)]]
        out = x
        for A, B in zip(factors, Bs):
            out = add(matmul(A, out), matmul(B, w))
        delta = energy(out, metrics[-1])-energy(x, metrics[0])
        linked = -energy(x, result['Q'])+2*matmul(transpose(x), matmul(result['X'], w))[0][0]+energy(w, result['Z'])
        self.assertEqual(delta, linked)
        self.assertNotEqual(out, matmul(result['M'], x))
        bound = storage_inequality(result, metrics[0], F(1, 2))
        self.assertLessEqual(energy(out, metrics[-1]), bound['rho']*energy(x, metrics[0])+energy(w, bound['supply_matrix']))

    def test_gain_residual_force_cannot_be_omitted(self):
        P = [[F(2), F(1, 3)], [F(1, 3), F(1)]]
        r = optimal_correction_storage(P, [[F(1), F(-1, 2)]], [[F(3, 5)]],
                                       [[F(1)], [F(2)]], [[F(1, 7)], [F(-1, 9)]])
        self.assertEqual(r['final'], r['initial']-r['loss']+r['cross']+r['forcing'])
        self.assertNotEqual(r['final'], r['initial']-r['loss'])

    def test_noncontracting_matrix_fails_closed(self):
        J = identity(2)
        record = storage_identity([J], [[[F(0)], [F(0)]]], J, J)
        with self.assertRaisesRegex(ValueError, 'D_SUFFICIENT_BOUND_FAILURE'):
            storage_inequality(record, J, F(1, 100))
        self.assertIsNone(certificate()['uniform_epsilon'])
        self.assertFalse(certificate()['strict_complete_word_storage_verified'])

    def test_covariance_residual_is_charged_to_same_operation_loss(self):
        P = [[F(2), F(1, 3)], [F(1, 3), F(1)]]
        H = [[F(1), F(-1, 2)]]
        dP = [[F(1, 5), F(-2, 7)], [F(-2, 7), F(-1, 3)]]
        c = covariance_residual_charge(P, H, [[F(3, 5)]], dP, [[F(2, 9)]])
        self.assertGreater(c['covariance_storage_loss'], 0)
        self.assertLessEqual(c['gain_residual_energy'], c['linked_charge_upper'])
        zero = covariance_residual_charge(P, H, [[F(3, 5)]], dP, [[F(0)]])
        self.assertEqual(zero['gain_residual_energy'], 0)

    def test_partial_joint_bound_retains_the_residual_square_charge(self):
        c = covariance_residual_charge([[F(1)]], [[F(1)]], [[F(100)]],
                                       [[F(1)]], [[F(101)]])
        # Zero mean tangent, unit covariance tangent and lambda=1. This is
        # an algebra regression, not a claimed reachable shipping operation.
        actual = c['gain_residual_energy']+1-c['covariance_storage_loss']
        correct = 1-(1-c['NIS'])*c['covariance_storage_loss']/2
        omitted_square = 1-c['covariance_storage_loss']/2
        self.assertLessEqual(actual, correct)
        self.assertGreater(actual, omitted_square)


if __name__ == '__main__':
    unittest.main()
