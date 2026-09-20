from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.stability.ou3_theorem.corrected_word import (
    certificate, conditional_scalar_margin, coupled_example, nuisance_root_bounds,
    polynomial_injection_defect, reset_remainder_bound, shipping_reset_remainder_bound,
)
from tools.stability.ou3_theorem.lin_path_certificate import inverse
from tools.stability.ou3_theorem.matrix_certificates import (
    add, congruence, identity, is_psd, matmul, transpose,
)
from tools.stability.ou3_theorem.word_energy import word_identity


class CoupledCovarianceTests(unittest.TestCase):
    def test_conditional_upper_bound_keeps_full_21_state_cross_covariance(self):
        c = [[F((3*i+2*j) % 5-2, 5) for j in range(15)] for i in range(6)]
        a = add(identity(6), matmul(c, transpose(c)))
        n = [[2*x for x in row] for row in identity(15)]
        p = [a[i]+c[i] for i in range(6)] + [transpose(c)[i]+n[i] for i in range(15)]
        conditional = add(a, matmul(c, transpose(c)), F(-1, 2))
        # alpha=1/2, U=2I, eta=1: retain the full conditional matrix.
        ag_upper = [[4*x for x in row] for row in conditional]
        upper = [ag_upper[i]+[F(0)]*15 for i in range(6)]
        upper += [[F(0)]*6+[4*x for x in row] for row in identity(15)]
        self.assertTrue(is_psd(add(upper, p, F(-1))))
        self.assertTrue(any(x for row in c for x in row))

    def test_marginal_floor_is_not_an_embedded_floor(self):
        # Conditional AG precision is exactly one and nuisance variance is one,
        # but AG variance is arbitrarily large. A marginal floor cannot replace
        # the full Loewner embedded comparison used in the proof.
        p = [[F(10001), F(100)], [F(100), F(1)]]
        self.assertEqual(inverse(p)[0][0], 1)
        self.assertFalse(is_psd(add(p, [[0, 0], [0, 1]], F(-1))))

    def test_source_pre_prediction_nuisance_floor(self):
        lower, ba, upper, alpha = nuisance_root_bounds()
        self.assertTrue(is_psd(add(lower, [[upper[i] if i == j else 0
                                         for j in range(4)] for i in range(4)], -alpha)))
        self.assertGreaterEqual(ba, alpha*upper[4])
        self.assertGreater(alpha, 0)
        self.assertLess(alpha, 1)

    def test_exact_prediction_closes_cross_coupled_example(self):
        example = coupled_example()
        self.assertTrue(example['exact_full_loss_check'])
        self.assertFalse(example['shipping_history'])
        self.assertGreater(F(example['strict_full_loss_margin']), 0)

    def test_zero_process_or_missing_six_column_margin_cannot_promote(self):
        for mu, q in ((0, F('.01')), (F('.008'), 0)):
            with self.assertRaises(ValueError):
                conditional_scalar_margin(mu, F(1, 3), 3, q, 1)
        report = certificate()
        for key in ('six_column_source_uniform_loss_verified',
                    'full_21_covariance_upper_verified',
                    'source_uniform_A21_linear_dissipativity', 'theorem_closed'):
            self.assertFalse(report[key])

    def test_actual_correction_supply_without_gain_norm_bound(self):
        p = [[F(100), F(9)], [F(9), F(1)]]
        h, r = [[F(1), F(2)]], [[F(3)]]
        s = add(congruence(p, transpose(h)), r)
        k = matmul(matmul(p, transpose(h)), inverse(s))
        a = add(identity(2), matmul(k, h), F(-1))
        post = add(congruence(p, transpose(a)), congruence(r, transpose(k)))
        self.assertTrue(is_psd(add(inverse(r), congruence(inverse(post), k), F(-1))))

    def test_joint_input_gramian_with_correlated_noise_and_nonorthogonal_reset(self):
        p = [[F(2), F(1)], [F(1), F(3)]]
        f = [[F(1), F(1, 5)], [F(0), F(1)]]
        q = [[F(1), F(1, 3)], [F(1, 3), F(1)]]
        h, r = [[F(1), F(2)]], [[F(3, 2)]]
        g = [[F(1), F(1, 10)], [F(-1, 10), F(1)]]
        word = word_identity(p, [
            {'kind': 'prediction', 'F': f, 'Q': q},
            {'kind': 'correction', 'H': h, 'R': r},
            {'kind': 'reset', 'G': g}])
        pm = add(congruence(p, transpose(f)), q)
        s = add(congruence(pm, transpose(h)), r)
        k = matmul(matmul(pm, transpose(h)), inverse(s))
        a = add(identity(2), matmul(k, h), F(-1))
        accumulated = congruence(add(congruence(q, transpose(a)),
                                    congruence(r, transpose(k))), transpose(g))
        propagated_root = congruence(p, transpose(word['end_transport']))
        self.assertEqual(add(propagated_root, accumulated), word['end_covariance'])
        self.assertTrue(is_psd(add(word['end_covariance'], accumulated, F(-1))))


class ResetRemainderTests(unittest.TestCase):
    def test_fixed_injection_remainder_is_not_quadratic_in_error(self):
        d = np.array([.3, 0, 0])
        direction = np.array([0, 1, 0])
        exact_derivative = np.array([0, np.sin(.3)/.3, (1-np.cos(.3))/.3])
        reset_derivative = direction+.5*np.cross(d, direction)
        limit = np.linalg.norm(exact_derivative-reset_derivative)
        self.assertGreater(limit, .014)
        for step in (1e-4, 1e-6):
            v = step*direction
            actual = (Rotation.from_rotvec(d+v)*Rotation.from_rotvec(-d)).as_rotvec()
            ratio = np.linalg.norm(actual-step*reset_derivative)/step
            self.assertAlmostEqual(ratio, limit, delta=1e-5)

    def test_noncommuting_and_large_realized_injections(self):
        # Independent SO(3) implementation; a large d is allowed even though
        # its resulting remainder ceiling need not be useful for retention.
        for d in ([.001, -.002, .003], [.2, .1, -.3], [2.3, -1.7, .8], [6.4, .1, -.2]):
            d = np.array(d)
            for v in ([.01, -.02, .03], [.3, -.1, .2], [1.7, .2, -.1]):
                v = np.array(v)
                actual = (Rotation.from_rotvec(d+v)*Rotation.from_rotvec(-d)).as_rotvec()
                comparison = v+.5*np.cross(d, v)
                bound = float(reset_remainder_bound(str(np.linalg.norm(d)),
                              str(np.linalg.norm(v)), '1.9'))
                self.assertLessEqual(np.linalg.norm(actual-comparison), bound+1e-14)

    def test_source_polynomial_defect_and_log_chart_guards(self):
        self.assertEqual(polynomial_injection_defect(0), 0)
        self.assertGreater(polynomial_injection_defect(F('.009')), 0)
        self.assertEqual(polynomial_injection_defect(F('.01')), 0)
        base = reset_remainder_bound(F('.009'), F('.1'), F('.2'))
        shipping = shipping_reset_remainder_bound(F('.009'), F('.1'), F('.2'))
        self.assertGreater(shipping, base)
        with self.assertRaises(ValueError):
            reset_remainder_bound(0, 1, F('.5'))
        with self.assertRaises(ValueError):
            shipping_reset_remainder_bound(0, 1, 2)


if __name__ == '__main__':
    unittest.main()
