from fractions import Fraction as F
import unittest
from tools.stability.ou3_theorem.matrix_certificates import add, transpose
from tools.stability.ou3_theorem.planar_linked_riccati_mean import product
from tools.stability.ou3_theorem.planar_causal_calculus import (
    certificate, covariance_prediction_differential, word_differential,
    information_differential, information_perturbation_bound, probe_word_differential,
)
from tests.validation.test_ou3_planar_linked_riccati_mean import jet, deriv, inverse2


class CausalCalculusTests(unittest.TestCase):
    def test_exact_word_and_prediction_differentials(self):
        a = [[F(1), F(1, 3)], [F(-1, 5), F(2)]]
        da = [[F(1, 10), F(0)], [F(1, 7), F(-1, 20)]]
        p = [[F(2), F(1, 4)], [F(1, 4), F(3)]]
        dp = [[F(1, 11), F(1, 13)], [F(1, 13), F(1, 17)]]
        _, dw = word_differential([a, p, transpose(a)], [da, dp, transpose(da)])
        self.assertEqual(dw, deriv(product(transpose(jet(a, da)), jet(p, dp), jet(a, da))))
        root, droot = [[F(1)], [F(2)]], [[F(1, 3)], [F(-1, 5)]]
        _, dprobe = probe_word_differential([a, p], [da, dp], root, droot)
        self.assertEqual(dprobe, deriv(product(jet(p, dp), jet(a, da), jet(root, droot))))
        dq = [[F(1, 100), F(0)], [F(0), F(1, 200)]]
        self.assertEqual(covariance_prediction_differential(a, p, da, dp, dq),
                         add(deriv(product(jet(a, da), jet(p, dp), transpose(jet(a, da)))), dq))

    def test_exact_information_differential(self):
        y = [[F(1), F(1, 3)], [F(1, 5), F(2)]]
        dy = [[F(1, 10), F(0)], [F(1, 7), F(-1, 20)]]
        s = [[F(2), F(1, 4)], [F(1, 4), F(3)]]
        ds = [[F(1, 11), F(1, 13)], [F(1, 13), F(1, 17)]]
        self.assertEqual(information_differential(y, inverse2(s), dy, ds),
                         deriv(product(transpose(jet(y, dy)), inverse2(jet(s, ds)), jet(y, dy))))

    def test_conditional_bound_does_not_invent_a_floor(self):
        self.assertEqual(information_perturbation_bound(center_stack_norm=F(2),
                         stack_difference_norm=F(1, 10), relative_innovation_radius=F(1, 5)), F(121, 80))
        with self.assertRaises(ValueError):
            information_perturbation_bound(center_stack_norm=1, stack_difference_norm=0, relative_innovation_radius=1)
        c = certificate()
        self.assertIsNone(c['source_uniform_storage_contraction'])
        self.assertIsNone(c['certified_radius'])
        self.assertIsNone(c['center_service_floor'])
        self.assertFalse(c['all_time_magnetic_service_verified'])


if __name__ == '__main__':
    unittest.main()
