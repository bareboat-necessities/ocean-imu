"""Rational identity regressions; none of these operands is a shipping replay."""
from fractions import Fraction as F
from pathlib import Path
import unittest

from tools.stability.ou3_theorem.information_shear_word import (
    certificate, congruence_shear, covariance_increment_shear, joint_congruence,
    joint_metric, linked_action_balance, mean_shift_shear, precision_projector,
    projector_differential, zeros,
)
from tools.stability.ou3_theorem.lin_path_certificate import inverse
from tools.stability.ou3_theorem.matrix_certificates import add, identity, is_psd, transpose
from tools.stability.ou3_theorem.planar_complete_word_storage import product


class InformationShearWordTests(unittest.TestCase):
    def setUp(self):
        self.P = [[F(2), F(1, 3)], [F(1, 3), F(1)]]
        self.dP = [[F(1, 5), F(-2, 7)], [F(-2, 7), F(1, 4)]]
        self.e, self.de = [[F(2, 3)], [F(-1, 5)]], [[F(1, 7)], [F(3, 8)]]
        self.w = F(3, 5)

    def test_optimal_action_kernel_includes_covariance_rows(self):
        P, H, R = identity(2), [[F(1), F(0)]], [[F(1)]]
        S = add(product(H, P, transpose(H)), R)
        A = add(identity(2), product(P, transpose(H), inverse(S), H), -1)
        C = product(A, P)
        L = joint_congruence(A)
        result = linked_action_balance([L], [joint_metric(P, self.w), joint_metric(C, self.w)], [identity(5), L])
        Q = result['action']
        # Coordinates eta1,eta2,dP11,dP12,dP22; raw symmetric off diagonal.
        self.assertEqual([Q[i][i] for i in range(5)], [F(1, 2), 0, 3*self.w/4, self.w, 0])
        self.assertEqual(result['signed_port_work'], zeros(5, 5))
        for v in ([[0], [1], [0], [0], [1]], [[0], [0], [0], [1], [0]]):
            loss = product(transpose(v), Q, v)[0][0]
            self.assertEqual(loss == 0, v[3][0] == 0)

    def test_rectangular_causal_root_preserves_signed_cross_terms(self):
        # Root has one auxiliary coordinate. No false mean/P Markov closure.
        J0 = joint_metric(self.P, self.w)
        Qprocess = [[F(1, 4), F(0)], [F(0), F(1, 5)]]
        C = add(self.P, Qprocess)
        J1 = joint_metric(C, self.w)
        G = [[F(1), F(1, 3)], [F(0), F(1)]]
        J2 = joint_metric(product(G, C, transpose(G)), self.w)
        L0, L1 = identity(5), joint_congruence(G)
        T0 = [row+[F(0)] for row in identity(5)]
        E0, E1 = zeros(5, 6), zeros(5, 6)
        E0[0][2], E0[1][5] = F(1, 7), F(-2, 5)
        E1[0][5], E1[2][0] = F(2, 9), F(-1, 8)
        T1 = add(T0, E0)
        T2 = add(product(L1, T1), E1)
        result = linked_action_balance([L0, L1], [J0, J1, J2], [T0, T1, T2])
        self.assertEqual(result['root_port_maps'], [E0, E1])
        self.assertNotEqual(result['signed_port_work'], zeros(6, 6))
        self.assertEqual(result['gap'], add(result['initial'], result['terminal'], -1))
        self.assertEqual(result['local_losses'][1], zeros(5, 5))

    def test_positive_base_action_does_not_prove_actual_gap(self):
        # Exact logical regression: an uncancelled endogenous port can reverse
        # a strictly positive base loss. This is NOT an admitted counterexample.
        result = linked_action_balance([[[F(1, 2)]]], [[[F(1)]], [[F(1)]]], [[[F(1)]], [[F(2)]]])
        self.assertEqual(result['action'], [[F(3, 4)]])
        self.assertEqual(result['gap'], [[F(-3)]])
        with self.assertRaisesRegex(ValueError, 'D_SUFFICIENT_BOUND_FAILURE'):
            linked_action_balance([[[F(2)]]], [[[F(1)]], [[F(1)]]], [[[F(1)]], [[F(2)]]])

    def test_reset_shear_retains_actual_chart_discrepancy(self):
        G = [[F(1), F(1, 3)], [F(-1, 3), F(1)]]
        dG = [[F(0), F(1, 7)], [F(-1, 7), F(0)]]
        t, dt = [[F(1, 13)], [F(-1, 17)]], [[F(2, 19)], [F(1, 23)]]
        out = congruence_shear(self.P, G, self.e, t, self.dP, dG, self.de, dt)
        self.assertNotEqual(out['eta_plus'], out['base'])
        L = joint_congruence(G)
        self.assertEqual(product(transpose(L), joint_metric(out['covariance_plus'], self.w), L), joint_metric(self.P, self.w))

    def test_aw_derivative_and_projection_have_linked_ports(self):
        U = [[F(0), F(0)], [F(0), F(2, 5)]]
        dU = [[F(0), F(1, 7)], [F(1, 7), F(-1, 5)]]
        out = covariance_increment_shear(self.P, U, self.e, self.dP, dU, self.de)
        without = covariance_increment_shear(self.P, U, self.e, self.dP, zeros(2, 2), self.de)
        self.assertNotEqual(out['eta_plus'], without['eta_plus'])
        self.assertTrue(is_psd(add(joint_metric(self.P, self.w), joint_metric(add(self.P, U), self.w), -1)))
        shift, dshift = [[F(0)], [F(-1, 3)]], [[F(0)], [F(1, 11)]]
        proj = mean_shift_shear(self.P, self.e, shift, self.dP, self.de, dshift)
        self.assertNotEqual(proj['port'], dshift)

    def test_moving_projector_derivative_and_endpoint_gauge_balance(self):
        J = inverse(self.P)
        dJ = [[-x for x in row] for row in product(J, self.dP, J)]
        R, dR = [[F(1)], [F(2)]], [[F(1, 7)], [F(-1, 3)]]
        Pi = precision_projector(J, R)
        dPi = projector_differential(J, R, dJ, dR)
        self.assertEqual(add(product(dPi, Pi), product(Pi, dPi)), dPi)
        self.assertEqual(add(product(dPi, R), product(Pi, dR)), dR)
        self.assertEqual(add(product(transpose(dPi), J), product(transpose(Pi), dJ)),
                         add(product(dJ, Pi), product(J, dPi)))
        N = add(identity(2), Pi, -1)
        self.assertEqual(add(product(transpose(N), J, N), product(transpose(Pi), J, Pi)), J)
        with self.assertRaises((ValueError, ArithmeticError, ZeroDivisionError)):
            precision_projector(J, [[F(0)], [F(0)]])

    def test_source_scope_and_fail_closed_status(self):
        repo = Path(__file__).resolve().parents[2]
        source = (repo/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h').read_text()
        injection = source.split('::applyQuaternionCorrectionFromErrorState()')[1].split('/*')[0]
        self.assertLess(injection.index('project_gyro_bias_();'), injection.index('qref = corr * qref;'))
        self.assertLess(injection.index('apply_error_state_reset_jacobian_(dtheta);'), injection.index('project_acc_bias_();'))
        self.assertIn('Pi_+(Sigma_aw_stat - P_awaw^-)', source)
        c = certificate()
        self.assertIsNone(c['uniform_c'])
        self.assertFalse(c['general_theorem_waits_for_planar_admission'])
        self.assertFalse(c['physical_compatibility_equals_homogeneous_kernel'])
        self.assertFalse(c['radius_solve_performed'])
        self.assertFalse(c['theorem_closed'])


if __name__ == '__main__':
    unittest.main()
