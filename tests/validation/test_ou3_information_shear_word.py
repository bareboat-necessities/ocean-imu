"""Rational identity regressions; none of these operands is a shipping replay."""
from fractions import Fraction as F
from pathlib import Path
import unittest

from tools.stability.ou3_theorem.information_shear_word import (
    certificate, comparison_word_budget, congruence_shear, covariance_increment_shear, joint_congruence,
    joint_metric, linked_action_balance, mean_shift_shear, precision_projector,
    projector_differential, score_loss_charge, word_score_normal_form, zeros,
)
from tools.stability.ou3_theorem.lin_path_certificate import inverse
from tools.stability.ou3_theorem.matrix_certificates import add, identity, is_psd, transpose
from tools.stability.ou3_theorem.planar_complete_word_storage import product
from tools.stability.ou3_theorem.planar_linked_riccati_mean import (
    gain_differential, optimal_covariance_differential,
)


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
        # Coordinates eta1,eta2,dP11,dP12,dP_end2; raw symmetric off diagonal.
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

    def test_complete_score_word_retains_generated_covariance_work(self):
        # Exact operation algebra, not a reached state or simulated trajectory.
        Ps, es = [self.P], [self.e]
        Bs, scores, Us, vs = [], [], [], []
        corrections = {}

        def append(B, C, ep, d, U, v):
            Bs.append(B)
            Ps.append(C)
            es.append(ep)
            scores.append(d)
            Us.append(U)
            vs.append(v)

        # Actual prediction differential, including physical mismatch and dF/dQ.
        B = [[F(1), F(1, 5)], [F(0), F(1)]]
        dB = [[F(0), F(1, 13)], [F(0), F(0)]]
        Q, dQ = [[F(1, 3), F(0)], [F(0), F(1, 7)]], [[F(1, 19), F(0)], [F(0), F(-1, 23)]]
        s, ds = [[F(1, 17)], [F(-1, 11)]], [[F(1, 29)], [F(1, 31)]]
        U = add(add(product(dB, Ps[-1], transpose(B)), product(B, Ps[-1], transpose(dB))), dQ)
        append(B, add(product(B, Ps[-1], transpose(B)), Q), add(product(B, es[-1]), s),
               zeros(2, 1), U, add(product(dB, es[-1]), ds))

        def correction(H, dH, R, dR, r, residual_mismatch):
            P, e = Ps[-1], es[-1]
            Sinv = inverse(add(product(H, P, transpose(H)), R))
            K = product(P, transpose(H), Sinv)
            A = add(identity(2), product(K, H), -1)
            C = product(A, P)
            dKcoeff = gain_differential(P, H, K, Sinv, zeros(2, 2), dH, dR)
            U = optimal_covariance_differential(P, H, K, zeros(2, 2), dH, dR)
            v = add(product(dKcoeff, r), product(K, residual_mismatch))
            corrections[len(Bs)] = (H, R, r)
            append(A, C, add(e, product(K, r)), product(transpose(H), Sinv, r), U, v)

        correction([[F(1), F(2, 3)]], [[F(1, 7), F(-1, 11)]], [[F(3, 5)]],
                   [[F(1, 17)]], [[F(2, 7)]], [[F(1, 23)]])
        # AW adds a PSD covariance increment; its face/target tangent is retained.
        append(identity(2), add(Ps[-1], [[F(0), F(0)], [F(0), F(1, 4)]]), es[-1],
               zeros(2, 1), [[F(0), F(1, 29)], [F(1, 29), F(1, 31)]], zeros(2, 1))
        correction([[F(-1, 2), F(1)]], [[F(0), F(0)]], [[F(4, 5)]],
                   [[F(0)]], [[F(-1, 3)]], [[F(0)]])
        # Covariance reset with nonzero chart discrepancy and reset derivative.
        G = [[F(1), F(1, 9)], [F(-1, 9), F(1)]]
        dG = [[F(0), F(1, 37)], [F(-1, 37), F(0)]]
        U = add(product(dG, Ps[-1], transpose(G)), product(G, Ps[-1], transpose(dG)))
        append(G, product(G, Ps[-1], transpose(G)), add(product(G, es[-1]), s),
               zeros(2, 1), U, add(product(dG, es[-1]), ds))
        append(identity(2), Ps[-1], add(es[-1], s), zeros(2, 1), zeros(2, 2), ds)
        out = word_score_normal_form(Bs, Ps, es, scores, Us, vs, self.dP, self.de)
        self.assertEqual(out['local_scores'][1], zeros(2, 1))
        self.assertEqual(out['local_scores'][3], zeros(2, 1))
        # A covariance tangent born at correction 1 contributes downstream.
        self.assertNotEqual(product(Us[1], out['suffix_scores'][2]), zeros(2, 1))
        self.assertNotEqual(out['local_scores'][2], zeros(2, 1))
        # Check the telescoped score independently of its backward recurrence.
        M, accumulated = identity(2), zeros(2, 1)
        for B, d in zip(Bs, scores):
            accumulated = add(accumulated, product(transpose(M), d))
            M = product(B, M)
        dual = add(add(product(inverse(Ps[0]), es[0]), accumulated),
                   product(transpose(M), inverse(Ps[-1]), es[-1]), -1)
        self.assertEqual(out['suffix_scores'][0], dual)
        budget = comparison_word_budget(Bs, Ps, es, corrections)
        self.assertEqual(add(budget['loss_generated_score'], budget['source_chart_score']), dual)
        self.assertLessEqual(budget['loss_score_charge'], budget['root_plus_supply_budget'])

    def test_score_charge_keeps_singular_action_kernel(self):
        P, C, M = identity(2), [[F(2), F(0)], [F(0), F(1)]], identity(2)
        out = score_loss_charge(P, C, M, self.dP, [[F(1, 3)], [F(0)]])
        self.assertEqual(out['score_charge'], F(2, 9))
        self.assertLessEqual(out['image_energy'], out['upper_bound'])
        with self.assertRaisesRegex(ValueError, 'kernel forcing'):
            score_loss_charge(P, C, M, self.dP, [[F(0)], [F(1)]])
        out = score_loss_charge(self.P, add(self.P, identity(2)), M, self.dP, self.e)
        self.assertGreater(out['covariance_loss'], 0)
        with self.assertRaisesRegex(ValueError, 'base mean action'):
            score_loss_charge(C, P, M, self.dP, self.e)

    def test_fixed_row_correction_word_cancels_arbitrary_residuals(self):
        Ps, es, Bs, scores = [self.P], [self.e], [], []
        for H, R, r in (([[F(1), F(1, 3)]], [[F(2, 5)]], [[F(7)]]),
                        ([[F(-1, 4), F(1)]], [[F(3, 7)]], [[F(-11)]])):
            P, e = Ps[-1], es[-1]
            Sinv = inverse(add(product(H, P, transpose(H)), R))
            K = product(P, transpose(H), Sinv)
            B = add(identity(2), product(K, H), -1)
            Bs.append(B)
            Ps.append(product(B, P))
            es.append(add(e, product(K, r)))
            scores.append(product(transpose(H), Sinv, r))
        out = word_score_normal_form(Bs, Ps, es, scores,
                                     [zeros(2, 2)]*2, [zeros(2, 1)]*2,
                                     self.dP, self.de)
        self.assertEqual(out['suffix_scores'], [zeros(2, 1)]*3)
        eta0 = add(self.de, product(self.dP, inverse(self.P), self.e), -1)
        self.assertEqual(out['eta_terminal'], product(Bs[1], Bs[0], eta0))

    def test_nonmeasurement_loss_score_has_linked_energy_budget(self):
        P0 = self.P
        B0 = [[F(1), F(1, 5)], [F(0), F(1)]]
        P1 = add(product(B0, P0, transpose(B0)), identity(2))
        P_end = add(P1, [[F(0), F(0)], [F(0), F(1, 3)]])
        J0, J1, J2 = inverse(P0), inverse(P1), inverse(P_end)
        Q0 = add(J0, product(transpose(B0), J1, B0), -1)
        Q1 = add(J1, J2, -1)
        e0, e1 = self.e, self.de  # Any two linked comparison operands.
        score = add(product(Q0, e0), product(transpose(B0), Q1, e1))
        out = score_loss_charge(P0, P_end, B0, self.dP, score)
        budget = product(transpose(e0), Q0, e0)[0][0]+product(transpose(e1), Q1, e1)[0][0]
        self.assertLessEqual(out['score_charge'], budget)

    def test_comparison_budget_retains_signed_work_and_physical_S(self):
        # Projection/shift work is signed, not independently squared noise.
        out = comparison_word_budget([identity(1)], [identity(1)]*2,
                                      [[[F(2)]], [[F(1)]]], {})
        self.assertEqual(out['signed_supply'], F(-3))
        self.assertEqual(out['source_chart_score'], [[F(1)]])
        # S=0 acts on nominal S=5; physical S=3 is inherited, not reset.
        H, R, r = [[F(1)]], [[F(2)]], [[F(-5)]]
        A = [[F(2, 3)]]
        out = comparison_word_budget([A], [identity(1), A],
                                      [[[F(2)]], [[F(1, 3)]]], {0: (H, R, r)})
        self.assertEqual(out['signed_supply'], F(9, 2))
        self.assertEqual(out['innovation_dissipation'], F(25, 3))
        self.assertEqual(out['root_plus_supply_budget'], F(1, 6))
        # A false reset of physical S must fail the literal additive boundary.
        with self.assertRaisesRegex(ValueError, 'literal optimal additive'):
            comparison_word_budget([A], [identity(1), A],
                                   [[[F(2)]], [[F(10, 3)]]], {0: (H, R, r)})


if __name__ == '__main__':
    unittest.main()
