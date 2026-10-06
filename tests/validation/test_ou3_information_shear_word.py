"""Rational identity regressions; none of these operands is a shipping replay."""
from fractions import Fraction as F
from pathlib import Path
import unittest

from tools.stability.ou3_theorem.information_shear_word import (
    certificate, comparison_word_budget, congruence_shear, covariance_increment_shear, joint_congruence,
    joint_metric, linked_action_balance, mean_shift_shear, precision_projector,
    projector_differential, score_loss_charge, word_score_normal_form, zeros,
    conditional_mixed_coefficients, conditional_word_mixed_bound,
    process_source_score, source_qualified_word_score, process_augmented_shear,
)
from tools.stability.ou3_theorem.lin_path_certificate import inverse
from tools.stability.ou3_theorem.matrix_certificates import add, identity, is_psd, ldlt, transpose
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

    def test_conditional_score_reader_matches_directional_closed_form(self):
        # Exact diagonal-coordinate identity, not sampled shipping covariances.
        r, q = [F(9, 10), F(1, 2), F(1, 3)], [[F(1, 7)], [F(-1, 11)], [F(1, 13)]]
        T = [[r[i]*F(i == j) for j in range(3)] for i in range(3)]
        out = conditional_mixed_coefficients(identity(3), T, q, F(2))
        reader = [[q[i][0]*q[j][0]/(2*(1-r[i]*r[j])) for j in range(3)] for i in range(3)]
        for i in range(3):
            reader[i][i] += sum(q[j][0]**2/(2*(1-r[i]*r[j])) for j in range(3))
        self.assertEqual(out['directional_score_reader'], reader)

    def test_linked_process_score_cancels_small_loss_in_isotropic_limit(self):
        # q=(1-r)e preserves the score/loss link. Exact threshold is
        # lambda > r/(1+r)|e|^2, not a product of independent extrema.
        r, q = F(9, 10), [[F(1, 30)], [F(0)], [F(0)]]
        T = [[r*x for x in row] for row in identity(3)]
        boundary = F(1, 19)
        at = conditional_mixed_coefficients(identity(3), T, q, boundary)
        below = conditional_mixed_coefficients(identity(3), T, q, boundary-F(1, 1000))
        self.assertTrue(is_psd(at['three_row_margin']))
        self.assertEqual(at['three_row_margin'][0][0], 0)
        self.assertFalse(is_psd(below['three_row_margin']))
        self.assertFalse(is_psd(below['net_conditional_Fisher_operator']))
        PN = [[F(10, 9)*x for x in row] for row in identity(3)]
        with self.assertRaisesRegex(ValueError, 'D_SUFFICIENT_BOUND_FAILURE'):
            conditional_word_mixed_bound(identity(3), PN, identity(3), q,
                zeros(3, 1), identity(3), q, identity(3), identity(3), boundary, F(0))

    def test_mixed_word_bound_keeps_generated_covariance_and_endpoint_work(self):
        # Formal rational word: prediction, optimal non-AW correction, AW
        # addition. Tests the exact complete tangent; not an admission replay.
        P = [[F(2), F(1, 5), 0, 0], [F(1, 5), F(3), F(1, 7), 0],
             [0, F(1, 7), F(2), F(1, 9)], [0, 0, F(1, 9), F(1)]]
        Fmap = identity(4)
        Fmap[0][1], Fmap[1][1], Fmap[2][2], Fmap[3][3] = F(1, 5), F(4, 5), F(4, 5), F(4, 5)
        predicted = add(product(Fmap, P, transpose(Fmap)), identity(4))
        H, R, residual = [[F(1), 0, 0, 0]], [[F(2)]], [[F(1, 7)]]
        S = add(product(H, predicted, transpose(H)), R)
        K = product(predicted, transpose(H), inverse(S))
        A = add(identity(4), product(K, H), -1)
        corrected = product(A, predicted)
        terminal = add(corrected, identity(4))
        e = [[F(1, 10)], [F(-1, 11)], [F(1, 13)], [F(1, 17)]]
        ep = add(product(Fmap, e), [[F(1, 19)], [0], [0], [0]])
        ec = add(ep, product(K, residual))
        D, de = zeros(4, 4), [[F(1, 23)], [F(1, 29)], [0], [0]]
        D[0][1] = D[1][0] = F(1, 31)
        D[1][2] = D[2][1] = F(1, 37)
        dF = zeros(4, 4)
        dF[0][1] = F(1, 41)
        U = add(product(dF, P, transpose(Fmap)), product(Fmap, P, transpose(dF)))
        Ua = zeros(4, 4)
        Ua[1][1] = F(1, 43)
        out = word_score_normal_form([Fmap, A, identity(4)], [P, predicted, corrected, terminal],
            [e, ep, ec, ec], [zeros(4, 1), product(transpose(H), inverse(S), residual), zeros(4, 1)],
            [U, zeros(4, 4), Ua], [product(dF, e), zeros(4, 1), zeros(4, 1)], D, de)
        selector = [[F(i == j+1) for j in range(3)] for i in range(4)]
        eta0 = add(de, product(D, inverse(P), e), -1)
        bound = conditional_word_mixed_bound(P, terminal, out['base_suffixes'][0], out['suffix_scores'][0],
            eta0, D, out['eta_terminal'], out['covariance_tangent_terminal'], selector, F(2), F(3, 47))
        self.assertGreaterEqual(bound['actual_signed_gap'], bound['linked_lower_bound'])
        self.assertGreaterEqual(bound['actual_signed_gap'], bound['half_loss_lower_bound'])
        self.assertEqual(bound['half_loss_lower_bound'], bound['net_conditional_root_form']/2+
                         bound['complementary_root_storage']-bound['packet_self_energy']-
                         2*bound['linked_mixed_coupling_charge']+bound['endpoint_adjustment'])
        self.assertEqual(bound['actual_signed_gap'], bound['linked_lower_bound']+
                         bound['retained_mean_square']+bound['retained_covariance_square'])
        self.assertEqual(bound['endpoint_adjustment'], F(3, 47))
        self.assertNotEqual(bound['aggregate_mean_packet'], zeros(4, 1))
        self.assertNotEqual(product(out['base_suffixes'][1], U, out['suffix_scores'][1]), zeros(4, 1))
        self.assertFalse(bound['uniform_margin_verified'])


    def test_anisotropic_linked_score_has_cross_direction_threshold(self):
        # CP1 exact rational substitution of the symbolic formula; no search,
        # independent-score maximization, or shipping reachability assertion.
        eps, amplitude = F(1, 100), F(2, 3)
        T = [[1-eps, 0, 0], [0, F(1, 2), 0], [0, 0, F(1, 2)]]
        e = [[0], [amplitude], [0]]
        q = product(add(identity(3), T, -1), e)
        out = conditional_mixed_coefficients(identity(3), T, q, F(1))
        critical = product(out['completed_score_weight'], out['directional_score_reader'])
        expected = [amplitude**2*(1-eps)/(4*eps*(1+eps)),
                    amplitude**2/3, amplitude**2/6]
        self.assertEqual(critical, [[expected[i]*F(i == j) for j in range(3)] for i in range(3)])
        self.assertEqual(product(transpose(q), inverse(add(identity(3), T, -1)), q)[0][0], amplitude**2/2)

    def test_base_schur_absorbs_complement_where_fixed_half_charge_fails(self):
        # CP3--CP6: exact formal base map, NOT an admitted shipping execution.
        M = [[F(3, 4), F(1, 5)], [F(1, 5), F(3, 4)]]
        ldlt(add(identity(2), product(transpose(M), M), -1))
        eta = [[0], [F(1)]]
        out = conditional_word_mixed_bound(identity(2), identity(2), M,
            zeros(2, 1), eta, zeros(2, 2), product(M, eta), zeros(2, 2),
            [[F(1)], [0]], F(1), F(0))
        self.assertEqual(out['actual_signed_gap'], F(159, 400))
        self.assertEqual(out['packet_self_energy'], F(241, 400))
        self.assertEqual(out['linked_mixed_coupling_charge'], F(12, 53))
        self.assertEqual(out['half_loss_lower_bound'], -F(1173, 21200))
        self.assertEqual(out['linked_lower_bound'], F(3627, 21200))
        # Check the full joint mean/Fisher base Schur identity, including
        # off-diagonal covariance variation, not only the displayed mean slice.
        H = joint_metric(identity(2), F(1))
        transport = joint_congruence(M)
        gap = add(H, product(transpose(transport), H, transport), -1)
        c, o = [0, 2], [1, 3, 4]  # eta_0,D_00 and their metric complement
        take = lambda rows, cols: [[gap[i][j] for j in cols] for i in rows]
        A, cross, complement = take(c, c), take(c, o), take(o, o)
        schur = add(complement, product(transpose(cross), inverse(A), cross), -1)
        ldlt(A)
        ldlt(schur)
        self.assertEqual(A[0][0], F(159, 400))
        self.assertEqual(schur[0][0], F(3627, 21200))

    def test_combined_process_score_preserves_supported_physical_defect(self):
        P = [[F(2), F(1, 3)], [F(1, 3), F(1)]]
        B = [[F(1), F(1, 5)], [0, F(4, 5)]]
        Q = [[F(1, 2), 0], [0, 0]]
        e, s = [[F(1, 7)], [F(-1, 11)]], [[F(2, 9)], [0]]
        out = process_source_score(P, B, Q, e, s)
        self.assertEqual(out['physical_process_action'], F(8, 81))
        self.assertEqual(product(out['mean_loss'], out['range_witness']), out['combined_score'])
        loss_only = product(out['mean_loss'], e)
        self.assertNotEqual(out['combined_score'], loss_only)
        with self.assertRaisesRegex(ValueError, 'outside process-noise range'):
            process_source_score(P, B, Q, e, [[0], [F(1, 19)]])

    def test_combined_word_score_keeps_uncovered_chart_and_S_defect(self):
        # Exact formal prediction/S-correction/chart/AW word. Physical S is
        # nonzero: delta=r+H e is retained, never reset by the pseudo-update.
        P, B, Q = [[F(2)]], [[F(4, 5)]], [[F(1, 3)]]
        e, s = [[F(1, 4)]], [[F(-1, 7)]]
        proc = process_source_score(P, B, Q, e, s)
        C, ep = proc['covariance_next'], proc['comparison_next']
        H, R, r = [[F(1)]], [[F(2)]], [[F(-1, 5)]]
        K = product(C, transpose(H), inverse(add(product(H, C, transpose(H)), R)))
        A = add(identity(1), product(K, H), -1)
        PC = product(A, C)
        ec = add(ep, product(K, r))
        shifted = add(ec, [[F(1, 13)]])
        PN = add(PC, [[F(1, 17)]])
        maps, covs = [B, A, identity(1), identity(1)], [P, C, PC, PC, PN]
        errors = [e, ep, ec, shifted, shifted]
        out = source_qualified_word_score(maps, covs, errors, {1: (H, R, r)}, [0, 3])
        self.assertNotEqual(out['uncovered_score'], zeros(1, 1))
        self.assertEqual(out['physical_process_action'], F(3, 49))
        self.assertEqual(out['correction_defect_action'], product(transpose(add(r, ep)), inverse(R), add(r, ep))[0][0])
        self.assertGreaterEqual(out['linked_budget'], out['covered_score_charge'])
        direct = word_score_normal_form(maps, covs, errors,
            [zeros(1, 1), product(transpose(H), inverse(add(C, R)), r), zeros(1, 1), zeros(1, 1)],
            [zeros(1, 1)]*4, [zeros(1, 1)]*4, zeros(1, 1), zeros(1, 1))
        self.assertEqual(out['actual_total_score'], direct['suffix_scores'][0])

    def test_augmented_prediction_retains_all_causal_derivative_ports(self):
        # Exact noncommuting derivative identity, not a sampled coefficient box.
        P = [[F(2), F(1, 5)], [F(1, 5), F(1)]]
        B = [[F(1), F(1, 4)], [0, F(4, 5)]]
        Q = [[F(1, 2), F(1, 9)], [F(1, 9), F(1, 3)]]
        e, s = [[F(2, 7)], [F(-1, 11)]], [[F(1, 13)], [F(1, 17)]]
        D = [[F(1, 19), F(1, 23)], [F(1, 23), F(-1, 29)]]
        de, ds = [[F(1, 31)], [F(1, 37)]], [[F(-1, 41)], [F(1, 43)]]
        dB = [[0, F(1, 47)], [0, F(1, 53)]]
        dQ = [[F(1, 59), F(-1, 61)], [F(-1, 61), F(1, 67)]]
        out = process_augmented_shear(P, B, Q, e, s, D, de, dB, dQ, ds, F(2))
        self.assertLessEqual(out['feedback_square'], out['feedback_square_upper'])
        self.assertEqual(out['actual_storage_change'], out['signed_auxiliary_work']-
                         out['hidden_mean_loss']-2*out['augmented_Fisher_loss']+
                         out['signed_feedback_cross']+out['feedback_square'])
        self.assertNotEqual(out['signed_feedback_cross'], 0)
        no_aux = process_augmented_shear(P, B, Q, e, s, D, de,
                                         zeros(2, 2), zeros(2, 2), zeros(2, 1), F(2))
        self.assertEqual(no_aux['signed_auxiliary_work'], 0)
        self.assertNotEqual(out['shear_next'], no_aux['shear_next'])
        self.assertFalse(out['uniform_absorption_verified'])



    def test_supported_process_score_is_invariant_in_endpoint_frames(self):
        P = [[F(2), F(1, 3)], [F(1, 3), F(1)]]
        B = [[F(1), F(1, 5)], [0, F(4, 5)]]
        Q = [[F(1, 2), 0], [0, 0]]
        e, s = [[F(1, 7)], [F(-1, 11)]], [[F(2, 9)], [0]]
        root, end = [[1, 0], [F(1, 13), 1]], [[1, 0], [F(-1, 17), 1]]
        out = process_source_score(P, B, Q, e, s)
        moved = process_source_score(product(root, P, transpose(root)),
            product(end, B, inverse(root)), product(end, Q, transpose(end)),
            product(root, e), product(end, s))
        self.assertEqual(moved['combined_score'], product(transpose(inverse(root)), out['combined_score']))
        self.assertEqual(moved['combined_score_charge'], out['combined_score_charge'])
        self.assertEqual(moved['physical_process_action'], out['physical_process_action'])


if __name__ == '__main__':
    unittest.main()
