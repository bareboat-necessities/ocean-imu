"""Exact rational identities, not trajectory or contraction experiments."""
from fractions import Fraction as F
from pathlib import Path
import unittest

from tools.stability.ou3_theorem.information_shear_word import zeros, word_score_normal_form
from tools.stability.ou3_theorem.information_shear_word import conditional_mixed_coefficients
from tools.stability.ou3_theorem.lin_path_certificate import inverse
from tools.stability.ou3_theorem.matrix_certificates import add, identity, transpose
from tools.stability.ou3_theorem.measurement_frame import (
    certificate, connection, correction, covariance_row_charge, frame,
    magnetic_connection_charge, planar_magnetic_loss_margin, planar_pitch_domain_bounds,
    pullback_differentials, rank_one_magnetic_balance,
    row_differentials, world_rows, aw_connection, aw_shear,
    aw_shear_differentials, row_connection_coboundary, fixed_row_joint_balance,
    planar_acc_mismatch_charge,
    moving_frame_word_ports, aw_frame_boundary_work,
    aw_conditional_sync,
    conditional_aw_storage, aw_precision_process_balance, qualified_aw_precision_ceiling,
    conditional_aw_process_coercivity,
    planar_pitch_prediction_calculus,
    scalar_aw_face_fisher_balance, covariance_word_signed_matrix,
    odd_covariance_gap_scope_check,
)
from tools.stability.ou3_theorem.planar_innovation_storage import information_shear_correction
from tools.stability.ou3_theorem.planar_linked_riccati_mean import product
from tools.stability.ou3_theorem.planar_linked_riccati_mean import gain_differential
from tools.stability.ou3_theorem.world_frame import quaternion_rotation, skew


class MeasurementFrameTests(unittest.TestCase):
    def test_actual_scalar_face_retains_signed_regression_work(self):
        # Identity operands, not a reached shipping covariance or target.
        P = [[F(1), F(1)], [F(1), F(2)]]
        D = [[F(1), F(1)], [F(1), F(1)]]
        out = scalar_aw_face_fisher_balance(P, D, F(4), 1)
        self.assertEqual(out['dC'], 0)
        self.assertEqual(out['dT'], [[F(0)]])
        self.assertEqual(out['positive_face_action'], 0)
        self.assertEqual(out['adverse_face_work'], F(1, 9))
        self.assertEqual(out['gap'], -F(1, 9))
        self.assertEqual(out['D_next'], [[F(1), F(1)], [F(1), F(0)]])
        with self.assertRaisesRegex(ValueError, 'strictly active'):
            scalar_aw_face_fisher_balance(P, D, F(2), 1)
        # A cross-covariance replacement is not this shipping increment.
        with self.assertRaisesRegex(ValueError, 'retain cross covariance'):
            covariance_word_signed_matrix(
                [P, [[F(1), F(0)], [F(0), F(4)]]],
                [[D], [out['D_next']]], [0], 1)

    def test_complete_signed_word_fails_automatic_face_absorption_on_fixed_AA_fibre(self):
        out = odd_covariance_gap_scope_check()
        self.assertEqual(out['signed_gap'], '-691428110973/567390082009')
        self.assertGreater(F(out['relative_reader_exact']), 1)
        self.assertTrue(out['root_AW_marginal_tangent_zero'])
        self.assertTrue(out['strict_positive_process'])
        self.assertEqual(out['classification'], 'D_SUFFICIENT_BOUND_FAILURE')
        self.assertFalse(out['shipping_counterexample'])
        self.assertIsNone(certificate()['odd_covariance_complete_gap']['uniform_odd_covariance_gap'])

    def test_odd_covariance_gain_port_cannot_change_planar_mean(self):
        # Both covariance parities retained. Literal acc/mag/S row pattern
        # and r_y=0 give the source proof's exact dK*r cancellation.
        odd = (0, 2, 3, 5, 7, 10, 13, 16, 19)
        D = zeros(21, 21)
        for j, i in enumerate(odd):
            D[i][i] = F(j+1, 17)
            D[i][odd[(j+1) % 9]] = D[odd[(j+1) % 9]][i] = F(1, 19)
        acc, mag = world_rows([F(1, 5), 0, -10], [7, 0, F(1, 9)])
        Hs = zeros(3, 21)
        for j in range(3):
            Hs[j][12+j] = F(1)
        P, R, r = self.P, self.noise, [[F(2, 5)], [F(0)], [F(3, 8)]]
        for H in (acc, mag, Hs):
            S = add(product(H, P, transpose(H)), R)
            K = product(P, transpose(H), inverse(S))
            dK = gain_differential(P, H, K, inverse(S), D, zeros(3, 21), zeros(3, 3))
            self.assertEqual(product(dK, r), zeros(21, 1))
            A = add(identity(21), product(K, H), -1)
            Dnext = product(A, D, transpose(A))
            self.assertTrue(all(Dnext[i][j] == 0 for i in range(21) if i not in odd
                                for j in range(21)))

    def test_same_operation_directional_mag_S_Fisher_loss(self):
        P = [[F(2), F(1, 3)], [F(1, 3), F(1)]]
        D = [[F(1, 5), F(-2, 7)], [F(-2, 7), F(1, 4)]]
        h, R = [[F(1), F(2, 3)]], [[F(4, 5)]]
        S = add(product(h, P, transpose(h)), R)
        A = add(identity(2), product(P, transpose(h), inverse(S), h), -1)
        Pn, Dn = product(A, P), product(A, D, transpose(A))
        out = covariance_word_signed_matrix([P, Pn], [[D], [Dn]], [], 1)
        J, s = inverse(P), S[0][0]
        expected = (2*product(h, D, J, D, transpose(h))[0][0]/s
                    - product(h, D, transpose(h))[0][0]**2/s**2)
        self.assertEqual(out['positive_action'][0][0], expected)
        self.assertEqual(out['signed_gap'], out['positive_action'])
        self.assertGreater(expected, 0)

    def test_literal_planar_pitch_process_and_quaternion_defect(self):
        h, x, qg, qb = F(3, 500), F(7, 1000), F(135, 100000)**2, F(1, 10**10)
        out = planar_pitch_prediction_calculus(h, x, qg, qb)
        self.assertEqual(out['pitch_BG_transition'], [[1, h], [0, 1]])
        self.assertEqual(out['quaternion_norm_squared'],
                         1+x**6/F(23040)-x**8/F(245760)+x**10/F(14745600))
        self.assertLess(out['literal_angle_derivative'], 1)
        self.assertLess(out['auxiliary_charge_upper'], F(1, 10**27))
        # Symbolically derived cap evaluated with exact rationals, not a grid.
        self.assertLess(F(3, 500)*F(1, 7680*10**12)**2/F(1, 10**6), F(1, 10**27))
        at_zero = planar_pitch_prediction_calculus(h, 0, qg, qb)
        self.assertEqual(at_zero['literal_angle_derivative'], 1)
        with self.assertRaisesRegex(ValueError, 'small-angle'):
            planar_pitch_prediction_calculus(h, F(1, 100), qg, qb)

    def test_planar_axis_annihilates_all_literal_integral_coefficients(self):
        # Arbitrary rational coefficients test the coefficient-free identity
        # W e_y=W^2 e_y=0, including either literal coefficient branch.
        h, w = F(3, 500), F(7, 13)
        W, ey = skew([0, w, 0]), [[0], [1], [0]]
        W2 = product(W, W)
        g0, g1, g2, g3 = F(2, 7), F(3, 11), F(-1, 13), F(5, 17)
        scale = lambda A, a: [[a*x for x in row] for row in A]
        R = add(add(identity(3), scale(W, -h*g0)), scale(W2, h*h*g1))
        B = add(add(scale(identity(3), h), scale(W, -h*h*g1)), scale(W2, h**3*g2))
        IB = add(add(scale(identity(3), h*h/2), scale(W, -h**3*g2)), scale(W2, h**4*g3))
        self.assertEqual(product(R, ey), ey)
        self.assertEqual(product(B, ey), scale(ey, h))
        self.assertEqual(product(IB, ey), scale(ey, h*h/2))
        self.assertNotEqual(R[0][2], 0)  # odd block is NOT constant in w

    def test_planar_process_calculus_is_bound_to_literal_source(self):
        root = Path(__file__).resolve().parents[2]
        core = (root/'src/kalman_ou_common/KalmanOUCoreMath.h').read_text()
        shipping = (root/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h').read_text()
        for text in ('w = std::fma(-t2, T(1)/T(8), w);',
                     'w = std::fma( t4, T(1)/T(384), w);',
                     'k = std::fma(-t2, T(1)/T(48), k);',
                     'k = std::fma( t4, T(1)/T(3840), k);', 'q.normalize();'):
            self.assertIn(text, core)
        for text in ('F_AA.template block<3,3>(0,3) = Bstep;',
                     'Qbg = Qbase.template bottomRightCorner<3,3>();',
                     'I_BB = simpson_B_Q_BT_(w, Ts, Qbg);'):
            self.assertIn(text, shipping)

    def setUp(self):
        self.R = quaternion_rotation([5, 0, 1, 0])
        self.T = frame(self.R)
        self.P = identity(21)
        self.P[1][18] = self.P[18][1] = F(1, 5)
        self.P[4][15] = self.P[15][4] = F(-1, 7)
        self.noise = [[F(i == j, 3) for j in range(3)] for i in range(3)]
        self.r = [[F(2, 5)], [F(-1, 7)], [F(3, 8)]]
        self.acc, self.mag = world_rows([F(1, 5), 0, -10], [7, 0, F(1, 9)])

    def test_full_21_state_congruence_including_literal_held_mask(self):
        body_P = product(self.T, self.P, transpose(self.T))
        body_H = product(self.R, self.acc, transpose(self.T))
        body_r = product(self.R, self.r)
        for held in (False, True):
            body = correction(body_P, body_H, self.noise, body_r, held)
            world = correction(self.P, self.acc, self.noise, self.r, held)
            self.assertEqual(product(transpose(self.T), body['C'], self.T), world['C'])
            self.assertEqual(product(transpose(self.T), body['K'], self.R), world['K'])
            self.assertEqual(product(transpose(self.T), body['increment']), world['increment'])
        # Held cross blocks are deliberately nonzero here: congruence applies,
        # but no optimal-Riccati/Fisher identity is asserted for that mask.

    def test_differentiated_frame_cancels_attitude_rows_not_aw_or_noise(self):
        omega = [F(1, 7), F(-1, 5), F(2, 9)]
        w, o = connection(omega)
        dHw, _ = row_differentials([F(1, 3), F(1, 8), 0], [0, 0, 0])
        body_H = product(self.R, self.acc, transpose(self.T))
        body_P = product(self.T, self.P, transpose(self.T))
        body_dH = product(self.R, add(add(product(w, self.acc), dHw),
                                     product(self.acc, o), -1), transpose(self.T))
        dPw = zeros(21, 21)
        dPw[1][15] = dPw[15][1] = F(1, 11)
        body_dP = product(self.T, add(add(dPw, product(o, self.P)),
                                     product(self.P, o), -1), transpose(self.T))
        pulled = pullback_differentials(self.R, omega, body_P, body_H, self.noise,
            product(self.R, self.r), body_dP, body_dH, zeros(3, 3), zeros(3, 1))
        self.assertEqual(pulled['dP'], dPw)
        self.assertEqual(pulled['dH'], dHw)
        self.assertEqual(pulled['dR'], zeros(3, 3))
        anisotropic = [[F(i+1)*F(i == j) for j in range(3)] for i in range(3)]
        p = pullback_differentials(self.R, omega, body_P, body_H, anisotropic,
            product(self.R, self.r), body_dP, body_dH, zeros(3, 3), zeros(3, 1))
        self.assertNotEqual(p['dR'], zeros(3, 3))

    def test_shear_connection_and_linked_magnetic_port(self):
        # Arbitrary comparison is retained; this checks the moving-frame
        # derivative rather than pretending it is a fixed orthogonal change.
        e, de, dP = zeros(21, 1), zeros(21, 1), zeros(21, 21)
        e[1][0], e[18][0], de[1][0] = F(1, 9), F(2, 11), F(1, 7)
        dP[1][18] = dP[18][1] = F(1, 13)
        omega = [0, de[1][0], 0]
        w, o = connection(omega)
        dr = [[-x for x in row] for row in product(self.mag, de)]
        dew = add(de, product(o, e), -1)
        dPw = add(add(dP, product(o, self.P), -1), product(self.P, o))
        drw = add(dr, product(w, self.r), -1)
        eta = add(de, product(dP, inverse(self.P), e), -1)
        etaw = add(dew, product(dPw, inverse(self.P), e), -1)
        self.assertEqual(etaw, add(eta, product(self.P, o, inverse(self.P), e), -1))
        result = information_shear_correction(self.P, self.mag, self.noise, e,
            self.r, dPw, dew, zeros(3, 21), zeros(3, 3), drw)
        charge = magnetic_connection_charge(self.P, self.mag, self.noise, self.r, e, omega)
        self.assertEqual(result['linked_row_noise_chart_force'], charge['port'])
        self.assertEqual(charge['energy'], charge['linked_energy'])

    def test_row_port_charge_is_same_operation_fisher_loss(self):
        dH, _ = row_differentials([F(1, 5), F(-1, 7), F(1, 11)], [0, 0, 0])
        result = covariance_row_charge(self.P, self.acc, self.noise, dH)
        self.assertLessEqual(result['Fisher_charge'], result['linked_upper'])

    def test_radial_reduction_is_exact_in_planar_stratum(self):
        b = [7, 0, F(1, 9)]
        r = [[F(1, 8)], [F(0)], [F(-1, 5)]]
        de = zeros(21, 1)
        de[1][0] = F(2, 7)
        # Block-diagonal parity ensures the planar output has only the pitch row.
        P = identity(21)
        c = correction(P, self.mag, self.noise, r)
        kappa = sum(b[i]*r[i][0] for i in range(3))/sum(x*x for x in b)
        lhs = product(c['K'], skew([0, de[1][0], 0]), r)
        rhs = [[kappa*x for x in row] for row in product(c['K'], self.mag, de)]
        self.assertEqual(lhs, rhs)

    def test_rank_one_balance_and_square_charge(self):
        P = [[F(2), F(1, 3)], [F(1, 3), F(1)]]
        dP = [[F(1, 5), F(-2, 7)], [F(-2, 7), F(1, 4)]]
        e, eta = [[F(2, 3)], [F(-1, 5)]], [[F(1, 7)], [F(3, 8)]]
        for kappa in (F(3, 500), F(-3, 500)):
            result = rank_one_magnetic_balance(P, [[F(1), F(0)]], F(1), e, eta, dP, kappa)
            self.assertLessEqual(result['final'], result['initial']-F(9, 10)*(
                result['mean_loss']+result['covariance_loss']))
        m = planar_magnetic_loss_margin(F(3, 500), 6, 100)
        self.assertEqual(m['strict_two_by_two_determinant'], '4673/1250000')
        self.assertEqual(m['covariance_square_charge_upper'], '189/10000')
        self.assertTrue(m['verified_on_stated_domain'])
        failed = planar_magnetic_loss_margin(F(1), 6, 100)
        self.assertFalse(failed['verified_on_stated_domain'])
        self.assertIsNone(failed['retained_fraction'])

    def test_two_observation_reader_cancels_inherited_root(self):
        delta, tail = F(7, 5), F(4, 7)
        rows = [[F(1), F(0)], [F(1), delta]]
        reader = [[-tail/delta, 1+tail/delta], [-1/delta, 1/delta]]
        self.assertEqual(product(reader, rows), [[F(1), delta+tail], [F(0), F(1)]])
        bounds = planar_pitch_domain_bounds()
        self.assertEqual(bounds['reader_entry_time_seconds'], '753/250')
        self.assertLess(F(bounds['pitch_variance_upper_exact']), F(3, 5000))
        self.assertLess(F(bounds['BG_y_variance_upper_exact']), F(1, 4000))
        self.assertLess(F(bounds['radial_fraction_abs_upper_exact']), F(3, 500))
        self.assertFalse(bounds['initial_covariance_upper_assumed'])

    def test_certificate_does_not_promote_operation_loss_to_word_gap(self):
        c = certificate()
        self.assertIsNone(c['uniform_complete_word_gap'])
        self.assertFalse(c['all_time_planar_magnetic_service_verified'])
        self.assertFalse(c['theorem_closed'])
        self.assertEqual(c['planar_magnetic_loss']['retained_fraction'], '9/10')

    def test_literal_profile_and_prediction_bindings(self):
        root = Path(__file__).resolve().parents[2]
        core = (root/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h').read_text()
        config = (root/'src/kalman_common/ProxyStartupFusion.h').read_text()
        wrapper = (root/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h').read_text()
        probe = (root/'tools/stability/ou3_theorem/planar_service_probe.cpp').read_text()
        for literal in ('F_AA.template block<3,3>(0,3) = Bstep;',
                        'const Matrix3 J_att = -skew_symmetric_matrix(v2hat);',
                        'S_mat = Rmag;', 'project_psd_ou_iii<T,6>(Q_AA, T(1e-12))'):
            self.assertIn(literal, core)
        self.assertIn('float b0  = 1e-10f;', config)
        self.assertIn('c.Pq0, c.Pb0, c.b0, c.R_S_noise', wrapper)
        self.assertIn('cfg.sigma_g.setConstant(.00135f);cfg.sigma_m.setConstant(.8f);', probe)
        for literal in ('Matrix3 Delta = aw_covariance_floor_target_ - P_aw;',
                        'evals(i) = std::max(T(0), evals(i));',
                        'Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;'):
            self.assertIn(literal, core)

    def test_aw_shear_fixes_acc_row_and_preserves_full_held_joseph(self):
        aw, g = [F(1, 5), 0, F(-1, 7)], [0, 0, 10]
        L, Linv = aw_shear(aw), aw_shear([-x for x in aw])
        H, Hmag = world_rows([a-b for a, b in zip(aw, g)], [7, 0, F(1, 9)])
        Hfixed, _ = world_rows([-x for x in g], [7, 0, F(1, 9)])
        self.assertEqual(product(L, Linv), identity(21))
        self.assertEqual(product(H, Linv), Hfixed)
        self.assertEqual(product(Hmag, Linv), Hmag)
        # L acts only inside the active block: even arbitrary held mask
        # congruence holds, without claiming optimality for nonzero held cross.
        transformed = product(L, self.P, transpose(L))
        for held in (False, True):
            old = correction(self.P, H, self.noise, self.r, held)
            new = correction(transformed, Hfixed, self.noise, self.r, held)
            self.assertEqual(new['K'], product(L, old['K']))
            self.assertEqual(new['C'], product(L, old['C'], transpose(L)))
        da = [F(1, 11), F(-1, 13), F(2, 17)]
        Gamma = aw_connection(da)
        self.assertEqual(product(Gamma, aw_connection(aw)), zeros(21, 21))
        dH, _ = row_differentials(da, [0, 0, 0])
        self.assertEqual(product(H, Gamma), dH)
        row_connection_coboundary(self.P, H, self.noise, Gamma)

    def test_aw_shear_correction_and_changed_endpoint_connection(self):
        aw, da = [F(1, 5), 0, F(-1, 7)], [F(1, 11), 0, F(2, 17)]
        H, _ = world_rows([aw[0], 0, aw[2]-10], [7, 0, 0])
        dH, _ = row_differentials(da, [0, 0, 0])
        e, de, dp = zeros(21, 1), zeros(21, 1), zeros(21, 21)
        e[1][0], e[15][0], e[18][0] = F(1, 20), F(1, 9), F(-1, 12)
        de[1][0], dp[1][15], dp[15][1] = F(1, 13), F(1, 19), F(1, 19)
        for i in range(3):
            de[15+i][0] = da[i]
        dr = [[-x for x in row] for row in product(H, de)]
        old = information_shear_correction(self.P, H, self.noise, e, self.r,
                                          dp, de, dH, zeros(3, 3), dr)
        op = correction(self.P, H, self.noise, self.r)
        dk = gain_differential(self.P, H, op['K'], inverse(op['S']), dp, dH, zeros(3, 3))
        ep = add(e, product(op['K'], self.r))
        dep = add(add(de, product(dk, self.r)), product(op['K'], dr))
        before = aw_shear_differentials(aw, da, self.P, dp, e, de)
        frozen = aw_shear_differentials(aw, da, op['C'], old['posterior_covariance_tangent'], ep, dep)
        Hfixed, _ = world_rows([0, 0, -10], [7, 0, 0])
        new = correction(before['P'], Hfixed, self.noise, self.r)
        A = add(identity(21), product(new['K'], Hfixed), -1)
        mismatch = add(dr, product(Hfixed, before['de']))
        self.assertEqual(frozen['dP'], product(A, before['dP'], transpose(A)))
        self.assertEqual(frozen['eta'], add(product(A, before['eta']), product(new['K'], mismatch)))
        fixed_row_joint_balance(before['P'], Hfixed, self.noise,
                                before['eta'], before['dP'], mismatch)
        aw_next = [aw[i]+op['increment'][15+i][0] for i in range(3)]
        da_next = [dep[15+i][0] for i in range(3)]
        Lnext = aw_shear(aw_next)
        jump = product(Lnext, aw_shear([-x for x in aw]))
        self.assertEqual(jump, aw_shear([op['increment'][15+i][0] for i in range(3)]))
        after = aw_shear_differentials(aw_next, da_next, op['C'],
                                       old['posterior_covariance_tangent'], ep, dep)
        djump = aw_connection([da_next[i]-da[i] for i in range(3)])
        expected = add(add(product(jump, frozen['dP'], transpose(jump)),
                           product(djump, frozen['P'], transpose(jump))),
                       product(jump, frozen['P'], transpose(djump)))
        self.assertEqual(after['dP'], expected)
        self.assertNotEqual(after['dP'], product(jump, frozen['dP'], transpose(jump)))

    def test_aw_increment_and_ou_prediction_keep_literal_structure(self):
        aw, da = [F(1, 5), 0, F(-1, 7)], [F(1, 11), 0, F(2, 17)]
        L, Linv, Gamma = aw_shear(aw), aw_shear([-x for x in aw]), aw_connection(da)
        increment = zeros(21, 21)
        increment[15][15], increment[17][17] = F(1, 9), F(2, 11)
        self.assertEqual(product(L, increment, transpose(L)), increment)
        self.assertEqual(product(Gamma, increment), zeros(21, 21))
        # Exact block algebra, not a frozen complete derivative or orbit.
        phi, Fmap = F(4, 5), identity(21)
        W, B = self.R, [[F(i == j, 7) for j in range(3)] for i in range(3)]
        for i in range(3):
            Fmap[i][:3], Fmap[i][3:6] = W[i], B[i]
            Fmap[15+i][15+i] = phi
            Fmap[6+i][15+i] = F(2, 7)
        transformed = product(aw_shear([phi*x for x in aw]), Fmap, Linv)
        expected = [[phi*x for x in row] for row in product(skew(aw), add(identity(3), W, -1))]
        self.assertEqual([row[:3] for row in transformed[15:18]], expected)
        bg = [[-phi*x for x in row] for row in product(skew(aw), B)]
        self.assertEqual([row[3:6] for row in transformed[15:18]], bg)
        self.assertNotEqual([row[:3] for row in transformed[6:9]], zeros(3, 3))

    def test_magnetic_loss_domain_survives_aw_shear(self):
        aw = [F(1, 5), 0, F(-1, 7)]
        L = aw_shear(aw)
        P = product(L, self.P, transpose(L))
        self.assertEqual([r[:6] for r in P[:6]], [r[:6] for r in self.P[:6]])
        e = zeros(21, 1)
        e[1][0], e[15][0] = F(1, 20), F(1, 9)
        ep = product(L, e)
        self.assertEqual(product(transpose(ep), inverse(P), ep),
                         product(transpose(e), inverse(self.P), e))

    def test_physical_acc_substitution_keeps_covariance_and_source_cross(self):
        # Formal exact operands for the algebraic identity, not a physical
        # trajectory or an independent domain enclosure of angle/rotation.
        e = zeros(21, 1)
        t, omega, da = F(1, 10), F(1, 7), [F(1, 11), 0, F(-1, 13)]
        e[1][0], e[15][0], e[17][0] = t, F(1, 9), F(-2, 17)
        df, nu = [F(3, 5), 0, F(-1, 8)], [F(1, 19), 0, F(1, 23)]
        P = [r[:] for r in self.P]
        P[1][15] = P[15][1] = F(1, 4)
        out = planar_acc_mismatch_charge(P, e, da, omega, df, nu, self.noise)
        residual_plus_ba = [[-e[15+i][0]+t*df[i]+nu[i]] for i in range(3)]
        Jy = skew([0, 1, 0])
        direct = add(product(Jy, [[t*x] for x in da]),
                     [[omega*x for x in r] for r in product(Jy, residual_plus_ba)], -1)
        self.assertEqual(out['mismatch'], direct)
        self.assertLessEqual(out['curvature_energy'], out['curvature_upper'])
        self.assertNotEqual(out['source_cross'], 0)
        uncoupled = [r[:] for r in P]
        uncoupled[1][15] = uncoupled[15][1] = F(0)
        other = planar_acc_mismatch_charge(uncoupled, e, da, omega, df, nu, self.noise)
        self.assertNotEqual(out['linked_covariance_form'], other['linked_covariance_form'])

    def test_complete_word_connections_cancel_only_with_generated_score(self):
        # Formal rational operands for the universal port-conjugacy identity.
        # Neither this word nor its covariance is asserted shipping-reachable.
        B = [[[F(1), F(1, 3)], [0, F(4, 5)]],
             [[F(2, 3), 0], [F(1, 7), F(1)]]]
        P = [identity(2), [[F(2), F(1, 5)], [F(1, 5), F(3)]], identity(2)]
        e = [[[F(1, 4)], [F(-1, 5)]], [[F(2, 7)], [F(1, 3)]], [[F(1, 9)], [F(2, 5)]]]
        d = [[[F(1, 6)], [F(2, 7)]], [[F(-1, 8)], [F(1, 5)]]]
        U = [[[F(1, 7), F(1, 8)], [F(1, 8), F(-1, 9)]], zeros(2, 2)]
        v = [[[F(1, 11)], [F(2, 13)]], [[F(1, 17)], [F(-1, 19)]]]
        D = [[F(1, 3), F(1, 5)], [F(1, 5), F(-1, 7)]]
        de = [[F(2, 9)], [F(-1, 4)]]
        L = [[[F(1), 0], [F(1, 3), F(1)]],
             [[F(1), F(1, 5)], [0, F(1)]], [[F(2), 0], [F(1, 7), F(1)]]]
        G = [[[0, F(1, 7)], [F(1, 11), 0]],
             [[0, F(1, 13)], [F(1, 17), 0]], [[0, F(1, 19)], [F(1, 23), 0]]]
        old = word_score_normal_form(B, P, e, d, U, v, D, de)
        ports = [moving_frame_word_ports(B[i], P[i], P[i+1], e[i], e[i+1],
                 d[i], U[i], v[i], L[i], L[i+1], G[i], G[i+1]) for i in range(2)]
        connectionD = add(product(G[0], P[0]), product(P[0], transpose(G[0])))
        D0 = product(L[0], add(D, connectionD), transpose(L[0]))
        de0 = product(L[0], add(de, product(G[0], e[0])))
        pp = [product(l, p, transpose(l)) for l, p in zip(L, P)]
        ee = [product(l, x) for l, x in zip(L, e)]
        new = word_score_normal_form([x['base'] for x in ports], pp, ee,
            [x['gain_score'] for x in ports], [x['covariance_port'] for x in ports],
            [x['mean_port'] for x in ports], D0, de0)
        expected_eta = product(L[-1], add(old['eta_terminal'],
            product(P[-1], transpose(G[-1]), inverse(P[-1]), e[-1]), -1))
        expected_D = product(L[-1], add(add(old['covariance_tangent_terminal'],
            product(G[-1], P[-1])), product(P[-1], transpose(G[-1]))), transpose(L[-1]))
        self.assertEqual(new['eta_terminal'], expected_eta)
        self.assertEqual(new['covariance_tangent_terminal'], expected_D)
        # Covariance created at the first step has a nonzero future-score term.
        omitted = product(new['base_suffixes'][1], ports[0]['covariance_port'], new['suffix_scores'][1])
        self.assertNotEqual(omitted, zeros(2, 1))

    def test_aw_boundary_work_exact_completion_and_conditional_precision(self):
        e, eta, D = zeros(21, 1), zeros(21, 1), zeros(21, 21)
        e[1][0], e[17][0], eta[1][0] = F(1, 20), F(1, 9), F(1, 7)
        D[1][17] = D[17][1] = F(1, 11)
        da = [F(1, 13), 0, F(-1, 17)]
        P = [r[:] for r in self.P]
        P[1][17] = P[17][1] = F(1, 3)
        out = aw_frame_boundary_work(P, e, eta, D, da)
        de = add(eta, product(D, inverse(P), e))
        moved = aw_shear_differentials([F(2, 7), 0, F(-1, 5)], da, P, D, e, de)
        Jm = inverse(moved['P'])
        from tools.stability.ou3_theorem.information_shear_word import trace
        actual = product(transpose(moved['eta']), Jm, moved['eta'])[0][0]
        actual += trace(product(Jm, moved['dP'], Jm, moved['dP']))
        self.assertEqual(actual, out['transformed_storage'])
        zero = aw_frame_boundary_work(P, e, eta, D, [0, 0, 0])
        self.assertEqual(zero['signed_boundary_work'], 0)
        # Q is controlled by conditional precision, not the AW marginal inverse.
        J = inverse(P)
        self.assertNotEqual([r[15:18] for r in J[15:18]],
                            inverse([r[15:18] for r in P[15:18]]))
        ja = product(J, e)[15:18]
        jaa = [r[15:18] for r in J[15:18]]
        ptt = [r[:3] for r in P[:3]]
        B = [r[15:18] for r in product(D, J)[:3]]
        V = add(product(ja, transpose(ja)), [[2*x for x in r] for r in jaa])
        for j in range(3):
            Sj = skew([F(i == j) for i in range(3)])
            hj = product(transpose(ja), Sj, eta[:3])[0][0]-2*trace(product(B, Sj))
            self.assertEqual(hj, out['linear_coefficient'][j][0])
            for k in range(3):
                Sk = skew([F(i == k) for i in range(3)])
                self.assertEqual(trace(product(V, Sj, ptt, transpose(Sk))),
                                 out['quadratic_coefficient'][j][k])
        self.assertGreaterEqual(out['completed_square_lower'], -out['initial_storage'])

    def test_actual_aw_increment_decreases_conditional_connection_coefficient(self):
        e = zeros(21, 1)
        e[1][0], e[15][0], e[17][0] = F(1, 7), F(1, 9), F(-1, 11)
        # A rank-deficient PSD increment checks the zero-face scope. Its
        # reachability/target is NOT inferred from this identity regression.
        increment = [[F(1, 5), 0, F(1, 10)], [0, 0, 0], [F(1, 10), 0, F(1, 20)]]
        out = aw_conditional_sync(self.P, e, increment, [F(1, 3), 0, F(1, 7)])
        Pnext = [r[:] for r in self.P]
        for i in range(3):
            for j in range(3):
                Pnext[15+i][15+j] += increment[i][j]
        self.assertEqual(out['next_conditional_precision'],
                         [r[15:18] for r in inverse(Pnext)[15:18]])
        drop = product(transpose(e), add(inverse(self.P), inverse(Pnext), -1), e)[0][0]
        self.assertEqual(drop, out['conditional_comparison_energy']-out['next_conditional_comparison_energy'])
        self.assertGreater(out['Fisher_connection_square_decrease'], 0)
        self.assertEqual(aw_conditional_sync(self.P, e, zeros(3, 3), [1, 0, 0])[
                         'Fisher_connection_square_decrease'], 0)

    def test_mag_S_conditional_storage_invariant_with_row_noise_ports(self):
        e, de, D = zeros(21, 1), zeros(21, 1), zeros(21, 21)
        e[1][0], e[15][0], e[12][0] = F(1, 7), F(2, 9), F(1, 5)
        de[1][0], de[15][0] = F(1, 11), F(2, 13)
        D[1][15] = D[15][1] = F(1, 17)
        eta = add(de, product(D, inverse(self.P), e), -1)
        old = conditional_aw_storage(self.P, e, eta, D)
        Hs = zeros(3, 21)
        for i in range(3):
            Hs[i][12+i] = F(1)
        for H in (self.mag, Hs):
            # dH_aw=0: arbitrary linked non-AW row/noise/residual ports are
            # allowed by the identity. These operands are not a reached cell.
            dH, dR = zeros(3, 21), zeros(3, 3)
            dH[0][1], dR[0][0] = F(1, 19), F(1, 23)
            dr = [[F(1, 29)], [F(-1, 31)], [F(2, 37)]]
            op = correction(self.P, H, self.noise, self.r)
            out = information_shear_correction(self.P, H, self.noise, e, self.r,
                D, de, dH, dR, dr)
            dk = gain_differential(self.P, H, op['K'], inverse(op['S']), D, dH, dR)
            ep = add(e, op['increment'])
            dep = add(add(de, product(dk, self.r)), product(op['K'], dr))
            Dp = out['posterior_covariance_tangent']
            etap = add(dep, product(Dp, inverse(op['C']), ep), -1)
            new = conditional_aw_storage(op['C'], ep, etap, Dp)
            for key in ('regression', 'conditional_covariance', 'd_regression',
                        'd_conditional_covariance', 'conditional_information',
                        'conditional_score', 'conditional_storage'):
                self.assertEqual(old[key], new[key], key)
            beforeQ = aw_frame_boundary_work(self.P, e, eta, D, [0, 0, 0])['quadratic_coefficient']
            afterQ = aw_frame_boundary_work(op['C'], ep, etap, Dp, [0, 0, 0])['quadratic_coefficient']
            from tools.stability.ou3_theorem.matrix_certificates import is_psd
            self.assertTrue(is_psd(add(beforeQ, afterQ, -1)))

    def test_literal_OU_row_keeps_integrated_chain_in_precision_balance(self):
        phi, Fmap, Q = F(4, 5), identity(21), identity(21)
        for i in range(3):
            Fmap[15+i][15+i] = phi
            Fmap[i][3+i] = F(1, 9)
            for row, coefficient in ((6, F(2, 7)), (9, F(1, 11)), (12, F(1, 13))):
                Fmap[row+i][15+i] = coefficient
        out = aw_precision_process_balance(self.P, Fmap, Q, phi, [F(1, 7), 0, F(-1, 9)])
        self.assertGreater(out['positive_connection_loss'], 0)
        self.assertNotEqual(out['signed_chain_work'], 0)
        self.assertNotEqual(out['signed_AG_work'], 0)
        self.assertEqual(out['connection_after']-out['connection_before'],
            -out['positive_connection_loss']+out['signed_chain_work']+out['signed_AG_work'])
        Fmap[15][6] = F(1, 17)
        with self.assertRaises(ValueError):
            aw_precision_process_balance(self.P, Fmap, Q, phi, [1, 0, 0])

    def test_existing_path_certificate_bounds_conditional_AW_precision(self):
        out = qualified_aw_precision_ceiling()
        self.assertLess(F(out['exact_upper']), F(15942618))
        self.assertTrue(out['conditional_not_marginal_precision'])
        self.assertFalse(out['sufficient_for_endpoint_Schur_margin'])

    def test_conditional_AW_process_loss_is_positive_but_not_full_gap(self):
        out = conditional_aw_process_coercivity()
        q, m = F(out['shorted_process_noise_lower']), F(out['conditional_covariance_upper'])
        c = F(out['mean_loss_fraction_lower'])
        self.assertEqual(c, q/(m+q))
        self.assertGreater(c, F(1, 10**14))
        self.assertEqual(F(out['covariance_loss_fraction_lower']), 2*c-c*c)
        self.assertFalse(out['mixed_block_or_signed_work_absorption'])
        self.assertFalse(out['uniform_complete_word_gap'])

    def test_conditional_process_comparison_retains_full_cross_covariance(self):
        from tools.stability.ou3_theorem.information_shear_word import trace
        from tools.stability.ou3_theorem.matrix_certificates import is_psd
        # Formal exact operands verify CA8/CA9; no reached-domain claim.
        E = zeros(21, 3)
        for i in range(3):
            E[15+i][i] = F(1)
        Fmap = identity(21)
        for i in range(3):
            Fmap[15+i][15+i] = F(4, 5)
            Fmap[6+i][15+i], Fmap[9+i][15+i] = F(1, 7), F(1, 13)
            Fmap[12+i][15+i] = F(1, 19)
        J = inverse(self.P)
        C = inverse(product(transpose(E), J, E))
        q = F(1, 20)
        noise = add(product(Fmap, E, [[q*x for x in r] for r in identity(3)],
                            transpose(E), transpose(Fmap)), identity(21))
        Jnext = inverse(add(product(Fmap, self.P, transpose(Fmap)), noise))
        compressed = product(transpose(E), transpose(Fmap), Jnext, Fmap, E)
        upper = inverse(add(C, [[q*x for x in r] for r in identity(3)]))
        self.assertTrue(is_psd(add(upper, compressed, -1)))
        m = max(sum(abs(x) for x in r) for r in C)
        c = q/(m+q)
        sa = [[F(1, 7)], [F(-1, 11)], [F(2, 13)]]
        dC = [[F(1, 17), F(1, 19), 0], [F(1, 19), F(-1, 23), 0], [0, 0, F(1, 29)]]
        eta, D = product(E, C, sa), product(E, dC, transpose(E))
        mean0 = product(transpose(eta), J, eta)[0][0]
        mean1 = product(transpose(eta), transpose(Fmap), Jnext, Fmap, eta)[0][0]
        fisher0 = trace(product(J, D, J, D))
        Dnext = product(Fmap, D, transpose(Fmap))
        fisher1 = trace(product(Jnext, Dnext, Jnext, Dnext))
        self.assertGreaterEqual(mean0-mean1, c*mean0)
        self.assertGreaterEqual(fisher0-fisher1, (2*c-c*c)*fisher0)

    def test_conditional_mixed_reader_is_invariant_under_state_dependent_AW_frame(self):
        e, de, D = zeros(21, 1), zeros(21, 1), zeros(21, 21)
        e[1][0], e[17][0], de[1][0], de[15][0] = F(1, 7), F(1, 11), F(1, 13), F(1, 17)
        D[1][17] = D[17][1] = F(1, 19)
        aw, daw = [F(1, 5), 0, F(1, 7)], [F(1, 17), 0, F(-1, 23)]
        moved = aw_shear_differentials(aw, daw, self.P, D, e, de)
        eta = add(de, product(D, inverse(self.P), e), -1)
        original = conditional_aw_storage(self.P, e, eta, D)
        changed = conditional_aw_storage(moved['P'], moved['e'], moved['eta'], moved['dP'])
        for key in ('conditional_covariance', 'd_conditional_covariance',
                    'conditional_score', 'conditional_storage'):
            self.assertEqual(original[key], changed[key], key)
        boundary = aw_frame_boundary_work(self.P, e, eta, D, daw)
        self.assertEqual(changed['complementary_storage']-original['complementary_storage'],
                         boundary['signed_boundary_work'])
        # Endpoint congruence of the BASE word; state-dependent derivatives
        # above remain in the actual complementary packets, never frozen.
        E = zeros(21, 3)
        for i in range(3):
            E[15+i][i] = F(1)
        M, PN = identity(21), add(self.P, identity(21))
        LN, L0 = aw_shear([F(1, 3), F(1, 11), F(-1, 7)]), moved['L']
        Mt, PNt = product(LN, M, inverse(L0)), product(LN, PN, transpose(LN))
        q = zeros(21, 1)
        q[15][0], q[1][0] = F(1, 101), F(1, 103)
        qt = product(transpose(inverse(L0)), q)
        L, Lt = product(M, E), product(Mt, E)
        before = conditional_mixed_coefficients(original['conditional_covariance'],
            product(transpose(L), inverse(PN), L), product(transpose(E), q), F(2))
        after = conditional_mixed_coefficients(changed['conditional_covariance'],
            product(transpose(Lt), inverse(PNt), Lt), product(transpose(E), qt), F(2))
        self.assertEqual(before, after)


if __name__ == '__main__':
    unittest.main()
