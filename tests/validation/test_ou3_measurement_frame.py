"""Exact rational identities, not trajectory or contraction experiments."""
from fractions import Fraction as F
from pathlib import Path
import unittest

from tools.stability.ou3_theorem.information_shear_word import zeros, word_score_normal_form
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
)
from tools.stability.ou3_theorem.planar_innovation_storage import information_shear_correction
from tools.stability.ou3_theorem.planar_linked_riccati_mean import product
from tools.stability.ou3_theorem.planar_linked_riccati_mean import gain_differential
from tools.stability.ou3_theorem.world_frame import quaternion_rotation, skew


class MeasurementFrameTests(unittest.TestCase):
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


if __name__ == '__main__':
    unittest.main()
