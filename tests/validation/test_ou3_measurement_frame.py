"""Exact rational identities, not trajectory or contraction experiments."""
from fractions import Fraction as F
from pathlib import Path
import unittest

from tools.stability.ou3_theorem.information_shear_word import zeros
from tools.stability.ou3_theorem.lin_path_certificate import inverse
from tools.stability.ou3_theorem.matrix_certificates import add, identity, transpose
from tools.stability.ou3_theorem.measurement_frame import (
    certificate, connection, correction, covariance_row_charge, frame,
    magnetic_connection_charge, planar_magnetic_loss_margin, planar_pitch_domain_bounds,
    pullback_differentials, rank_one_magnetic_balance,
    row_differentials, world_rows,
)
from tools.stability.ou3_theorem.planar_innovation_storage import information_shear_correction
from tools.stability.ou3_theorem.planar_linked_riccati_mean import product
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


if __name__ == '__main__':
    unittest.main()
