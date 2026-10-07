import unittest
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_theorem.lin_path_certificate import inverse
from tools.stability.ou3_theorem.matrix_certificates import add, identity, transpose
from tools.stability.ou3_theorem.planar_innovation_storage import (
    certificate, energy, held_ba_reachable_correction, information_shear_correction,
    information_shear_prediction, inherited_handoff_bounds, innovation_domain_bound,
    innovation_variational_identity, linked_prediction_charge,
    physical_acceleration_remainder_form, shear_coercivity_factor,
)
from tools.stability.ou3_theorem.planar_linked_riccati_mean import product


class InnovationStorageTests(unittest.TestCase):
    """Exact rational algebra regressions, never reachable-word experiments."""

    def setUp(self):
        self.P = [[F(2), F(1, 3)], [F(1, 3), F(1)]]
        self.H = [[F(1), F(-1, 2)], [F(1, 4), F(2)]]
        self.R = [[F(3, 5), F(1, 7)], [F(1, 7), F(1)]]
        self.dP = [[F(1, 5), F(-2, 7)], [F(-2, 7), F(-1, 3)]]
        self.e = [[F(3, 7)], [F(-2, 5)]]
        self.de = [[F(2, 11)], [F(1, 3)]]
        self.r = [[F(2, 9)], [F(-3, 8)]]
        self.zero = [[F(0), F(0)], [F(0), F(0)]]

    def test_nis_has_exact_square_gap_without_triangle_cross_charge(self):
        result = innovation_variational_identity(self.P, self.H, self.R, self.r, self.e)
        self.assertGreater(result['square_gap'], 0)
        optimal = innovation_variational_identity(self.P, self.H, self.R, self.r, result['minimizer'])
        self.assertEqual(optimal['square_gap'], 0)
        self.assertEqual(optimal['NIS'], optimal['comparison_upper'])

    def test_fixed_rows_cancel_covariance_gain_residual_for_nonzero_residual(self):
        dr = [[-x for x in row] for row in product(self.H, self.de)]
        result = information_shear_correction(self.P, self.H, self.R, self.e,
            self.r, self.dP, self.de, self.zero, self.zero, dr)
        self.assertEqual(result['eta_plus'], result['homogeneous'])
        # Comparison reseeding to zero at the output would destroy the identity.
        self.assertNotEqual(self.r, [[F(0)], [F(0)]])

    def test_variable_rows_noise_and_residual_chart_port_are_retained(self):
        dH = [[F(1, 5), F(1, 2)], [F(-1, 7), F(2, 9)]]
        dR = [[F(1, 8), F(-1, 11)], [F(-1, 11), F(1, 4)]]
        result = information_shear_correction(self.P, self.H, self.R, self.e,
            self.r, self.dP, self.de, dH, dR, [[F(3, 5)], [F(-2, 3)]])
        self.assertNotEqual(result['eta_plus'], result['homogeneous'])

    def test_shear_is_coercive_only_with_comparison_domain(self):
        J = inverse(self.P)
        V = energy(self.e, J)
        weight = F(2, 3)
        eta = add(self.de, product(self.dP, J, self.e), -1)
        Z2 = product(J, self.dP, J, self.dP)
        cov = sum((row[i] for i, row in enumerate(Z2)), F(0))
        base = energy(self.de, J)+weight*cov
        shear = energy(eta, J)+weight*cov
        C = shear_coercivity_factor(V, weight)
        self.assertLessEqual(base/C, shear)
        self.assertLessEqual(shear, C*base)

    def test_prediction_retains_source_and_coefficient_variation(self):
        Fmap = [[F(1), F(1, 5)], [F(0), F(3, 4)]]
        source = [[F(1, 7)], [F(-1, 3)]]
        result = information_shear_prediction(self.P, Fmap, self.R, self.e,
            source, self.dP, self.H, self.dP, self.de, source)
        self.assertNotEqual(result['eta_plus'], result['homogeneous'])
        c = linked_prediction_charge(self.P, Fmap, self.R, self.e, self.dP)
        self.assertGreater(c['mean_process_loss'], 0)
        self.assertLessEqual(c['sheared_covariance_port_energy'], c['linked_charge_upper'])

    def test_held_ba_matches_literal_masked_joseph_without_changing_ba(self):
        # One active and one BA coordinate suffice to test the block identity.
        px, pb, h, r = [[F(2)]], [[F(3, 7)]], [[F(3)]], [[F(1, 5)]]
        result = held_ba_reachable_correction(px, pb, h, r)
        P = [[px[0][0], F(0)], [F(0), pb[0][0]]]
        Hs = [[F(3), F(1)]]
        S = add(product(Hs, P, transpose(Hs)), r)
        B = [[F(6)], [F(0)]]  # literal masked PCt, not P Hs'
        K = product(B, inverse(S))
        C = add(add(add(P, product(K, transpose(B)), -1),
                    product(B, transpose(K)), -1), product(K, S, transpose(K)))
        self.assertEqual(S, result['innovation'])
        self.assertEqual(C, [[result['active_C'][0][0], F(0)], [F(0), pb[0][0]]])
        self.assertNotEqual(K, product(P, transpose(Hs), inverse(S)))

    def test_source_contains_held_decoupling_and_no_held_ba_process_noise(self):
        source = (Path(__file__).resolve().parents[2]/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h').read_text()
        hold = source.split('void set_acc_bias_updates_enabled(bool en)')[1].split('// Velocity in world')[0]
        for block in ('<3,BASE_N>(OFF_BA, 0)', '<BASE_N,3>(0, OFF_BA)',
                      '<3,12>(OFF_BA, OFF_V)', '<12,3>(OFF_V, OFF_BA)'):
            self.assertIn(block+'.setZero()', hold)
        prediction = source.split('// Optional residual accel-bias OU and cross terms')[1].split('Eigen::Matrix<T,NA,NB> tmpAB')[0]
        self.assertIn('acc_bias_updates_enabled_ ? std::exp(-Ts / tau_b) : T(1)', prediction)
        self.assertIn('if (acc_bias_updates_enabled_) {', prediction)
        self.assertIn('P_BB.noalias() += Q_bacc_ * qd_scale;', prediction)
        self.assertIn('if (!use_ba) freeze_acc_bias_rows_(PCt);', source)
        self.assertIn('if (!use_ba) freeze_acc_bias_rows_(K);', source)

    def test_remainder_domination_fails_closed_and_status_stays_open(self):
        B = physical_acceleration_remainder_form(21, [0, 1, 2], [15, 16, 17], F(10), F(1))
        result = innovation_domain_bound(identity(21), B, F(6), F(1, 100), F(0), F(1, 25))
        self.assertEqual(result, F(1, 10))
        with self.assertRaisesRegex(ValueError, 'D_SUFFICIENT_BOUND_FAILURE'):
            innovation_domain_bound(identity(21), B, F(1), F(1, 100), F(0), F(1, 25))
        self.assertIsNone(certificate()['uniform_complete_word_epsilon'])
        self.assertIsNone(certificate()['uniform_NIS_cap'])
        self.assertFalse(certificate()['source_uniform_nominal_domain_verified'])

    def test_first_handoff_bound_is_exact_and_is_not_a_prefix_certificate(self):
        result = inherited_handoff_bounds()
        self.assertLess(F(result['precision_energy_upper_exact']), F(21, 1000))
        self.assertLess(F(result['accelerometer_NIS_upper_exact']), F(21, 1000))
        self.assertLess(F(result['magnetic_NIS_upper_exact']), F(13, 200))
        self.assertGreater(F(result['precision_energy_terms']['physical_S']), 0)
        self.assertGreater(F(result['precision_energy_terms']['slow_BA_central_fibre']), 0)
        self.assertFalse(result['post_handoff_prefix_invariance_verified'])
        self.assertFalse(result['finite_handoff_time_verified'])


if __name__ == '__main__':
    unittest.main()
