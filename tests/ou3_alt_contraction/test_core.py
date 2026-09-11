"""Host algebra and anti-shortcut tests; none is a universal stability gate."""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_alt_contraction.core import (
    measurement_lift, dissipation_form, sector_iqc, iqc_master,
    structural_findings, open_obligations, bias_prediction_lift, joint_bias_projection_iqc,
    measurement_increment_lift, masked_supply_certificate,
)


class WordAlgebraTests(unittest.TestCase):
    def test_masked_gain_is_not_replaced_with_PHt(self):
        H, N, S = np.array([[1., 1.]]), np.array([[1.], [0.]]), np.array([[3.]])
        lift = measurement_lift(H, N, S)
        e, w = np.array([2., -1.]), np.array([.25])
        q = np.linalg.solve(S, H @ e + w)
        chi = np.concatenate((e, w, q))
        np.testing.assert_allclose(lift.equality @ chi, 0, atol=1e-14)
        np.testing.assert_allclose(lift.after @ chi, e-N @ q)
        self.assertEqual((lift.after @ chi)[1], e[1])

    def test_dense_unmasked_graph_matches_solve(self):
        rng = np.random.default_rng(984)
        A = rng.normal(size=(6, 6))
        P = A @ A.T + np.eye(6)
        H = rng.normal(size=(3, 6))
        N, S = P @ H.T, H @ P @ H.T + 2*np.eye(3)
        e, w = rng.normal(size=6), rng.normal(size=3)
        lift = measurement_lift(H, N, S)
        q = np.linalg.solve(S, H @ e + w)
        chi = np.r_[e, w, q]
        np.testing.assert_allclose(lift.equality @ chi, 0, atol=1e-12)
        np.testing.assert_allclose(lift.after @ chi, e-N @ q)

    def test_dissipation_telescopes_on_graph(self):
        lift = measurement_lift([[1]], [[2]], [[5]])
        Q = dissipation_form(lift.before, lift.after, lift.disturbance,
                             [[1]], [[1]], .5, [[.6]])
        for e, w in [(1., 0.), (.7, -.2), (0., 1.)]:
            chi = np.array([e, w, (e+w)/5])
            expected = (.6*e-.4*w)**2-.5*e*e-.6*w*w
            self.assertAlmostEqual(float(chi @ Q @ chi), expected)
            self.assertLess(expected, 0)

    def test_joint_cross_terms_and_bounded_supply_retained(self):
        M = np.array([[2., .4], [.4, 1.]])
        X0, X1 = np.eye(2), np.diag([.5, 1.])
        Q = dissipation_form(X0, X1, np.zeros((0, 2)), M, M, .9,
                             np.zeros((0, 0)), [[0., 1.]], [[.3]])
        expected = X1.T @ M @ X1-.9*M-np.diag([0., .3])
        np.testing.assert_allclose(Q, expected)
        self.assertNotEqual(Q[0, 1], 0)

    def test_equality_multiplier_vanishes_on_graph(self):
        lift = measurement_lift([[1]], [[2]], [[5]])
        Q = np.diag([1., -2., 3.])
        Y = np.array([[.1], [-.7], [.3]])
        master = iqc_master(Q, lift.equality, Y)
        chi = np.array([.8, -.4, .08])
        self.assertAlmostEqual(float(chi @ master @ chi), float(chi @ Q @ chi))

    def test_projection_iqc_sign(self):
        Q = sector_iqc(3)
        def project(x):
            return x/max(1., np.linalg.norm(x))
        rng = np.random.default_rng(83)
        for _ in range(100):
            a, b = rng.normal(size=3)*3, rng.normal(size=3)*3
            dx, dy = a-b, project(a)-project(b)
            chi = np.r_[dx, dy]
            self.assertGreaterEqual(float(chi @ Q @ chi), -1e-12)

    def test_s_procedure_uses_plus_sign(self):
        # q>=0 and Q+q<=0 imply Q<=0. Minus would certify Q=q>0 incorrectly.
        iqc = sector_iqc(1)
        E, Y = np.zeros((0, 2)), np.zeros((2, 0))
        np.testing.assert_allclose(iqc_master(-iqc, E, Y, [(1., iqc)]), 0)
        chi = np.array([1., .5])
        self.assertGreater(float(chi @ iqc @ chi), 0)
        self.assertGreater(float(chi @ iqc_master(iqc, E, Y, [(1., iqc)]) @ chi), 0)

    def test_negative_multiplier_rejected(self):
        with self.assertRaises(ValueError):
            iqc_master(np.eye(2), np.zeros((0, 2)), np.zeros((2, 0)), [(-1, np.eye(2))])

    def test_invalid_inputs_rejected(self):
        with self.assertRaises(np.linalg.LinAlgError):
            measurement_lift([[1]], [[1]], [[-1]])
        with self.assertRaises(ValueError):
            measurement_lift([[1]], [[float('nan')]], [[1]])
        with self.assertRaises(ValueError):
            sector_iqc(1, 2, 1)

    def test_exact_neutral_obstruction(self):
        d = structural_findings()
        self.assertTrue(d['held_mode_neutral_witness']['A_v_equals_v'])
        self.assertFalse(d['masked_information_identity']['unqualified_identity_valid'])
        self.assertEqual(d['masked_information_identity']['actual_minus_unmasked_information'],
                         [['-1/2', '-1'], ['-1', '-1']])

    def test_joint_bias_driver_is_shared(self):
        for ph in (1., .9999):
            for pt in (1., .9998):
                F21 = np.eye(21)
                F21[18:21, 18:21] *= ph
                A, B = bias_prediction_lift(F21, ph, pt)
                z = np.arange(24, dtype=float)/100
                w = np.array([.001, -.002, .003])
                nxt = A @ z + B @ w
                np.testing.assert_allclose(nxt[18:21], ph*z[18:21]+(pt-ph)*z[21:]+w)
                np.testing.assert_allclose(nxt[21:], pt*z[21:]+w)
                np.testing.assert_allclose(B[18:21], B[21:])

    def test_projection_retains_true_bias_cross_information(self):
        Q = joint_bias_projection_iqc()
        def projected_error(e, beta):
            estimate = beta-e
            return beta-estimate/max(1., np.linalg.norm(estimate)/.4)
        rng = np.random.default_rng(74)
        for _ in range(100):
            e1, e2, b1, b2 = rng.normal(size=(4, 3))
            dy = projected_error(e1, b1)-projected_error(e2, b2)
            chi = np.r_[e1-e2, b1-b2, dy]
            self.assertGreaterEqual(float(chi @ Q @ chi), -1e-12)
        self.assertGreater(np.linalg.norm(Q[:3, 3:6]), 0)

    def test_finite_gain_increment_retains_nonzero_residual(self):
        rng = np.random.default_rng(176)
        for _ in range(10):
            B0, B1 = rng.normal(size=(2, 2)), rng.normal(size=(2, 2))
            S0, S1 = B0 @ B0.T+np.eye(2), B1 @ B1.T+np.eye(2)
            N0, N1 = rng.normal(size=(4, 2)), rng.normal(size=(4, 2))
            r0, r1 = rng.normal(size=2), rng.normal(size=2)
            q0, q1 = np.linalg.solve(S0, r0), np.linalg.solve(S1, r1)
            lift = measurement_increment_lift(N1, S1, q0)
            chi = np.r_[r1-r0, (N1-N0).ravel(), (S1-S0).ravel(), q1-q0]
            np.testing.assert_allclose(lift.equality @ chi, 0, atol=1e-12)
            np.testing.assert_allclose(lift.correction_difference @ chi,
                                       N1 @ q1-N0 @ q0, atol=1e-12)
            frozen = N1 @ np.linalg.solve(S1, r1-r0)
            self.assertGreater(np.linalg.norm(frozen-(N1 @ q1-N0 @ q0)), 1e-6)

    def test_positive_exact_masked_supply_certificate(self):
        d = masked_supply_certificate()
        self.assertTrue(d['negative_definite_exact'])
        self.assertEqual(d['Q_determinant'], '1027/2880')
        self.assertEqual(d['M_determinant'], '15/16')
        self.assertFalse(d['shipping_certificate'])

    def test_no_proof_promotion(self):
        status = open_obligations()
        for key in ('ALT_LIVE_PASS', 'ALT_STARTUP_PASS', 'ALT_END_TO_END_PASS',
                    'P4_promoted', 'P5_promoted'):
            self.assertIs(status[key], False)
        self.assertGreaterEqual(len(status['missing']), 8)


if __name__ == '__main__':
    unittest.main()
