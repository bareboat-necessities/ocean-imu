"""Exact equality for masked/unmasked Joseph: no optimal-gain assumptions."""
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(Path(__file__).parent)]
from test_finite_measurement_graph import Poly, solve3, dense_covariance
from tools.stability.ou3_alt_contraction import finite_measurement_graph as G


class FiniteCovarianceRank3Tests(unittest.TestCase):
    def test_cancellation_for_all_independent_row_and_symmetric_core_coefficients(self):
        # One arbitrary i,j covariance entry is enough: every matrix entry is
        # this same polynomial with independent row variables, not point tests.
        with patch.object(Poly, 'N', 12):
            x = [Poly.var(i) for i in range(12)]
            ki, kj = x[:3], x[3:6]
            S = [[x[6], x[7], x[8]], [x[7], x[9], x[10]], [x[8], x[10], x[11]]]
            ni, nj = G.mv(S, ki), G.mv(S, kj)
            raw_delta = -G.dot(ki, nj)-G.dot(ni, kj)+G.dot(ki, G.mv(S, kj))
            self.assertEqual(raw_delta, -G.dot(ki, nj))

    def test_held_and_active_full_covariance_keep_actual_numerator(self):
        P = dense_covariance()
        for active in (False, True):
            H, _ = G.residual_factors('accelerometer', [0]*3,
                                      f_hat=[F(1, 5), F(-1, 4), -9], R_hat=G.eye(3))
            N, S = G.measurement_operands(P, G.eye(3), H, active_bias=active)
            K = [solve3(S, row) for row in N]
            raw = G.plus(G.plus(G.plus(P, G.mm(K, G.transpose(N)), -1),
                                G.mm(N, G.transpose(K)), -1), G.mm(G.mm(K, S), G.transpose(K)))
            self.assertEqual(G.solved_joseph_covariance(P, N, S, K), raw)
            if not active:
                self.assertNotEqual(N, G.mm(P, G.transpose(H)))
                self.assertEqual(N[18:21], G.zeros(3, 3))

    def test_numerical_solve_defect_cannot_be_cancelled(self):
        P, N, S = G.eye(21), G.zeros(21, 3), G.eye(3)
        K = G.zeros(21, 3); K[0][0] = F(1, 10**12)
        with self.assertRaisesRegex(ValueError, 'gain relation'):
            G.solved_joseph_covariance(P, N, S, K)


if __name__ == '__main__':
    unittest.main()
