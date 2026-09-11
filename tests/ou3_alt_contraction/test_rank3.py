import unittest
import numpy as np
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_alt_contraction import rank3


class Rank3Tests(unittest.TestCase):
    def test_state_apply_matches_dense(self):
        rng = np.random.default_rng(17)
        for n in (18, 21, 24):
            Psi = rng.normal(size=(n, n))
            K = rng.normal(size=(n, 3))
            H = rng.normal(size=(3, n))
            np.testing.assert_allclose(
                rank3.state_apply(Psi, K, H),
                (np.eye(n) - K @ H) @ Psi,
                rtol=2e-13, atol=2e-13,
            )

    def test_joseph_matches_dense(self):
        rng = np.random.default_rng(23)
        n = 24
        X = rng.normal(size=(n, n))
        P = X @ X.T + np.eye(n)
        H = rng.normal(size=(3, n))
        R = np.diag([.7, .9, 1.1])
        PCt = P @ H.T
        S = H @ PCt + R
        K = np.linalg.solve(S, PCt.T).T
        A = np.eye(n) - K @ H
        np.testing.assert_allclose(
            rank3.joseph(P, K, S, PCt),
            A @ P @ A.T + K @ R @ K.T,
            rtol=2e-12, atol=2e-12,
        )

    def test_storage_delta_rank_at_most_six(self):
        rng = np.random.default_rng(29)
        n = 24
        X = rng.normal(size=(n, n))
        M = X @ X.T + np.eye(n)
        K = rng.normal(size=(n, 3))
        H = rng.normal(size=(3, n))
        A = np.eye(n) - K @ H
        d = rank3.storage_delta(M, K, H)
        np.testing.assert_allclose(d, A.T @ M @ A - M, rtol=2e-12, atol=2e-12)
        self.assertLessEqual(np.linalg.matrix_rank(d, tol=1e-9), 6)

    def test_product_ports_match_exact_finite_increment(self):
        rng = np.random.default_rng(41)
        n, m = 24, 3
        B0, B1 = rng.normal(size=(m, m)), rng.normal(size=(m, m))
        S0, S1 = B0 @ B0.T + np.eye(m), B1 @ B1.T + np.eye(m)
        N0, N1 = rng.normal(size=(n, m)), rng.normal(size=(n, m))
        r0, r1 = rng.normal(size=m), rng.normal(size=m)
        q0, q1 = np.linalg.solve(S0, r0), np.linalg.solve(S1, r1)
        dS, dN = S1 - S0, N1 - N0
        chi = np.r_[r1 - r0, q1 - q0, dS @ q0, dN @ q0]
        lift = rank3.increment_product_ports(N1, S1)
        np.testing.assert_allclose(lift.equality @ chi, 0, atol=1e-12)
        np.testing.assert_allclose(
            lift.correction_difference @ chi,
            N1 @ q1 - N0 @ q0,
            atol=1e-12,
        )
        self.assertEqual(chi.size, 33)
        self.assertLess(chi.size, 105)


if __name__ == '__main__':
    unittest.main()
