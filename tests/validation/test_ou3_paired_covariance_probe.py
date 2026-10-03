import unittest

import numpy as np

from tools.stability.ou3_theorem.paired_covariance_probe import (
    certificate, correction, prior_to_posterior_drop, storage,
)


class PairedCovarianceProbeTests(unittest.TestCase):
    def test_exact_correction_storage_identity(self):
        c = certificate()
        self.assertLess(c["identity_defect"], 1e-10)
        self.assertGreater(c["measurement_loss_norm"], .4)
        self.assertFalse(c["shipping_service_proved"])

    def test_scalar_informative_correction_is_not_storage_invariant(self):
        P = X = H = R = np.ones((1, 1))
        Pp, Xp = correction(P, X, H, R)
        np.testing.assert_allclose(Pp, [[.5]])
        np.testing.assert_allclose(Xp, [[.5]])
        np.testing.assert_allclose(storage(Pp, Xp), [[.5]])
        result = prior_to_posterior_drop(P, X, H, R)
        np.testing.assert_allclose(result["storage_drop"], [[.5]])
        np.testing.assert_allclose(result["service"], [[.5]])
        self.assertLess(result["storage_defect"], 1e-14)

    def test_complete_cross_covariance_and_multiple_probe_columns(self):
        rng = np.random.default_rng(21)
        for n, m, columns in ((3, 1, 2), (6, 3, 4), (21, 3, 2)):
            with self.subTest(states=n):
                B = rng.normal(size=(n, n))
                P = B @ B.T + np.eye(n)
                B = rng.normal(size=(m, m))
                R = B @ B.T + np.eye(m)
                H = rng.normal(size=(m, n))
                X = rng.normal(size=(n, columns))
                result = prior_to_posterior_drop(P, X, H, R)
                np.testing.assert_allclose(result["storage_drop"], result["service"], atol=1e-12)
                self.assertGreaterEqual(np.linalg.eigvalsh(result["storage_drop"])[0], -1e-12)

    def test_zero_measurement_leaves_storage_unchanged(self):
        P = np.array([[2., .3], [.3, 1.]])
        result = prior_to_posterior_drop(P, np.eye(2), np.zeros((1, 2)), np.eye(1))
        np.testing.assert_allclose(result["storage_drop"], np.zeros((2, 2)), atol=1e-14)
        np.testing.assert_allclose(result["service"], np.zeros((2, 2)))


if __name__ == "__main__":
    unittest.main()
