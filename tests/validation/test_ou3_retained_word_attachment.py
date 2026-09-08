from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/kalman_ou_iii"))
import ou3_p4_retained_word_attachment as A  # noqa: E402


class RetainedWordAttachmentTests(unittest.TestCase):
    def covariance(self, n):
        v = np.linspace(-0.5, 0.8, n)
        return np.diag(np.linspace(0.3, 1.2, n)) + np.outer(v, v)

    def event(self, n):
        H = np.zeros((3, n))
        H[:, :3] = -A.FAST._skew([0.1, 0.2, -9.8])
        H[:, 15:18] = np.eye(3)
        if n == 21:
            H[:, 18:21] = np.eye(3)
        P = self.covariance(n)
        R = np.diag([0.1, 0.2, 0.3])
        Pj, K, _ = A.BASE._joseph(P, H, R)
        d = np.array([0.3, -0.2, 0.1])
        G = A.reset_matrix(d, n)
        return {"type": A.BASE.EV_ACC, "H": H, "R": R, "dtheta": d,
                "Pbefore_shipping": P, "Pafter_shipping": G @ Pj @ G.T}, K, Pj

    def test_dense_joseph_covariance_and_gain_require_transformed_H(self):
        for n in (18, 21):
            ev, K, _ = self.event(n)
            T = np.linalg.inv(A.reset_matrix([0.4, 0.1, -0.2], n))
            d = A.transport_cell(ev, T, n)
            Pnew, Knew, _ = A.BASE._joseph(d["Pbefore"], d["H"], d["R"])
            np.testing.assert_allclose(Pnew, d["Pafter"], atol=1e-14)
            np.testing.assert_allclose(Knew, T @ K, atol=1e-14)
            wrong, _, _ = A.BASE._joseph(d["Pbefore"], ev["H"], ev["R"])
            self.assertGreater(np.linalg.norm(wrong - d["Pafter"]), 1e-3)
            np.testing.assert_array_equal(d["R"], ev["R"])

    def test_prediction_transports_F_and_Q_together(self):
        n = 21
        P = self.covariance(n)
        F = np.eye(n)
        F[:3, 3:6] = 0.005 * np.eye(3)
        Q = 0.01 * self.covariance(n)
        T = np.linalg.inv(A.reset_matrix([0.4, -0.1, 0.2], n))
        ev = {"type": A.BASE.EV_PRED, "Pbefore_shipping": P,
              "Pafter_shipping": F @ P @ F.T + Q, "linear_shipping": F, "Q": Q}
        d = A.transport_cell(ev, T, n)
        np.testing.assert_allclose(d["F"] @ d["Pbefore"] @ d["F"].T + d["Q"], d["Pafter"], atol=1e-14)
        self.assertGreater(np.linalg.norm(d["Q"] - Q), 1e-4)

    def test_finite_pullback_preserves_energy_for_nonlinear_map(self):
        for n in (18, 21):
            ev, K, Pj = self.event(n)
            x = np.linspace(-0.03, 0.05, n)
            x[:3] = [0.4, -0.2, 0.1]
            T0 = np.linalg.inv(A.reset_matrix([0.1, 0.2, 0.3], n))
            T1 = T0 @ np.linalg.inv(A.reset_matrix(ev["dtheta"], n))
            def finite(e, event=ev, dimension=n, gain=K):
                return A.FAST._joseph_value(event, e, dimension, gain, 0.4)[0]
            out = finite(x)
            z = A.finite_map_pullback(T0 @ x, T0, T1, finite)
            J = np.linalg.inv(Pj)
            Jz = np.linalg.inv(T1 @ Pj @ T1.T)
            self.assertAlmostEqual(float(z @ Jz @ z), float(out @ J @ out), places=11)
            self.assertGreater(np.linalg.norm(z - finite(T0 @ x)), 1e-3)

    def test_finite_reset_signed_identity_includes_projection_defect(self):
        n = 21
        ev, K, Pj = self.event(n)
        x = np.zeros(n)
        x[:3] = [0.5, 0.1, -0.2]
        x[18:21] = [0.7, 0.1, 0.2]
        # Keep this algebra fixture's projection active after the correction.
        K[18:21] = 0
        y = A.FAST._physical_residual(ev, x, n)
        out, branch = A.FAST._joseph_value(ev, x, n, K, 0.4)
        self.assertEqual(branch, "active")
        d = A.finite_reset_energy(x, out, K, y, Pj)
        self.assertGreater(d["reset_defect_norm"], 0.1)
        self.assertLess(d["relative_signed_identity_residual"], 1e-13)

    def test_capture_identity_and_caps_do_not_supply_missing_realization(self):
        ev, _, _ = self.event(21)
        d = A.source_evidence_inventory({"events": [ev], "source_id": "BRMM", "caps_pass": True})
        self.assertEqual(d["membership_decision"], "UNDETERMINED")
        self.assertIn("joint_response_witness", d["missing_common_realization_inputs"])
        self.assertIn("true_bias_root_and_history", d["missing_common_realization_inputs"])


if __name__ == "__main__":
    unittest.main()
