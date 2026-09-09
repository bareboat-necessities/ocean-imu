"""Algebraic regressions for the nonzero physical-bias driver lift."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "kalman_ou_iii"
sys.path.insert(0, str(ROOT))
SPEC = importlib.util.spec_from_file_location(
    "joint_nonzero_driver", ROOT / "ou3_p4_joint_nonzero_bias_driver.py")
D = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(D)


class DeclaredDriverTests(unittest.TestCase):
    def setUp(self):
        self.root = {
            "bias_driver": "DETERMINISTIC_SINUSOIDAL_GM",
            "bias_root": [0.08, -0.05, 0.03],
            "bias_tau_true_s": 1200.0,
            "bias_driver_omega_rad_s": 2*np.pi/600,
            "bias_driver_amplitude": [0.015, 0.010, -0.008],
            "bias_driver_phase": [0.0, 1.1, -0.7],
        }

    def test_affine_increment_reconstructs_declared_bias(self):
        dt = .005
        for t in (.005, 1., 33., 1144.):
            phi, w = D.driver_increment(self.root, t, dt)
            before = D.declared_bias(self.root, t-dt)
            after = D.declared_bias(self.root, t)
            np.testing.assert_allclose(phi*before+w, after, rtol=0, atol=2e-17)
            self.assertGreater(np.linalg.norm(w), 0)

    def test_driver_requires_complete_declaration(self):
        bad = dict(self.root)
        bad["bias_driver"] = "ZERO"
        with self.assertRaises(ValueError):
            D.declared_bias(bad, 1.)
        bad = dict(self.root)
        del bad["bias_driver_phase"]
        with self.assertRaises(ValueError):
            D.declared_bias(bad, 1.)

    def test_prediction_injects_same_driver_into_error_and_truth(self):
        dt = .005
        t = 3.0
        phi_hat = .999
        a21 = np.eye(21)
        a21[18:21, 18:21] *= phi_hat
        step = {"A": a21, "B": np.zeros((21, 1)), "u": np.ones(1),
                "stage": "prediction", "kind": "none", "index": 7,
                "bias_cost": False}
        row = {"word": "A21", "stage": "prediction", "index": 7,
               "source_time": t}
        original = D.J.projection_scales
        D.J.projection_scales = lambda rows, mode: {}
        try:
            lifted, meta = D.augmented_steps(self.root, [row], [step], "A21", dt)
        finally:
            D.J.projection_scales = original
        a = lifted[0]["A24"]
        f = lifted[0]["forcing24"]
        phi_true, w = D.driver_increment(self.root, t, dt)
        np.testing.assert_allclose(a[18:21, 21:24], (phi_true-phi_hat)*np.eye(3))
        np.testing.assert_allclose(a[21:24, 21:24], phi_true*np.eye(3))
        np.testing.assert_allclose(f[18:21], w)
        np.testing.assert_allclose(f[21:24], w)
        self.assertTrue(meta["driver_affine_input_lift_present"])

    def test_projection_couples_true_bias_on_same_graph(self):
        a21 = np.eye(21)
        step = {"A": a21, "B": np.zeros((21, 1)), "u": np.ones(1),
                "stage": "projection", "kind": "accelerometer", "index": 8,
                "bias_cost": False}
        original = D.J.projection_scales
        D.J.projection_scales = lambda rows, mode: {(8, "accelerometer"): .75}
        try:
            lifted, meta = D.augmented_steps(self.root, [], [step], "A21", .005)
        finally:
            D.J.projection_scales = original
        a = lifted[0]["A24"]
        np.testing.assert_allclose(a[18:21, 21:24], .25*np.eye(3))
        np.testing.assert_allclose(a[21:24, 21:24], np.eye(3))
        self.assertEqual(meta["minimum_projection_scale"], .75)
        self.assertEqual(meta["saturated_projection_count"], 1)
        self.assertTrue(meta["projection_beta_coupling_present"])


if __name__ == "__main__":
    unittest.main()
