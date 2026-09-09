"""Regression tests for the joint bias/source compatible-storage diagnostic."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "kalman_ou_iii"
SPEC = importlib.util.spec_from_file_location(
    "joint_bias_source", ROOT / "ou3_p4_joint_bias_source_storage.py")
J = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(J)


def sample(a, g=None, d=None, index=0):
    full = np.eye(21)
    full[:18, :18] = a
    if g is not None:
        full[:18, 18:] = g
    forcing = np.zeros(21)
    if d is not None:
        forcing[:18] = d
    step = {"stage": "projection", "kind": "accelerometer", "index": index}
    return {"index": index, "full": full, "forcing": forcing,
            "prefixes": [(full.copy(), forcing.copy(), step)]}


class CompatibleMetricTests(unittest.TestCase):
    def test_periodic_metric_closes_and_uses_one_cycle(self):
        a0 = .8*np.eye(18)
        a1 = .7*np.eye(18)
        metrics, q, report = J.periodic_metrics([sample(a0), sample(a1, index=1)],
                                                np.ones(18))
        self.assertEqual(len(metrics), 3)
        np.testing.assert_allclose(metrics[0], metrics[-1], rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(metrics[0], a0.T@metrics[1]@a0+q, rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(metrics[1], a1.T@metrics[2]@a1+q, rtol=2e-12, atol=2e-12)
        self.assertAlmostEqual(report["word_spectral_radius"], .56)

    def test_individually_unstable_sample_is_allowed_if_cycle_is_schur(self):
        a0 = np.eye(18)
        a1 = np.eye(18)
        a0[0, 0] = 1.2
        a1[0, 0] = .5
        metrics, _, report = J.periodic_metrics([sample(a0), sample(a1, index=1)],
                                                np.ones(18))
        self.assertGreater(a0[0, 0], 1)
        self.assertLess(report["word_spectral_radius"], 1)
        self.assertGreater(np.linalg.eigvalsh(metrics[0])[0], 0)

    def test_unstable_cycle_is_rejected(self):
        with self.assertRaises(ValueError):
            J.periodic_metrics([sample(1.01*np.eye(18))], np.ones(18))

    def test_one_step_supply_bounds_correlated_input(self):
        a = .8*np.eye(18)
        g = np.zeros((18, 3)); g[0, 0] = .2
        d = np.zeros(18); d[1] = .1
        mnext = 2*np.eye(18)
        q = np.eye(18)
        bias = np.array([.3, 0, 0])
        rate, supply, v, _ = J.one_step_budget(a, g, d, bias, mnext, q)
        mi = a.T@mnext@a+q
        rng = np.random.default_rng(7)
        for _ in range(100):
            x = rng.normal(size=18)
            lhs = (a@x+v) @ mnext @ (a@x+v)
            rhs = rate*(x@mi@x)+supply
            self.assertLessEqual(lhs, rhs+1e-10*max(1., abs(rhs)))

    def test_prefix_bound_keeps_bias_and_source_in_one_vector(self):
        a = .5*np.eye(18)
        g = np.zeros((18, 3)); g[0, 0] = 2
        d = np.zeros(18); d[0] = -1
        s = sample(a, g, d)
        result = J.prefix_bound(s, np.eye(18), 1., np.array([.5, 0, 0]), np.ones(18))
        # G*b+d cancels exactly.  Splitting their norms would not.
        self.assertAlmostEqual(result[0]["groups"]["attitude"]["joint_bias_source_excursion"], 0.)
        self.assertAlmostEqual(result[0]["groups"]["attitude"]["bound"], .5)


class BiasRecurrenceTests(unittest.TestCase):
    def rows(self, value=(0., 0., 0.)):
        return [{"word": "H18", "stage": "prediction_enter", "true_bias": list(value)}
                for _ in range(600)]

    def test_literal_zero_driver_is_reconstructed_not_promoted(self):
        r = J.true_bias_recurrence({"bias_root": [0, 0, 0], "bias_driver": "ZERO"},
                                   self.rows(), "H18")
        self.assertTrue(r["recurrence_reconstruction_pass"])
        self.assertTrue(r["zero_driver_is_not_uniform_BIAS1_coverage"])
        self.assertFalse(r["nonzero_driver_coverage"])

    def test_zero_driver_nonzero_root_must_continue_without_reseed(self):
        r = J.true_bias_recurrence({"bias_root": [.1, -.2, .3], "bias_driver": "ZERO"},
                                   self.rows((.1, -.2, .3)), "H18")
        self.assertTrue(r["recurrence_reconstruction_pass"])
        self.assertEqual(r["recurrence_reconstruction_max_abs"], 0.)

    def test_unexplained_nonzero_driver_is_not_inferred_from_trace(self):
        r = J.true_bias_recurrence({"bias_root": [0, 0, 0], "bias_driver": "SOMETHING"},
                                   self.rows(), "H18")
        self.assertFalse(r["recurrence_reconstruction_pass"])
        self.assertFalse(r["nonzero_driver_coverage"])

    def test_incomplete_history_fails(self):
        with self.assertRaises(ValueError):
            J.true_bias_recurrence({"bias_root": [0, 0, 0], "bias_driver": "ZERO"},
                                   self.rows()[:-1], "H18")


if __name__ == "__main__":
    unittest.main()
