import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))

SPEC = importlib.util.spec_from_file_location(
    "ou3_information_completion", TOOLS / "ou3_information_completion.py"
)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


class Ou3InformationCompletionTests(unittest.TestCase):
    def test_information_energy_uses_full_inverse_covariance(self):
        e = np.array([1.0, 2.0])
        P = np.array([[2.0, 0.4], [0.4, 1.0]])
        expected = float(e @ np.linalg.solve(P, e))
        self.assertAlmostEqual(MOD.info_energy(e, P), expected, places=13)

    def test_capture_recursion_closes_for_strict_affine_contraction(self):
        lam = 0.8
        gamma = 0.1
        b = gamma / (1.0 - lam)
        n = MOD.finite_capture_steps(10.0, lam, gamma, b * 1.01)
        self.assertIsNotNone(n)
        self.assertGreater(n, 0)

    def test_capture_refuses_noncontracting_lambda(self):
        self.assertIsNone(MOD.finite_capture_steps(1.0, 1.0, 0.1, 1.0))
        self.assertIsNone(MOD.finite_capture_steps(1.0, 1.01, 0.1, 1.0))

    def test_nominal_replay_cannot_promote_neighborhood_theorem(self):
        text = (TOOLS / "ou3_information_completion.py").read_text()
        self.assertIn("EXECUTED_NOISY_REPLAY_ONLY", text)
        self.assertIn('"numerical_neighborhood_certificate": "NOT_ESTABLISHED"', text)
        self.assertIn('"deployment_theorem_certificate": "NOT_ESTABLISHED"', text)


if __name__ == "__main__":
    unittest.main()
