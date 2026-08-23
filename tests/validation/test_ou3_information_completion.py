import importlib.util
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

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

    def test_record_index_accepts_short_slug_and_exact_source_stem(self):
        index = MOD.record_name_index()
        source = "wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv"
        self.assertEqual(index["jonswap_0_27"], source)
        self.assertEqual(index[Path(source).stem], source)

    def test_funnel_uses_tightest_certified_horizon_not_first_pass(self):
        attempts = [
            {"horizon_s": 2.0, "information_pass": True, "lambda_worst_information": 0.9999},
            {"horizon_s": 4.0, "information_pass": True, "lambda_worst_information": 0.99},
            {"horizon_s": 8.0, "information_pass": False, "lambda_worst_information": 1.01},
        ]
        def fake_eval(_record_data, mode, horizon, lam):
            b = {2.0: 1000.0, 4.0: 100.0}[horizon]
            return {"mode": mode, "status": "PASS", "horizon_s": horizon,
                    "lambda_information_bound": lam, "invariant_level_b_replay": b}
        with patch.object(MOD, "evaluate_mode", side_effect=fake_eval):
            selected, rows = MOD.evaluate_contracting_horizons({}, "H", attempts)
        self.assertEqual(selected["horizon_s"], 4.0)
        self.assertEqual(len(rows), 2)
        self.assertEqual(selected["selection_basis"],
                         "MINIMUM_REPLAY_INVARIANT_LEVEL_B_OVER_CERTIFIED_HORIZONS")

    def test_nominal_replay_cannot_promote_neighborhood_theorem(self):
        text = (TOOLS / "ou3_information_completion.py").read_text()
        self.assertIn("EXECUTED_NOISY_REPLAY_ONLY", text)
        self.assertIn('"numerical_neighborhood_certificate": "NOT_ESTABLISHED"', text)
        self.assertIn('"deployment_theorem_certificate": "NOT_ESTABLISHED"', text)


if __name__ == "__main__":
    unittest.main()
