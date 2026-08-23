import importlib.util
import math
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, TOOLS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


BASE = load("ou3_numerical_certificate")
COMP = load("ou3_certificate_completion")
CONTRACT = load("ou3_information_enclosure_contract")
ENC = load("ou3_validate_enclosure")


class Ou3CertificateCompletionTests(unittest.TestCase):
    def test_group_metric_matches_exact_so3_energy(self):
        # Legacy path-metric completion remains a diagnostic utility. P00=1 ->
        # a_R=2; zero xi; theta=pi/2 gives V_R=1 and W=2.
        P = np.eye(21)
        e = np.zeros(21)
        e[0] = math.pi / 2
        self.assertAlmostEqual(COMP.group_metric_value(e, P, "A"), 2.0, places=12)

    def test_recurrence_has_finite_capture_when_lambda_lt_one(self):
        r = COMP.recurrence_capture(c0=10.0, lam=0.8, gamma=0.1, b=0.5)
        self.assertEqual(r["status"], "PASS")
        self.assertGreater(r["N_H_words"], 0)

    def test_recurrence_refuses_noncontracting_word_family(self):
        r = COMP.recurrence_capture(c0=10.0, lam=1.01, gamma=0.1, b=1.0)
        self.assertEqual(r["status"], "NOT_AVAILABLE")

    @staticmethod
    def contract_mode():
        return {
            "mode": "A",
            "recommended_word_horizon_s": 4.0,
            "executed_reference_only": {
                "relative_Riccati_injection_margin_worst": 0.01,
                "Sigma_endpoint_lambda_min": 1.0e-6,
                "Sigma_endpoint_lambda_max": 100.0,
            },
        }

    @staticmethod
    def valid_mode_payload():
        return {
            "source_complete": True,
            "outward_rounded": True,
            "word_horizon_s": 4.0,
            "relative_Riccati_injection_margin_lower": 0.005,
            "Sigma_lambda_min_lower": 5.0e-7,
            "Sigma_lambda_max_upper": 120.0,
            "prefix_information_gain_upper": 2.0,
            "theta_star": 1.0,
            "mu_W_lower": 1.0e-4,
            "all_word_prefixes_safe": True,
        }

    def test_information_mode_accepts_strict_continuous_bounds(self):
        ans = ENC.validate_mode("A", self.valid_mode_payload(), self.contract_mode())
        self.assertTrue(ans["linear_pass"])
        self.assertTrue(ans["nonlinear_pass"])
        self.assertAlmostEqual(ans["lambda_information_upper"], 0.995)
        self.assertAlmostEqual(ans["information_metric_lambda_min_lower"], 1.0 / 120.0)
        self.assertAlmostEqual(ans["information_metric_lambda_max_upper"], 2.0e6)

    def test_information_mode_rejects_zero_injection(self):
        payload = self.valid_mode_payload()
        payload["relative_Riccati_injection_margin_lower"] = 0.0
        ans = ENC.validate_mode("A", payload, self.contract_mode())
        self.assertFalse(ans["linear_pass"])
        self.assertTrue(any("Riccati" in x for x in ans["failures"]))

    def test_information_mode_rejects_optimistic_anchor_exclusion(self):
        payload = self.valid_mode_payload()
        # The continuous source family contains the executed reference points,
        # so its minimum cannot exceed the observed minimum.
        payload["relative_Riccati_injection_margin_lower"] = 0.02
        ans = ENC.validate_mode("A", payload, self.contract_mode())
        self.assertFalse(ans["linear_pass"])
        self.assertTrue(any("exceeds included executed minimum" in x for x in ans["failures"]))

    def test_contract_prefers_strongest_executed_information_horizon(self):
        info = {
            "status": "PASS",
            "held": {
                "selected": {"horizon_s": 2.0, "lambda_worst_information": 0.9999},
                "strongest_executed_margin": {
                    "horizon_s": 16.0,
                    "lambda_worst_information": 0.95,
                    "relative_Riccati_injection_margin_worst": 0.05,
                    "Sigma_endpoint_lambda_min": 1e-7,
                    "Sigma_endpoint_lambda_max": 10.0,
                },
            },
            "active": {
                "selected": {"horizon_s": 2.0, "lambda_worst_information": 0.999},
                "strongest_executed_margin": {
                    "horizon_s": 4.0,
                    "lambda_worst_information": 0.99,
                    "relative_Riccati_injection_margin_worst": 0.01,
                    "Sigma_endpoint_lambda_min": 1e-7,
                    "Sigma_endpoint_lambda_max": 10.0,
                },
            },
        }
        completion = {
            "status": "PASS_EXECUTED_REPLAY",
            "held": {"asymptotic_floor_b_star_replay": 10.0},
            "active": {"asymptotic_floor_b_star_replay": 20.0},
        }
        c = CONTRACT.build_contract(info, completion)
        self.assertEqual(c["modes"]["H"]["recommended_word_horizon_s"], 16.0)
        self.assertEqual(c["modes"]["A"]["recommended_word_horizon_s"], 4.0)
        self.assertEqual(c["metric"], "source-varying inverse estimator covariance")

    def test_completion_does_not_promote_sampled_replay(self):
        text = (TOOLS / "ou3_information_completion.py").read_text()
        self.assertIn('"numerical_neighborhood_certificate": "NOT_ESTABLISHED"', text)
        self.assertIn('"deployment_theorem_certificate": "NOT_ESTABLISHED"', text)

    def test_enclosure_gate_uses_information_metric_not_old_path_metrics(self):
        text = (TOOLS / "ou3_validate_enclosure.py").read_text()
        self.assertIn("relative_Riccati_injection_margin_lower", text)
        self.assertIn("Sigma_lambda_min_lower", text)
        self.assertIn("prefix_information_gain_upper", text)
        self.assertIn("mu_W_lower", text)
        self.assertIn("theta_star", text)
        self.assertNotIn("path_metrics.npz", text)
        self.assertNotIn("robust_box_lmi_upper", text)


if __name__ == "__main__":
    unittest.main()
