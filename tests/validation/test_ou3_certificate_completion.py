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
ENC = load("ou3_validate_enclosure")


class Ou3CertificateCompletionTests(unittest.TestCase):
    def test_group_metric_matches_exact_so3_energy(self):
        # P00=1 -> a_R=2; zero xi; theta=pi/2 gives V_R=1 and W=2.
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

    def test_robust_matrix_box_lmi_accepts_strict_contraction(self):
        P = np.eye(3)
        ans = ENC.robust_box_lmi_upper(0.8*np.eye(3), np.zeros((3,3)), P, P)
        self.assertTrue(ans["pass"])
        self.assertLess(ans["robust_difference_lambda_upper"], 0.0)

    def test_robust_matrix_box_lmi_rejects_expansion(self):
        P = np.eye(3)
        ans = ENC.robust_box_lmi_upper(1.05*np.eye(3), np.zeros((3,3)), P, P)
        self.assertFalse(ans["pass"])
        self.assertGreater(ans["robust_difference_lambda_upper"], 0.0)

    def test_interval_radius_is_penalized(self):
        P = np.eye(2)
        exact = ENC.robust_box_lmi_upper(0.7*np.eye(2), np.zeros((2,2)), P, P)
        box = ENC.robust_box_lmi_upper(0.7*np.eye(2), 0.05*np.ones((2,2)), P, P)
        self.assertGreater(box["robust_difference_lambda_upper"],
                           exact["robust_difference_lambda_upper"])

    def test_completion_does_not_promote_sampled_replay(self):
        text = (TOOLS / "ou3_certificate_completion.py").read_text()
        self.assertIn('"neighborhood_numerical_certificate": "NOT_ESTABLISHED"', text)
        self.assertIn('"deployment_theorem_certificate": "NOT_ESTABLISHED"', text)
        self.assertIn("enclosure_contract.json", text)

    def test_enclosure_gate_requires_validated_provenance(self):
        text = (TOOLS / "ou3_validate_enclosure.py").read_text()
        self.assertIn("validated_arithmetic", text)
        self.assertIn("outward_rounding", text)
        self.assertIn("source_generated_not_trajectory_fit", text)
        self.assertIn("robust_box_lmi_upper", text)
        self.assertIn("mu_W_lower", text)
        self.assertIn("theta_star", text)


if __name__ == "__main__":
    unittest.main()
