import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))

SPEC = importlib.util.spec_from_file_location(
    "ou3_information_certificate", TOOLS / "ou3_information_certificate.py"
)
INFO = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = INFO
SPEC.loader.exec_module(INFO)


class Ou3InformationCertificateTests(unittest.TestCase):
    def test_information_metric_contracts_when_covariance_has_positive_increment(self):
        A = np.diag([1.02, 0.95, 0.8])
        P0 = np.diag([2.0, 1.0, 0.5])
        Q = np.diag([0.3, 0.2, 0.1])
        P1 = A @ P0 @ A.T + Q
        lam = INFO.information_lambda(A, P0, P1)
        self.assertLess(lam, 1.0)
        inc = INFO.covariance_increment_margin(A, P0, P1)
        self.assertGreater(inc["omega_relative_lambda_min"], 0.0)

    def test_information_metric_detects_inconsistent_covariance(self):
        A = 1.1 * np.eye(2)
        P0 = np.eye(2)
        P1 = np.eye(2)
        inc = INFO.covariance_increment_margin(A, P0, P1)
        self.assertLess(inc["omega_relative_lambda_min"], 0.0)
        self.assertGreater(INFO.information_lambda(A, P0, P1), 1.0)

    def test_tool_uses_estimator_covariance_not_truth_error_covariance(self):
        text = (TOOLS / "ou3_information_certificate.py").read_text()
        self.assertIn("Sigma_KF", text)
        self.assertIn("Phi Sigma0 Phi' + Omega", text)
        self.assertNotIn("metric_from_samples", text)
        self.assertNotIn("X.T@X", text)

    def test_host_observer_records_full_covariance_at_map_boundaries(self):
        text = (ROOT / "tests" / "kalman_ou_iii" / "ou3-information-certificate-sim.cpp").read_text()
        self.assertIn("covariance_full()", text)
        self.assertIn("OU3COV1", text)
        self.assertIn("CertificateAdapter inner_", text)
        self.assertIn("map_block_started_", text)
        self.assertNotIn("enableLinearBlock(false)", text)
        self.assertNotIn("setFixedTuning", text)


if __name__ == "__main__":
    unittest.main()
