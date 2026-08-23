import importlib.util
import math
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))

SPEC = importlib.util.spec_from_file_location(
    "ou3_numerical_certificate", TOOLS / "ou3_numerical_certificate.py"
)
CERT = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = CERT
SPEC.loader.exec_module(CERT)


class Ou3NumericalCertificateTests(unittest.TestCase):
    def test_reference_inventory_is_exactly_eight_noisy_publication_seas(self):
        self.assertEqual(len(CERT.RECORDS), 8)
        self.assertEqual(
            [(family, hs) for family, hs, _ in CERT.RECORDS],
            [
                ("JONSWAP", 0.27), ("JONSWAP", 1.50),
                ("JONSWAP", 4.00), ("JONSWAP", 8.50),
                ("PM-Stokes", 0.27), ("PM-Stokes", 1.50),
                ("PM-Stokes", 4.00), ("PM-Stokes", 8.50),
            ],
        )

    def test_exact_so3_group_energy_is_used(self):
        self.assertAlmostEqual(CERT.group_energy(np.zeros(3)), 0.0, places=14)
        self.assertAlmostEqual(
            CERT.group_energy(np.array([math.pi / 2, 0.0, 0.0])), 1.0, places=12
        )
        self.assertAlmostEqual(
            CERT.group_energy(np.array([math.pi, 0.0, 0.0])), 2.0, places=12
        )

    def test_generalized_word_factor_distinguishes_contraction(self):
        P = np.eye(3)
        self.assertLess(CERT.generalized_lambda(0.8 * np.eye(3), P, P), 1.0)
        self.assertGreater(CERT.generalized_lambda(1.05 * np.eye(3), P, P), 1.0)

    def test_zu_ned_basis_matches_filter_contract(self):
        got = CERT.zu_to_ned(np.array([1.0, 2.0, 3.0]))
        np.testing.assert_allclose(got, np.array([2.0, 1.0, -3.0]))

    def test_replay_result_cannot_be_promoted_to_theorem_certificate(self):
        text = (TOOLS / "ou3_numerical_certificate.py").read_text()
        self.assertIn('"deployment_theorem_certificate": "NOT_ESTABLISHED"', text)
        self.assertIn("validated continuous-source word-family enclosure", text)
        self.assertNotIn("Markov", text)

    def test_certificate_simulator_observes_full_s_cross_covariance(self):
        text = (ROOT / "tests" / "kalman_ou_iii" / "ou3-certificate-sim.cpp").read_text()
        self.assertIn("covariance_full()", text)
        self.assertIn("P(:,S)", text)
        self.assertIn("get_integral_displacement()", text)
        self.assertIn("setPeriodicAwCovarianceSync(true)", text)
        self.assertNotIn("enableLinearBlock(false)", text)
        self.assertNotIn("setFixedTuning", text)


if __name__ == "__main__":
    unittest.main()
