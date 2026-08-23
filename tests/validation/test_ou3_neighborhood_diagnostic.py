import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))

SPEC = importlib.util.spec_from_file_location(
    "ou3_neighborhood_diagnostic", TOOLS / "ou3_neighborhood_diagnostic.py"
)
MOD = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


class Ou3NeighborhoodDiagnosticTests(unittest.TestCase):
    def test_information_normalization_hits_requested_energy(self):
        P = np.diag([4.0, 9.0, 16.0])
        d = MOD.information_normalize(P, 1, target_W=0.25, sign=-1)
        W = float(d @ np.linalg.solve(P, d))
        self.assertAlmostEqual(W, 0.25, places=13)
        self.assertLess(d[1], 0.0)
        self.assertEqual(np.count_nonzero(d), 1)

    def test_held_compact_basis_never_injects_accel_bias(self):
        self.assertTrue(all(i < 18 for i in MOD.COMPACT_H))
        self.assertIn(18, MOD.COMPACT_A)

    def test_pair_trace_requires_source_match_and_contraction(self):
        header = (
            "time_s,time_from_injection_s,endpoint,mode,source_match,"
            "acc_accept_match,mag_accept_match,W_nominal,W_perturbed,"
            "theta_rad,error_norm,covariance_rel_fro,nom_live,nom_active,"
            "nom_mag_lock,nom_mag_refined,nom_tau,nom_sigma,nom_rs,"
            "pert_tau,pert_sigma,pert_rs"
        )
        rows = [
            "300,0,0,A,1,1,1,1.0,1.0,0.01,0.01,0,1,1,1,1,2,1,4,2,1,4",
            "304,4,1,A,1,1,1,0.8,0.81,0.005,0.008,0.01,1,1,1,1,2,1,4,2,1,4",
        ]
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "trace.csv"
            p.write_text(header + "\n" + "\n".join(rows) + "\n")
            result = MOD.parse_trace(p)
        self.assertEqual(result["status"], "PASS_SAMPLED")
        self.assertTrue(result["pass_sampled"])
        self.assertAlmostEqual(result["relative_decrement"], 0.2)

    def test_pair_trace_rejects_source_word_divergence(self):
        header = (
            "time_s,time_from_injection_s,endpoint,mode,source_match,"
            "acc_accept_match,mag_accept_match,W_nominal,W_perturbed,"
            "theta_rad,error_norm,covariance_rel_fro,nom_live,nom_active,"
            "nom_mag_lock,nom_mag_refined,nom_tau,nom_sigma,nom_rs,"
            "pert_tau,pert_sigma,pert_rs"
        )
        rows = [
            "300,0,0,A,1,1,1,1.0,1.0,0.01,0.01,0,1,1,1,1,2,1,4,2,1,4",
            "304,4,1,A,0,0,1,0.7,0.7,0.005,0.008,0.01,1,1,1,1,2,1,4,2,1,4",
        ]
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "trace.csv"
            p.write_text(header + "\n" + "\n".join(rows) + "\n")
            result = MOD.parse_trace(p)
        self.assertEqual(result["status"], "FAIL_SAMPLED")
        self.assertFalse(result["source_match_all"])
        self.assertFalse(result["measurement_acceptance_match_all"])

    def test_sampled_diagnostic_cannot_promote_theorem(self):
        text = (TOOLS / "ou3_neighborhood_diagnostic.py").read_text()
        self.assertIn("SAMPLED_PAIRWISE_NONLINEAR_DIAGNOSTIC_ONLY", text)
        self.assertIn('"numerical_neighborhood_certificate": "NOT_ESTABLISHED"', text)
        self.assertIn('"deployment_theorem_certificate": "NOT_ESTABLISHED"', text)


if __name__ == "__main__":
    unittest.main()
