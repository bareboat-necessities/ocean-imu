"""Semantic contract for the OU--III Phase-2 Live-basin certificate."""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "doc" / "kalman_ou_iii"
TEST = ROOT / "tests" / "kalman_ou_iii"


class OUIIILiveBasinCertificateContractTests(unittest.TestCase):
    def test_certificate_uses_riccati_metric_not_fixed_euclidean_contraction(self):
        text = (DOC / "w3d-live-basin-certificate.tex-part").read_text(encoding="utf-8")
        for marker in (
            r"\label{eq:iss-cert-Pz}",
            r"\label{eq:iss-cert-Vnorm}",
            r"\chi_{V,H}",
            r"\label{thm:iss-cert-finite-horizon}",
            r"\sqrt{\frac{\overline p}{\underline p}}",
            "0.9354",
            "0.9443",
            "$6.03$",
            "finite-grid proof",
        ):
            self.assertIn(marker, text)

    def test_diagnostic_uses_fixed_size_safe_svd_and_metric_gate(self):
        src = (TEST / "live_basin_diagnostic.cpp").read_text(encoding="utf-8")
        self.assertIn("covariance_metric_norm", src)
        self.assertIn("all_reference_metric_contract", src)
        self.assertNotIn("ComputeThinU", src)
        self.assertNotIn("ComputeThinV", src)
        run = (TEST / "run_tests.sh").read_text(encoding="utf-8")
        self.assertIn("./live_basin_diagnostic", run)

    def test_reference_grid_is_evidence_not_the_uniform_theorem(self):
        text = (DOC / "w3d-live-basin-certificate.tex-part").read_text(encoding="utf-8")
        self.assertIn("not a grid theorem", text)
        self.assertIn("An unexamined Live trajectory", text)


if __name__ == "__main__":
    unittest.main()
