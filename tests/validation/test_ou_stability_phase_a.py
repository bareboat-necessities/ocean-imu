"""Semantic contract for OU--III stability hardening Phase A."""

import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"


def _read(name: str) -> str:
    return (DOC / name).read_text(encoding="utf-8")


class OUIIIStabilityPhaseAContractTests(unittest.TestCase):
    def test_heading_error_is_relative_to_adopted_gauge(self):
        proof = _read("w3d-iss-stability.tex-part")
        startup = _read("w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\label{eq:iss-gauge-reference}",
            r"\label{eq:iss-tilt-yaw-split}",
            r"\label{eq:iss-true-yaw-error}",
            r"b_{\psi,\mathrm{ref}}",
        ):
            self.assertIn(marker, proof)
        self.assertIn(r"\label{eq:semiglobal-handoff-yaw-internal}", startup)
        self.assertIn("not part of the MEKF capture error", startup)

    def test_accelerometer_bias_has_full_pe_observability_route(self):
        proof = _read("w3d-iss-stability.tex-part")
        for marker in (
            r"\label{eq:iss-eta-state}",
            r"\label{eq:iss-eta-gramian}",
            r"\label{eq:iss-eta-pe}",
            r"\label{lem:iss-21-uco-fullpe}",
            r"\delta b_{a,k+1}=\delta b_{a,k}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("random-walk limit", proof)
        self.assertIn("finite $\\tau_b$ remains a valid fallback", proof)

    def test_vector_information_is_windowed_and_asynchronous(self):
        proof = _read("w3d-iss-stability.tex-part")
        self.assertIn(r"\mathcal A_{k,N}", proof)
        self.assertIn(r"\mathcal M_{k,N}", proof)
        self.assertIn("They need not coincide", proof)
        self.assertIn("unequal sensor rates", proof)
        self.assertNotIn(r"\overline\omega\Delta_{\max}<\frac{\pi}{2}", proof)
        self.assertNotIn("two accepted accelerometer--magnetometer", proof)

    def test_s_observability_uses_any_four_usable_updates(self):
        proof = _read("w3d-iss-stability.tex-part")
        self.assertIn("four \\emph{usable} $S$", proof)
        self.assertIn(r"t_3\le T_{S,W}", proof)
        self.assertIn(r"\Delta_{S,-}", proof)
        self.assertIn(r"\frac{\Delta_{S,-}^6}{12}", proof)
        self.assertNotIn("four consecutive $S$", proof)

    def test_process_covariance_is_finite_horizon_not_pointwise_full_rank(self):
        proof = _read("w3d-iss-stability.tex-part")
        self.assertIn(r"0\preceq\mat Q_k^{\rm eff}\preceq q_+\mat I", proof)
        self.assertIn(r"\label{eq:iss-process-gramian}", proof)
        self.assertIn("Individual $\\mat Q_k^{\\rm eff}$ may therefore be singular", proof)
        self.assertNotIn(r"0<q_-\mat I\preceq\mat Q_k", proof)

    def test_no_magnetometer_mode_is_a_quotient_theorem(self):
        proof = _read("w3d-iss-stability.tex-part")
        startup = _read("w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\label{sec:iss-nomag-quotient}",
            r"\label{eq:iss-nomag-pe}",
            r"\label{thm:iss-nomag-quotient}",
            r"\label{eq:iss-nomag-ues}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("No convergence of absolute yaw is claimed", proof)
        self.assertIn(r"\ref{thm:iss-nomag-quotient}", startup)

    def test_full_heading_ues_keeps_fallback_and_relaxed_route(self):
        proof = _read("w3d-iss-stability.tex-part")
        self.assertIn(r"\label{thm:iss-21-ues}", proof)
        self.assertIn("full-PE route", proof)
        self.assertIn("OU-fallback route", proof)
        self.assertIn(r"\label{eq:iss-eta6-pe}", proof)
        self.assertNotIn("finite upper bound on $\\tau_b$ is essential", proof)


if __name__ == "__main__":
    unittest.main()
