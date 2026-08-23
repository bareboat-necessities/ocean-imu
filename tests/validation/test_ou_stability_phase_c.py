"""Semantic/source contract for exact nonlinear Live and startup composition."""

import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
SRC = REPO_ROOT / "src"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


class OUIIIStabilityPhaseCContractTests(unittest.TestCase):
    def test_obsolete_local_iss_certificate_is_removed(self):
        self.assertFalse((DOC / "w3d-block-local-iss.tex-part").exists())
        main = _read(DOC / "kalman_ou-w3d.tex")
        semiglobal = _read(DOC / "w3d-semiglobal-stability.tex-part")
        self.assertNotIn("w3d-block-local-iss", main + semiglobal)

    def test_exact_large_angle_group_correction_is_primary(self):
        proof = _read(DOC / "w3d-stability-widening-phase-c.tex-part")
        for marker in (
            r"\label{eq:widen-SO3-energy}",
            r"\label{eq:widen-exact-source-correction}",
            r"\label{eq:widen-exact-energy-change}",
            r"\label{eq:widen-large-angle-sector}",
            r"\label{thm:widen-large-angle-dissipation}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("full $S\\to\\theta$", proof)
        self.assertIn(r"\theta_e\le\theta_{*,i}<\pi", proof)

    def test_explicit_mahony_proxy_bound_is_retained(self):
        proof = _read(DOC / "w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\label{eq:semiglobal-proxy-gains}",
            r"\label{eq:semiglobal-proxy-force-bound}",
            r"\label{eq:semiglobal-proxy-reduced-dynamics}",
            r"\label{eq:semiglobal-proxy-lyapunov}",
            r"\label{eq:semiglobal-proxy-iss-differential}",
            r"\label{thm:semiglobal-proxy-explicit}",
            r"\delta_P=\SI{5}{s}",
            r"q_P=0.5",
        ):
            self.assertIn(marker, proof)
        self.assertIn("almost-global", proof)
        self.assertIn("antipodal", proof)

    def test_quality_gate_initializes_funnel_not_outer_box(self):
        proof = _read(DOC / "w3d-semiglobal-stability.tex-part")
        init = _read(DOC / "w3d-init.tex-part")
        for marker in (
            r"\label{eq:semiglobal-quality-entry-time}",
            r"\label{eq:semiglobal-aligned-branch}",
            r"\label{eq:semiglobal-capture-entry}",
            r"\mathcal H_{k_H}\subset\mathcal D_{k_H}",
            r"\label{thm:semiglobal-proxy-live}",
            r"\eqref{eq:capture-initial-level}",
        ):
            self.assertIn(marker, proof)
        self.assertIn(r"\mathcal H_{k_H}\subset\mathcal D_{k_H}", init)
        self.assertNotIn(r"R_C", proof)
        self.assertIn("radius summary is not propagated", proof)

    def test_timeout_has_gauged_and_quotient_routes(self):
        proof = _read(DOC / "w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\label{eq:semiglobal-timeout-150}",
            r"\label{eq:semiglobal-timeout-energy}",
            r"\label{eq:semiglobal-timeout-tilt}",
            r"\label{eq:semiglobal-timeout-capture-entry}",
            r"\label{eq:semiglobal-timeout-quotient-entry}",
            r"\ref{thm:hybrid-quotient-capture}",
            r"\label{thm:semiglobal-timeout-live}",
        ):
            self.assertIn(marker, proof)

    def test_proxy_source_gains_and_timeout_still_match_source(self):
        wrapper = _read(SRC / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h")
        self.assertIn("constexpr float STARTUP_PROXY_TWO_KP_DEFAULT = 0.2f;", wrapper)
        self.assertIn("constexpr float STARTUP_PROXY_TWO_KI_DEFAULT = 0.02f;", wrapper)
        self.assertIn("float proxy_startup_timeout_sec = 150.0f;", wrapper)
        self.assertIn("/*allow_acc_bias=*/false", wrapper)


if __name__ == "__main__":
    unittest.main()
