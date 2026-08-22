"""Semantic contract for OU--III finite-step regional stability."""

import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
SRC = REPO_ROOT / "src"


def _read(name: str) -> str:
    return (DOC / name).read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return " ".join(text.split())


def _source(relative: str) -> str:
    return (SRC / relative).read_text(encoding="utf-8")


class OUIIIFiniteCaptureContractTests(unittest.TestCase):
    def test_main_article_wires_current_stability_chain(self):
        main = _read("kalman_ou-w3d.tex")
        local = r"\input{w3d-iss-stability.tex-part}"
        regional = r"\input{w3d-semiglobal-stability.tex-part}"
        evidence = r"\input{w3d-sim-charts.tex-part}"
        self.assertIn(local, main)
        self.assertIn(regional, main)
        self.assertLess(main.index(local), main.index(regional))
        self.assertLess(main.index(regional), main.index(evidence))

    def test_regional_capture_precedes_mahony_composition(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        self.assertTrue(proof.startswith(r"\input{w3d-block-local-iss.tex-part}"))
        capture = r"\input{w3d-finite-live-capture.tex-part}"
        self.assertIn(capture, proof)
        self.assertLess(capture and proof.index(capture), proof.index(r"\section{Almost-Global Proxy Entry"))

    def test_capture_uses_complete_nonlinear_map_and_exact_envelopes(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        for marker in (
            r"\label{eq:capture-exact-lift}",
            r"\mathcal F_{k,m}",
            r"\label{eq:capture-comparison-envelope}",
            r"\mathcal T_{k,m}",
            r"\label{eq:capture-comparison-monotone}",
            r"\label{eq:capture-comparison-sequence}",
            r"\label{thm:finite-live-capture}",
        ):
            self.assertIn(marker, capture)
        self.assertIn("full $S=0$ Kalman correction", capture)
        self.assertIn("No Schmidt, block-only, or gain-truncated", capture)

    def test_capture_allows_nonuniform_finite_step_transients(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        for marker in (
            r"\label{eq:capture-lift-domain}",
            r"\label{eq:capture-prefix-domain}",
            r"\label{eq:capture-finite-entry}",
            r"\label{eq:capture-inner-invariance}",
            r"\label{eq:capture-inner-contraction}",
            r"m_n\ge1",
        ):
            self.assertIn(marker, capture)
        self.assertIn("early post-handoff lifts may expand", _flat(capture))

    def test_capture_uses_structured_sensitive_kinematic_regions(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        for marker in (
            r"\label{eq:capture-structured-radii}",
            r"r_{\xi,k}",
            r"r_{\ell,k}",
            r"\label{eq:capture-structured-box}",
            r"0\prec\vct R_B\prec\vct R_C\preceq\vct R_D",
            r"\|\vct a\|_{\vct w,\infty}",
        ):
            self.assertIn(marker, capture)

    def test_capture_time_starts_from_handoff_envelope(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        for marker in (
            r"\label{eq:capture-handoff-sequence}",
            r"N_H:=\min",
            r"\label{eq:capture-count}",
            r"\label{eq:capture-time}",
        ):
            self.assertIn(marker, capture)

    def test_operational_contract_uses_structured_handoff_margin(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        startup = _read("w3d-init.tex-part")
        for marker in (
            r"\label{eq:capture-handoff-margin}",
            r"\label{eq:capture-certificate-outputs}",
            r"\label{eq:capture-certificate-pass}",
            r"\vct R_H\prec\vct R_C\preceq\vct R_D",
        ):
            self.assertIn(marker, capture)
        self.assertIn(r"\vct R_H^{\max}\prec\vct R_C", startup)

    def test_proxy_targets_structured_capture_region(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        startup = _read("w3d-init.tex-part")
        for marker in (
            r"\label{eq:semiglobal-handoff-metric-radius}",
            r"\label{eq:semiglobal-capture-entry}",
            r"\vct R_H^{\max}\prec\vct R_C",
            r"\ref{thm:finite-live-capture}",
            r"\label{thm:semiglobal-proxy-live}",
        ):
            self.assertIn(marker, proof)
        self.assertIn(r"\vct R_H^{\max}\prec\vct R_C", startup)

    def test_schedule_is_part_of_exact_regional_family(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\label{eq:semiglobal-mse-schedule}",
            r"\widehat\sigma_{a,B}^{\,6/7}",
            r"\tau^{41/14}",
            r"\label{eq:semiglobal-mse-powers}",
            r"\Pi_{k,m}",
        ):
            self.assertIn(marker, proof)

    def test_proxy_almost_global_scope_is_explicit(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\mathcal K_{\epsilon,B_\beta}",
            r"\mathcal N_P:=W^s(\mathcal U_P)",
            "antipodal",
            r"\cite{Mahony2008NonlinearComplementary}",
            r"\label{thm:semiglobal-proxy-live}",
        ):
            self.assertIn(marker, proof)

    def test_block_local_iss_is_retained_as_analytical_result(self):
        block = _read("w3d-block-local-iss.tex-part")
        for marker in (
            r"\label{lem:iss-block-remainder}",
            r"\label{thm:iss-block-local}",
            r"\label{eq:iss-block-basin}",
            r"\label{rem:iss-block-scope}",
            r"M_{\xi\ell}",
        ):
            self.assertIn(marker, block)

    def test_gaussian_noise_has_separate_stochastic_theorem(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        proof = _read("w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\label{thm:capture-stochastic-live}",
            r"\label{eq:capture-stochastic-drift}",
            r"\label{eq:capture-stochastic-ms}",
            r"\label{eq:capture-stochastic-finite-hp}",
        ):
            self.assertIn(marker, capture)
        self.assertIn(r"\ref{thm:capture-stochastic-live}", proof)

    def test_handoff_and_timeout_keep_aligned_branch_gate(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        startup = _read("w3d-init.tex-part")
        self.assertIn("0.075", proof)
        self.assertIn(r"\SI{150}{s}", proof)
        self.assertIn(r"\label{eq:semiglobal-aligned-branch}", proof)
        self.assertIn(r"\eqref{eq:semiglobal-aligned-branch}", startup)

        common = _source("kalman_common/SeaStateFusionFilterCommon.h")
        self.assertIn("inline bool gravityAlignedBranch", common)
        self.assertIn("unitVecAlignCos(s_obs, g_slow) > 0.0f", common)
        self.assertIn("if (!aligned_branch) return false;", common)

        for header in (
            "kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
            "kalman_ou_ii/SeaStateFusionFilter_OU_II.h",
        ):
            gate = _flat(_source(header))
            self.assertIn("seastate::common::gravityAlignedBranch", gate)
            self.assertIn("const bool tilt_trusted = mag_gravity_aligned_branch_ &&", gate)
            self.assertIn(
                "const bool ready_by_timeout = proxy_ready && (t_ >= timeout_sec) "
                "&& mag_gravity_aligned_branch_;",
                gate,
            )


if __name__ == "__main__":
    unittest.main()
