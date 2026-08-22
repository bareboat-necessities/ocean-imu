"""Semantic contract for OU--III finite capture and proxy-to-Live stability."""

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
    def test_main_article_wires_composite_stability_after_live_iss(self):
        main = _read("kalman_ou-w3d.tex")
        local = r"\input{w3d-iss-stability.tex-part}"
        semiglobal = r"\input{w3d-semiglobal-stability.tex-part}"
        evidence = r"\input{w3d-sim-charts.tex-part}"
        self.assertIn(semiglobal, main)
        self.assertLess(main.index(local), main.index(semiglobal))
        self.assertLess(main.index(semiglobal), main.index(evidence))
        self.assertNotIn(r"\input{w3d-computer-assisted-live-basin.tex-part}", main)

    def test_capture_is_inserted_before_mahony_composition(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        self.assertTrue(proof.startswith(r"\input{w3d-block-local-iss.tex-part}"))
        self.assertIn(r"\input{w3d-finite-live-capture.tex-part}", proof)
        self.assertLess(
            proof.index(r"\input{w3d-finite-live-capture.tex-part}"),
            proof.index(r"\section{Almost-Global Proxy Entry"),
        )

    def test_operational_capture_uses_complete_nonlinear_map(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        for marker in (
            r"\label{eq:capture-exact-lift}",
            r"\label{eq:capture-comparison-envelope}",
            r"\label{eq:capture-comparison-monotone}",
            r"\label{eq:capture-comparison-sequence}",
            r"\label{thm:finite-live-capture}",
            r"\label{eq:capture-practical-iss}",
        ):
            self.assertIn(marker, capture)
        flat = _flat(capture)
        self.assertIn("complete implemented nonlinear error map", flat)
        self.assertIn("No individual sample, and no individual lifted interval, is required to contract", flat)
        self.assertIn("preserves any cancellation produced by the Kalman correction", flat)
        self.assertIn("No Schmidt, block-only, or gain-truncated approximation", flat)

    def test_capture_allows_nonuniform_expanding_early_lifts(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        for marker in (
            r"\label{eq:capture-lift-domain}",
            r"\label{eq:capture-prefix-domain}",
            r"\label{eq:capture-finite-entry}",
            r"\label{eq:capture-inner-invariance}",
            r"\label{eq:capture-inner-contraction}",
        ):
            self.assertIn(marker, capture)
        flat = _flat(capture)
        self.assertIn("the first post-handoff lifts may expand", flat)
        self.assertIn("finite-step Lyapunov/comparison result rather than a one-step contraction theorem", flat)

    def test_operational_capture_is_not_bottlenecked_by_block_local_basin(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        proof = _read("w3d-semiglobal-stability.tex-part")
        block = _read("w3d-block-local-iss.tex-part")
        self.assertIn(
            "No inclusion in the separate block-local entrance condition is required",
            _flat(capture),
        )
        self.assertNotIn("is contained in the block-local ISS entrance set", _flat(capture))
        self.assertIn(
            "It is not an additional set-inclusion test in the operational startup chain",
            _flat(proof),
        )
        self.assertIn(
            "startup proof therefore does not need to show that its regional inner set is contained",
            _flat(block),
        )

    def test_capture_uses_structured_two_block_regions(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        for marker in (
            r"\label{eq:capture-structured-radii}",
            r"\label{eq:capture-structured-box}",
            r"0\prec\vct R_B\prec\vct R_C\preceq\vct R_D",
            r"\|\vct a\|_{\vct w,\infty}",
        ):
            self.assertIn(marker, capture)
        flat = _flat(capture)
        self.assertIn("does not force the admissible attitude/$S$/bias errors and the kinematic errors into one spherical radius", flat)

    def test_capture_time_iterates_from_actual_handoff(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        for marker in (
            r"\label{eq:capture-handoff-sequence}",
            r"N_H:=\min",
            r"\label{eq:capture-count}",
            r"\label{eq:capture-time}",
        ):
            self.assertIn(marker, capture)
        self.assertNotIn(r"\frac{R_C-R_B}{\mu_C}", capture)
        self.assertNotIn(r"\mu_C", capture)
        self.assertIn(
            "startup claim is not penalized by pretending every accepted Mahony handoff begins on the extreme boundary",
            _flat(capture),
        )

    def test_no_arbitrary_factor_five_is_a_stability_gate(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        proof = _read("w3d-semiglobal-stability.tex-part")
        startup = _read("w3d-init.tex-part")
        for marker in (
            r"\label{eq:capture-handoff-margin}",
            r"\label{eq:capture-certificate-pass}",
            r"\vct R_H\prec\vct R_C\preceq\vct R_D",
        ):
            self.assertIn(marker, capture)
        self.assertNotIn(r"\Gamma_{\rm req}=5", capture)
        self.assertNotIn(r"\Gamma_{\rm cap}\ge5", capture)
        self.assertNotIn(r"R_C/R_B\ge5", proof)
        self.assertIn("No prescribed ratio between outer and inner radii", startup)

    def test_proxy_targets_structured_capture_region_not_local_basin(self):
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
        self.assertNotIn(r"\label{eq:semiglobal-basin-entry}", proof)
        self.assertIn(
            "It does not assert direct entry into either the regional inner set or the separate block-local ISS neighborhood",
            _flat(proof),
        )
        self.assertIn(r"\vct R_H^{\max}\prec\vct R_C", startup)
        self.assertIn(r"\ref{thm:finite-live-capture}", startup)
        self.assertIn("not required to satisfy the block-local invariant-tube", _flat(startup))

    def test_extension_uses_spectral_mse_schedule_not_cubic_adaptation(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\label{eq:semiglobal-mse-schedule}",
            r"\widehat\sigma_{a,B}^{\,6/7}",
            r"\tau^{41/14}",
            r"\label{eq:semiglobal-mse-powers}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("No $\\sigma\\tau^3$ adaptation relation", proof)
        self.assertIn(r"\Pi_{k,m}", proof)

    def test_extension_is_almost_global_only_through_proxy(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\mathcal K_{\epsilon,B_\beta}",
            r"\mathcal N_P:=W^s(\mathcal U_P)",
            "antipodal",
            r"\cite{Mahony2008NonlinearComplementary}",
            r"\label{thm:semiglobal-proxy-live}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn("does not claim semiglobal convergence", flat)
        self.assertIn("attitude-startup qualification is almost global", flat)

    def test_block_local_analysis_is_retained_as_independent_result(self):
        block = _read("w3d-block-local-iss.tex-part")
        proof = _read("w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\label{lem:iss-block-remainder}",
            r"\label{thm:iss-block-local}",
            r"\label{eq:iss-block-basin}",
            r"\delta\vct S",
            r"\delta\vct v",
            r"\delta\vct p",
            r"M_{\xi\ell}",
        ):
            self.assertIn(marker, block)
        self.assertIn(r"B_{\xi,H}", proof)
        self.assertIn(r"B_{\ell,H}", proof)
        self.assertIn("independent analytical statement", _flat(proof))
        self.assertNotIn("Theorem~\\ref{thm:iss-block-local} applies after capture", _flat(proof))

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
        flat = _flat(capture)
        self.assertIn("Gaussian process and measurement noise are treated separately", flat)
        self.assertIn("mean-square practical", flat)
        self.assertIn(r"\ref{thm:capture-stochastic-live}", proof)

    def test_handoff_and_timeout_keep_aligned_branch_gate(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        startup = _read("w3d-init.tex-part")
        self.assertIn("0.075", proof)
        self.assertIn(r"\SI{150}{s}", proof)
        self.assertIn(r"\label{eq:semiglobal-aligned-branch}", proof)
        self.assertIn("sign of the inner product", proof)
        self.assertIn("withholds a timeout-forced handoff until", _flat(proof))
        self.assertIn(r"\eqref{eq:semiglobal-aligned-branch}", startup)
        self.assertIn(
            "delayed rather than accepted on the antipodal branch",
            _flat(startup),
        )

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
            self.assertIn(
                "const bool tilt_trusted = mag_gravity_aligned_branch_ &&", gate
            )
            self.assertIn(
                "const bool ready_by_timeout = proxy_ready && (t_ >= timeout_sec) "
                "&& mag_gravity_aligned_branch_;",
                gate,
            )


if __name__ == "__main__":
    unittest.main()
