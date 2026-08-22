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

    def test_capture_is_inserted_before_mahony_composition(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        self.assertTrue(proof.startswith(r"\input{w3d-block-local-iss.tex-part}"))
        self.assertIn(r"\input{w3d-finite-live-capture.tex-part}", proof)
        self.assertLess(
            proof.index(r"\input{w3d-finite-live-capture.tex-part}"),
            proof.index(r"\section{Almost-Global Proxy Entry"),
        )

    def test_finite_capture_theorem_allows_oscillatory_prefix(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        for marker in (
            r"\label{thm:finite-live-capture}",
            r"\label{eq:capture-linear-contraction}",
            r"\label{eq:capture-nonlinear-bound}",
            r"\label{eq:capture-disturbance-bound}",
            r"\label{eq:capture-radius-order}",
            r"\label{eq:capture-count}",
            r"\label{eq:capture-practical-iss}",
        ):
            self.assertIn(marker, capture)
        flat = _flat(capture)
        self.assertIn("Individual samples, and even short prefixes, need not contract", flat)
        self.assertIn("A finite prefix amplification larger than one is permitted", flat)
        self.assertIn("full Kalman correction", flat)
        self.assertIn("no Schmidt, block-only, or gain-truncated approximation", flat)

    def test_operational_capture_region_is_materially_wider_than_iss_entry(self):
        capture = _read("w3d-finite-live-capture.tex-part")
        for marker in (
            r"\Gamma_{\rm cap}:=\frac{R_C}{R_B}",
            r"\Gamma_{\rm req}=5",
            r"\Gamma_{\rm cap}\ge5",
            r"R_H^{\max}<R_C",
            r"\label{eq:capture-certificate-pass}",
        ):
            self.assertIn(marker, capture)
        self.assertIn("materially larger capture region", _flat(capture))
        self.assertIn("must not insert the \\SI{120}{s} equilibrium covariance warm-up", _flat(capture))

    def test_proxy_targets_capture_region_not_local_basin(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        startup = _read("w3d-init.tex-part")
        for marker in (
            r"\label{eq:semiglobal-handoff-metric-radius}",
            r"\label{eq:semiglobal-capture-entry}",
            r"R_H^{\max}<R_C",
            r"\ref{thm:finite-live-capture}",
            r"\label{thm:semiglobal-proxy-live}",
        ):
            self.assertIn(marker, proof)
        self.assertNotIn(r"\label{eq:semiglobal-basin-entry}", proof)
        self.assertIn(
            r"does \emph{not} assert entry into the final local ISS tube",
            _flat(proof),
        )
        self.assertIn(r"R_H^{\max}<R_C", startup)
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

    def test_block_local_structure_is_retained_after_capture(self):
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
        self.assertIn("Theorem~\\ref{thm:iss-block-local} applies", proof)

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
