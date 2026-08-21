"""Semantic contract for the OU--III semiglobal stability extension."""

import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
SRC = REPO_ROOT / "src"


def _read(name: str) -> str:
    return (DOC / name).read_text(encoding="utf-8")


def _flat(text: str) -> str:
    """Line breaks in the sources are wrapping, not content."""
    return " ".join(text.split())


def _source(relative: str) -> str:
    return (SRC / relative).read_text(encoding="utf-8")


class OUIIISemiglobalStabilityContractTests(unittest.TestCase):
    def test_main_article_wires_semiglobal_extension_after_live_iss(self):
        main = _read("kalman_ou-w3d.tex")
        local = r"\input{w3d-iss-stability.tex-part}"
        semiglobal = r"\input{w3d-semiglobal-stability.tex-part}"
        evidence = r"\input{w3d-sim-charts.tex-part}"
        self.assertIn(semiglobal, main)
        self.assertLess(main.index(local), main.index(semiglobal))
        self.assertLess(main.index(semiglobal), main.index(evidence))

    def test_extension_uses_spectral_mse_schedule_not_cubic_adaptation(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\label{eq:semiglobal-mse-schedule}",
            r"\widehat\sigma_{a,B}^{\,6/7}",
            r"\tau^{41/14}",
            r"\label{eq:semiglobal-mse-powers}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("no $\\sigma\\tau^3$ adaptation relation", proof)

    def test_extension_is_explicitly_almost_global_not_global_ekf(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\mathcal K_{\epsilon,B_\beta}",
            r"\mathcal N_P:=W^s(\mathcal U_P)",
            "antipodal",
            r"\cite{Mahony2008NonlinearComplementary}",
            r"\label{thm:semiglobal-proxy-live}",
            r"\ref{thm:iss-21-local}",
            r"\label{eq:semiglobal-basin-entry}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("does not claim semiglobal convergence", _flat(proof))

    def test_phase1_uses_block_structured_dimensionless_live_basin(self):
        block = _read("w3d-block-local-iss.tex-part")
        semiglobal = _read("w3d-semiglobal-stability.tex-part")

        self.assertTrue(
            semiglobal.startswith(r"\input{w3d-block-local-iss.tex-part}"),
            "block-local refinement must be inserted before Section XII",
        )
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

        self.assertIn("dimensionless", block)
        self.assertIn("post-prediction, pre-correction", block)
        self.assertIn("uniform bounds on both", block)
        self.assertIn("no direct dependence on", block)
        self.assertIn("does \\emph{not} certify arbitrarily large", block)
        self.assertIn("Nor does it remove $S$ from", block)

        for marker in (
            r"B_{\xi,H}",
            r"B_{\ell,H}",
            r"M_{\xi\ell}B_{\ell,H}",
            r"\overline d_L",
            r"\ref{thm:iss-block-local}",
        ):
            self.assertIn(marker, semiglobal)
        self.assertIn("replaces the earlier unscaled 21-state Euclidean-ball condition", semiglobal)
        self.assertIn("$S$ remains in the sensitive block", semiglobal)

    def test_implemented_handoff_and_timeout_are_part_of_scope(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        startup = _read("w3d-init.tex-part")
        self.assertIn("0.075", proof)
        self.assertIn(r"\SI{150}{s}", proof)
        self.assertIn(r"\label{eq:semiglobal-aligned-branch}", proof)
        self.assertIn("sine-only residual is nevertheless not by itself", proof)
        self.assertIn(r"\eqref{eq:semiglobal-aligned-branch}", startup)
        self.assertIn(r"\eqref{eq:semiglobal-basin-entry}", startup)
        self.assertIn(r"\ref{thm:semiglobal-proxy-live}", startup)

    def test_aligned_branch_is_gated_and_not_merely_assumed(self):
        """The sine residual is blind to a 180 deg flip, so the branch has to
        be tested where the handoff is accepted -- including on the forced
        timeout path, which never looks at the residual at all."""
        proof = _flat(_read("w3d-semiglobal-stability.tex-part"))
        startup = _flat(_read("w3d-init.tex-part"))

        self.assertIn("sign of the inner product", proof)
        self.assertIn("sign of the inner product", startup)
        self.assertIn("withholds a timeout-forced handoff until", proof)
        self.assertIn("delayed rather than accepted on the antipodal branch", startup)

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
            # The quality hold and the 150 s forced handoff are both branch-gated.
            self.assertIn(
                "const bool tilt_trusted = mag_gravity_aligned_branch_ &&", gate)
            self.assertIn(
                "const bool ready_by_timeout = proxy_ready && (t_ >= timeout_sec) "
                "&& mag_gravity_aligned_branch_;", gate)


if __name__ == "__main__":
    unittest.main()
