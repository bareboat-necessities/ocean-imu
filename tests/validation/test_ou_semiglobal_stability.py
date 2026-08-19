"""Semantic contract for the OU--III semiglobal stability extension."""

import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"


def _read(name: str) -> str:
    return (DOC / name).read_text(encoding="utf-8")


class OUIIISemiglobalStabilityContractTests(unittest.TestCase):
    def test_main_article_wires_semiglobal_extension_after_live_iss(self):
        main = _read("kalman_ou-w3d.tex")
        local = r"\input{w3d-iss-stability.tex-part}"
        semiglobal = r"\input{w3d-semiglobal-stability.tex-part}"
        adaptation = r"\input{w3d-adaptation-motivation.tex-part}"
        self.assertIn(semiglobal, main)
        self.assertLess(main.index(local), main.index(semiglobal))
        self.assertLess(main.index(semiglobal), main.index(adaptation))

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
            "antipodal",
            r"\cite{Mahony2008NonlinearComplementary}",
            r"\label{thm:semiglobal-proxy-live}",
            r"\ref{thm:iss-21-local}",
            r"\label{eq:semiglobal-basin-entry}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("does not claim semiglobal convergence", proof)

    def test_implemented_handoff_and_timeout_are_part_of_scope(self):
        proof = _read("w3d-semiglobal-stability.tex-part")
        startup = _read("w3d-init.tex-part")
        self.assertIn("0.075", proof)
        self.assertIn(r"\SI{150}{s}", proof)
        self.assertIn(r"\eqref{eq:semiglobal-basin-entry}", startup)
        self.assertIn(r"\ref{thm:semiglobal-proxy-live}", startup)


if __name__ == "__main__":
    unittest.main()
