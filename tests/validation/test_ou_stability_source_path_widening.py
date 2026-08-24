"""Publication contract for the sole source-reachable path OU--III certificate."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", text)


class OUIIISourcePathStabilityWideningContractTests(unittest.TestCase):
    def test_source_path_is_wired_before_semiglobal_composition(self):
        main = _read(DOC / "kalman_ou-w3d.tex")
        path_input = r"\input{w3d-stability-widening-source-path.tex-part}"
        semiglobal_input = r"\input{w3d-semiglobal-stability.tex-part}"
        phase_d_input = r"\input{w3d-stability-widening-phase-d.tex-part}"
        self.assertIn(path_input, main)
        self.assertLess(main.index(path_input), main.index(semiglobal_input))
        self.assertLess(main.index(semiglobal_input), main.index(phase_d_input))
        phase_d = _read(DOC / "w3d-stability-widening-phase-d.tex-part")
        self.assertNotIn("w3d-stability-widening-source-path", phase_d)

    def test_fixed_dimension_source_graph_is_explicit(self):
        proof = _read(DOC / "w3d-stability-widening-source-path.tex-part")
        flat = _flat(proof)
        for marker in (
            r"\label{eq:widen-src-graph}",
            r"\label{eq:widen-src-edge-family}",
            r"\label{eq:widen-src-word-family}",
            r"\label{eq:widen-src-word-length}",
        ):
            self.assertIn(marker, proof)
        self.assertIn("different active coordinates (18 and 21 states)", flat)
        self.assertIn("applied separately to each fixed-dimensional normal-Live mode", flat)
        self.assertIn("jointly source reachable", flat)
        self.assertIn("no square homogeneous transition is formed across a state-dimension change", flat.casefold())

    def test_path_metric_is_primary_linear_certificate(self):
        proof = _read(DOC / "w3d-stability-widening-source-path.tex-part")
        for marker in (
            r"\label{eq:widen-path-Pbar}",
            r"\label{eq:widen-path-P-bounds}",
            r"\label{eq:widen-path-lmi}",
            r"\label{thm:widen-path-ues}",
            r"\label{eq:widen-path-generalized-gain}",
            r"\label{eq:widen-path-prefix}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn("No decrease is required at an individual sample inside the word", flat)
        self.assertNotIn(r"\kappa_{N,m}^{\rm src}", proof)
        self.assertNotIn("special choice", flat)
        for retired in ("widen-cylinder-iss", "widen-component-iss", "widen-direct-iss"):
            self.assertNotIn(retired, proof)
        self.assertIn("sole quantitative normal-Live certificate path", flat)

    def test_linear_metric_is_local_quadratic_of_group_metric(self):
        proof = _read(DOC / "w3d-stability-widening-source-path.tex-part")
        flat = _flat(proof)
        self.assertIn(r"\frac{a_{R,i}}{2}\mat I_3", proof)
        self.assertIn(r"\mat P_{\xi,i}", proof)
        self.assertIn("precisely the local quadratic of the exact group metric", flat)
        self.assertIn("no attitude--linear cross term", flat.casefold())
        self.assertIn("using the same $a_{r,i}$ and $\\mat p_{\\xi,i}$", flat.casefold())

    def test_nonlinear_metric_uses_group_energy_and_exact_words(self):
        proof = _read(DOC / "w3d-stability-widening-source-path.tex-part")
        for marker in (
            r"\label{eq:widen-path-group-metric}",
            r"\label{eq:widen-path-nonlinear-lift}",
            r"\label{eq:widen-path-direct-margin}",
            r"\label{thm:widen-path-direct-iss}",
            r"\label{eq:widen-path-direct-condition}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn("full $S=0$ attitude cross-gain", flat)
        self.assertIn("exact finite-correction identity", flat)
        self.assertIn("positive one-sample changes", flat.casefold())

    def test_verification_contract_requires_complete_continuous_enclosures(self):
        proof = _read(DOC / "w3d-stability-widening-source-path.tex-part")
        self.assertIn(r"\label{eq:widen-path-contract}", proof)
        flat = _flat(proof)
        self.assertIn("every reachable same-mode transition is covered", flat)
        self.assertIn("complete continuous jointly source-reachable family", flat)
        self.assertIn("Monte Carlo and dense trajectory sampling are diagnostics only", flat)
        self.assertIn("sole quantitative normal-Live certificate path", flat)
        self.assertIn("changes neither the estimator nor the adaptation law", flat)


if __name__ == "__main__":
    unittest.main()
