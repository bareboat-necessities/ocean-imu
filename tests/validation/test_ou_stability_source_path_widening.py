"""Publication contract for source-reachable/path-dependent OU--III stability widening."""

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
    def test_source_path_refinement_is_wired_before_phase_d(self):
        phase_d = _read(DOC / "w3d-stability-widening-phase-d.tex-part")
        self.assertTrue(
            phase_d.startswith(r"\input{w3d-stability-widening-source-path.tex-part}")
        )

    def test_source_family_is_explicit_subset_of_broad_envelope(self):
        proof = _read(DOC / "w3d-stability-widening-source-path.tex-part")
        for marker in (
            r"\label{eq:widen-src-subset}",
            r"\label{eq:widen-src-kappa}",
            r"\label{eq:widen-src-monotone}",
            r"\label{thm:widen-src-kappa-ues}",
            r"\label{eq:widen-src-kappa-condition}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn(r"\Pi^{{\rm src},m}_{k,N}\subseteq\Pi^{{\rm env},m}_{k,N}", proof)
        self.assertIn(r"\kappa_{N,m}^{\rm src}\le\kappa_{N,m}^{\rm env}", proof)
        self.assertIn("Individual samples inside the $N$-step window may have norm gain greater than one", flat)
        self.assertIn("additive refinement", flat)

    def test_mode_dimension_change_is_not_hidden_in_matrix_product(self):
        flat = _flat(_read(DOC / "w3d-stability-widening-source-path.tex-part"))
        self.assertIn("different active coordinates (18 and 21 states)", flat)
        self.assertIn("applied separately to each fixed-dimensional normal-Live mode", flat)
        self.assertIn("transition that changes the active state dimension is a boundary between such graphs", flat)
        self.assertIn("held-to-active bias-release transition must place the post-transition error", flat)

    def test_path_complete_metric_retains_source_consistency(self):
        proof = _read(DOC / "w3d-stability-widening-source-path.tex-part")
        for marker in (
            r"\label{eq:widen-src-graph}",
            r"\label{eq:widen-src-edge-family}",
            r"\label{eq:widen-src-word-family}",
            r"\label{eq:widen-src-word-length}",
            r"\label{eq:widen-path-P-bounds}",
            r"\label{eq:widen-path-lmi}",
            r"\label{thm:widen-path-ues}",
            r"\label{eq:widen-path-generalized-gain}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn("source complete", flat)
        self.assertIn("complete tuple is source-reachable", flat)
        self.assertIn("No one-step decrease is required inside a word", flat)
        self.assertIn("checking finitely many nominal trajectories alone is not a proof", flat.lower())

    def test_path_metric_contains_old_routes_as_special_cases(self):
        flat = _flat(_read(DOC / "w3d-stability-widening-source-path.tex-part"))
        self.assertIn("Phase-B Euclidean $N$-step certificate is contained as the special choice", flat)
        self.assertIn("common quadratic metric is contained as the special choice", flat)
        self.assertIn("failure of a common one-step or arbitrary-envelope certificate does not imply failure", flat)

    def test_constructive_constants_avoid_old_M_relaxation(self):
        proof = _read(DOC / "w3d-stability-widening-source-path.tex-part")
        for marker in (
            r"\label{eq:widen-path-prefix}",
            r"\label{eq:widen-path-ues-bound}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn(r"\kappa_N\rightarrow\rho\rightarrow M\rightarrow p_+", proof)
        self.assertIn("does not require the conservative intermediate chain", flat)

    def test_lifted_nonlinear_route_allows_transient_positive_delta_v(self):
        proof = _read(DOC / "w3d-stability-widening-source-path.tex-part")
        for marker in (
            r"\label{eq:widen-path-nonlinear-lift}",
            r"\label{eq:widen-path-direct-margin}",
            r"\label{thm:widen-path-direct-iss}",
            r"\label{eq:widen-path-direct-condition}",
        ):
            self.assertIn(marker, proof)
        flat = _flat(proof)
        self.assertIn("Positive one-step Lyapunov increments are permitted inside a word", flat)
        self.assertIn("full $S=0$ cross-gain", flat)
        self.assertIn("additive widening route", flat)

    def test_verification_contract_requires_rigorous_enclosures(self):
        proof = _read(DOC / "w3d-stability-widening-source-path.tex-part")
        self.assertIn(r"\label{eq:widen-path-contract}", proof)
        flat = _flat(proof)
        self.assertIn("every reachable same-mode source transition is covered", flat)
        self.assertIn("complete continuous source-reachable family", flat)
        self.assertIn("Monte Carlo trajectories and dense sampling are diagnostics only", flat)
        self.assertIn("changes neither the estimator nor the adaptation law", flat)
        self.assertIn("does not weaken the arbitrary-envelope theorem", flat)


if __name__ == "__main__":
    unittest.main()
