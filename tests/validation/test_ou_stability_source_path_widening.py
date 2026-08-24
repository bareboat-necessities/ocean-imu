"""Publication contract for the sole source-reachable OU--III Live certificate."""

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
        main=_read(DOC/"kalman_ou-w3d.tex")
        path_input=r"\input{w3d-stability-widening-source-path.tex-part}"
        semiglobal_input=r"\input{w3d-semiglobal-stability.tex-part}"
        phase_d_input=r"\input{w3d-stability-widening-phase-d.tex-part}"
        self.assertIn(path_input,main); self.assertLess(main.index(path_input),main.index(semiglobal_input)); self.assertLess(main.index(semiglobal_input),main.index(phase_d_input))

    def test_fixed_dimension_joint_source_graph_is_explicit(self):
        proof=_read(DOC/"w3d-stability-widening-source-path.tex-part"); flat=_flat(proof)
        for marker in (r"\label{eq:widen-src-graph}",r"\label{eq:widen-src-edge-family}",r"\label{eq:widen-src-word-family}",r"\label{eq:widen-src-word-length}"):
            self.assertIn(marker,proof)
        self.assertIn("18 and 21",flat); self.assertIn("jointly source reachable",flat); self.assertIn("fixed dimensions",flat)

    def test_homogeneous_certificate_is_complete_word_information_not_one_step(self):
        proof=_read(DOC/"w3d-stability-widening-source-path.tex-part"); flat=_flat(proof)
        for marker in (r"\label{eq:widen-path-information-word}",r"\label{eq:widen-path-prefix}",r"\label{eq:widen-path-information-metric}",r"\label{eq:widen-path-P-bounds}",r"\label{thm:widen-path-ues}"):
            self.assertIn(marker,proof)
        self.assertIn("No decrease is required at an individual sample inside the word",flat)
        self.assertIn("four-$S$",proof); self.assertIn("not a theorem-promotion fallback",flat)

    def test_cayley_metric_is_exact_information_geometry_with_cross_terms(self):
        proof=_read(DOC/"w3d-stability-widening-source-path.tex-part"); flat=_flat(proof)
        for marker in (r"\label{eq:widen-path-cayley-coordinate}",r"\label{eq:widen-path-group-metric}",r"\label{eq:widen-path-cayley-product}"):
            self.assertIn(marker,proof)
        self.assertIn(r"2\tan\!\left(\frac{\theta_e}{2}\right)u_e",proof)
        self.assertIn(r"s_m z_C^\top\Sigma_i^{-1}z_C",proof)
        self.assertIn("retains all source information cross terms",flat)
        self.assertIn("single mode-global positive scalar",flat)
        self.assertNotIn(r"\operatorname{blkdiag}",proof)
        self.assertNotIn(r"P_{\xi,i}",proof)

    def test_exact_shipping_nonlinear_word_and_prefix_bootstrap_are_explicit(self):
        proof=_read(DOC/"w3d-stability-widening-source-path.tex-part"); flat=_flat(proof)
        for marker in (r"\label{eq:widen-path-source-series-quaternion}",r"\label{eq:widen-path-nonlinear-lift}",r"\label{eq:widen-path-quadratic-defect}",r"\label{eq:widen-path-word-defect}",r"\label{eq:widen-path-W-star}",r"\label{eq:widen-path-direct-decrease}"):
            self.assertIn(marker,proof)
        self.assertIn("due $S=0$ correction occurs",flat)
        self.assertIn("before the accelerometer correction",flat)
        self.assertIn("each accepted correction performs its own quaternion injection",flat)
        self.assertIn("full $S\\to\\theta$ gain is retained",flat)
        self.assertIn("without enumerating all rejection strings",flat)

    def test_direct_positive_gap_and_mu_are_not_lost_to_binary64_one_minus(self):
        proof=_read(DOC/"w3d-stability-widening-source-path.tex-part"); flat=_flat(proof)
        for marker in (r"\label{eq:widen-path-direct-margin}",r"\label{eq:widen-path-direct-condition}",r"\label{thm:widen-path-direct-iss}"):
            self.assertIn(marker,proof)
        self.assertIn("far below binary64 machine epsilon",flat)
        self.assertIn("does not evaluate $1-\\delta_m/2$",flat)
        self.assertIn(r"\mu_W\ge\frac{\delta_m}{2}m_->0",proof)
        self.assertIn("Positive one-sample changes",flat)

    def test_A_projection_surface_is_not_silently_linearized(self):
        flat=_flat(_read(DOC/"w3d-stability-widening-source-path.tex-part"))
        self.assertIn("0.45",flat); self.assertIn("0.5",flat); self.assertIn("identity interior branch",flat)
        self.assertIn("outer capture/hybrid obligation",flat)

    def test_verification_contract_requires_source_complete_outward_enclosure(self):
        proof=_read(DOC/"w3d-stability-widening-source-path.tex-part"); flat=_flat(proof)
        self.assertIn(r"\label{eq:widen-path-contract}",proof)
        self.assertIn("every reachable same-mode transition",flat)
        self.assertIn("outward-rounded",flat)
        self.assertIn("Monte Carlo and dense trajectory sampling are diagnostics only",flat)
        self.assertIn("changes neither the estimator nor the adaptation law",flat)
        self.assertIn("sole quantitative normal-Live certificate path",flat)


if __name__=="__main__":
    unittest.main()
