"""Semantic/source contract for direct regional and hybrid recovery proofs."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
SRC = REPO_ROOT / "src" / "kalman_ou_iii"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return re.sub(r"\s+", " ", text)


class OUIIIStabilityPhaseBContractTests(unittest.TestCase):
    def test_capture_propagates_actual_handoff_set(self):
        capture = _read(DOC / "w3d-finite-live-capture.tex-part")
        flat = _flat(capture)
        for marker in (
            r"\label{eq:capture-set-envelope}",
            r"X_0=\mathcal H_{k_H}",
            r"\label{eq:capture-comparison-sequence}",
            r"\label{eq:capture-prefix-domain}",
            r"\label{eq:capture-finite-entry}",
            r"\label{thm:finite-live-capture}",
        ):
            self.assertIn(marker, capture)
        self.assertIn("compact handoff set actually produced by the source", flat)
        self.assertIn("propagation of $x$ itself is the certificate object", flat.casefold())
        self.assertNotIn(r"\vct R_H\prec\vct R_C", capture)
        self.assertNotIn(r"\vct R_C", capture)

    def test_radius_box_is_optional_enclosure(self):
        capture = _read(DOC / "w3d-finite-live-capture.tex-part")
        flat = _flat(capture)
        for marker in (
            r"\label{eq:capture-structured-radii}",
            r"\label{eq:capture-radius-enclosure}",
        ):
            self.assertIn(marker, capture)
        self.assertIn("optional outer approximation", flat)
        self.assertIn("need not be attained by the same physical state", flat)

    def test_certificate_contract_starts_at_exact_go_live_mode(self):
        capture = _read(DOC / "w3d-finite-live-capture.tex-part")
        self.assertIn(r"\texttt{goLive()}", capture)
        self.assertIn("accelerometer-bias hold", capture)
        self.assertIn(r"\label{eq:capture-certificate-pass}", capture)
        self.assertIn(r"\mathcal H_{k_H}\subseteq\mathcal D_{k_H}", capture)

    def test_hybrid_uses_jump_then_recovery(self):
        hybrid = _read(DOC / "w3d-hybrid-stability.tex-part")
        for marker in (
            r"\label{eq:hybrid-jump-map}",
            r"\label{eq:hybrid-mag-envelope}",
            r"\label{eq:hybrid-tilt-force-bound}",
            r"\label{eq:hybrid-jump-entry}",
            r"\mathfrak J_j(X^-)\subseteq\mathcal D",
            r"\mathfrak T_{k^+,N_j}",
            r"\label{thm:hybrid-live-recovery}",
        ):
            self.assertIn(marker, hybrid)
        self.assertNotIn("returns directly to the outer capture set", hybrid)

    def test_repeated_relocks_have_dwell_small_gain_bound(self):
        hybrid = _read(DOC / "w3d-hybrid-stability.tex-part")
        for marker in (
            r"T_{\rm cool}=3\ \mathrm{s}",
            r"\label{eq:hybrid-cycle-bounds}",
            r"\mu_T:=\kappa_T\lambda_F^{N_{\rm cool}}<1",
            r"\label{eq:hybrid-recurrent-bound}",
            "excludes Zeno",
        ):
            self.assertIn(marker, hybrid)

    def test_quotient_has_regional_capture(self):
        hybrid = _read(DOC / "w3d-hybrid-stability.tex-part")
        for marker in (
            r"\label{sec:hybrid-quotient-capture}",
            r"\label{eq:hybrid-quotient-metric}",
            r"\label{thm:hybrid-quotient-capture}",
            r"\mathcal H_k^\perp",
            "No absolute-yaw convergence is asserted",
        ):
            self.assertIn(marker, hybrid)

    def test_source_tilt_relock_still_matches_proof(self):
        src = _read(SRC / "SeaStateFusionFilter_OU_III.h")
        for token in (
            "constexpr float TILT_RESET_DEG = 70.0f;",
            "constexpr float TILT_RESET_HOLD_SEC = 0.35f;",
            "constexpr float TILT_RESET_COOLDOWN_SEC = 3.0f;",
            "initialize_from_acc_preserve_yaw(acc)",
        ):
            self.assertIn(token, src)


if __name__ == "__main__":
    unittest.main()
