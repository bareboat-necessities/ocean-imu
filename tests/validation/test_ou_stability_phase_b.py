"""Semantic/source contract for funnel capture and hybrid recovery proofs."""

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
    def test_capture_starts_from_actual_handoff_set_and_funnel_level(self):
        capture = _read(DOC / "w3d-finite-live-capture.tex-part")
        flat = _flat(capture)
        for marker in (
            r"\label{eq:capture-set-envelope}",
            r"X_0:=\mathcal H_{k_H}",
            r"\label{eq:capture-initial-level}",
            r"\label{eq:capture-handoff-inclusion}",
            r"\label{eq:capture-prefix-domain}",
            r"\label{eq:capture-funnel-recursion}",
            r"\label{thm:finite-live-capture}",
        ):
            self.assertIn(marker, capture)
        self.assertIn("exact compact handoff", flat)
        self.assertIn("propagate only verified scalar funnel levels", flat.casefold())
        self.assertNotIn("radius-only implementation", flat)
        self.assertNotIn("optional outer approximation", flat)

    def test_radius_summary_is_diagnostic_only(self):
        capture = _read(DOC / "w3d-finite-live-capture.tex-part")
        flat = _flat(capture).casefold()
        self.assertIn(r"\label{eq:capture-radius-enclosure}", capture)
        self.assertIn("for reporting only", flat)
        self.assertIn("are not propagated and are not certificate state variables", flat)

    def test_certificate_contract_starts_at_exact_go_live_mode(self):
        capture = _read(DOC / "w3d-finite-live-capture.tex-part")
        self.assertIn(r"\texttt{goLive()}", capture)
        self.assertIn("accelerometer-bias hold", capture)
        self.assertIn(r"\label{eq:capture-certificate-pass}", capture)
        self.assertIn(r"\mathcal H_{k_H}\subseteq\bigcup_i\mathcal D_i", capture)
        self.assertIn(r"c_{N_H,i}\le b_i", capture)

    def test_hybrid_uses_exact_jump_to_same_funnel_metric(self):
        hybrid = _read(DOC / "w3d-hybrid-stability.tex-part")
        flat = _flat(hybrid)
        for marker in (
            r"\label{eq:hybrid-jump-map}",
            r"\label{eq:hybrid-mag-envelope}",
            r"\label{eq:hybrid-tilt-force-bound}",
            r"\label{eq:hybrid-jump-entry}",
            r"\mathfrak J_j(\mathcal F_i(c_i))\subseteq\mathcal D_j",
            r"W_j(\mathcal R_j(e,\pi;\mu_j))",
            r"\label{thm:hybrid-live-recovery}",
        ):
            self.assertIn(marker, hybrid)
        self.assertIn("same scalar funnel recursion", flat)
        self.assertNotIn(r"\mathfrak T_{k^+,N_j}", hybrid)
        self.assertNotIn(r"\lambda_F", hybrid)

    def test_recurrent_relocks_are_owned_by_phase_e(self):
        hybrid = _read(DOC / "w3d-hybrid-stability.tex-part")
        phase_e = _read(DOC / "w3d-stability-widening-phase-e.tex-part")
        self.assertIn(r"T_{\rm cool}=3\ \mathrm{s}", hybrid)
        self.assertIn("Recurrent-reset small gain is derived once", hybrid)
        for marker in (
            r"\label{eq:widen-reset-to-set}",
            r"\label{eq:widen-cooldown-factor}",
            r"\label{eq:widen-hybrid-mu}",
            r"\label{eq:widen-reset-ultimate}",
            r"\label{thm:widen-hybrid-small-gain}",
        ):
            self.assertIn(marker, phase_e)
        self.assertNotIn(r"\label{eq:hybrid-cycle-bounds}", hybrid)
        self.assertNotIn(r"\label{eq:hybrid-recurrent-bound}", hybrid)

    def test_quotient_uses_same_source_funnel_construction(self):
        hybrid = _read(DOC / "w3d-hybrid-stability.tex-part")
        flat = _flat(hybrid)
        for marker in (
            r"\label{sec:hybrid-quotient-capture}",
            r"\label{eq:hybrid-quotient-metric}",
            r"\label{eq:hybrid-quotient-word}",
            r"\label{thm:hybrid-quotient-capture}",
            r"\mathcal H^\perp",
            "No convergence of absolute yaw is asserted",
        ):
            self.assertIn(marker, hybrid)
        self.assertIn("source-complete funnel recursion", flat)
        self.assertNotIn("finite capture modulo yaw", flat)

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
