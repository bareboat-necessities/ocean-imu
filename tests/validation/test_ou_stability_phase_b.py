"""Semantic/source contract for OU--III stability hardening Phase B."""

import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
SRC = REPO_ROOT / "src" / "kalman_ou_iii"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return " ".join(text.split())


class OUIIIStabilityPhaseBContractTests(unittest.TestCase):
    def test_hybrid_proof_is_wired_into_stability_chain(self):
        proof = _read(DOC / "w3d-semiglobal-stability.tex-part")
        init = _read(DOC / "w3d-init.tex-part")
        conclusion = _read(DOC / "w3d-conclusion-summary.tex-part")

        self.assertIn(r"\input{w3d-hybrid-stability.tex-part}", proof)
        self.assertIn(r"\ref{thm:hybrid-live-recovery}", proof)
        self.assertIn(r"\ref{thm:hybrid-live-recovery}", init)
        self.assertIn(r"\ref{thm:hybrid-live-recovery}", conclusion)
        self.assertNotIn("deferred to the next proof phase", init)
        self.assertNotIn("deferred to the next stability-hardening phase", conclusion)

    def test_hybrid_theorem_contains_both_exact_jump_families(self):
        hybrid = _read(DOC / "w3d-hybrid-stability.tex-part")
        for marker in (
            r"\label{eq:hybrid-jump-map}",
            r"\label{eq:hybrid-mag-attitude-map}",
            r"\label{eq:hybrid-mag-heading-residual}",
            r"\label{eq:hybrid-mag-envelope}",
            r"\label{eq:hybrid-tilt-attitude-map}",
            r"\label{eq:hybrid-tilt-force-bound}",
            r"\label{eq:hybrid-tilt-reset-bound}",
            r"\label{eq:hybrid-tilt-yaw-chart}",
            r"\label{eq:hybrid-tilt-yaw-preserved}",
            r"\label{eq:hybrid-tilt-dwell}",
            r"\label{eq:hybrid-jump-entry}",
            r"\label{thm:hybrid-live-recovery}",
            r"\label{eq:hybrid-cooldown-capture-condition}",
        ):
            self.assertIn(marker, hybrid)

        folded = hybrid.casefold()
        self.assertIn("not asserted numerically here", folded)
        self.assertIn("not convergence to zero", folded)
        self.assertIn("continuous hard-iron correction", folded)

    def test_magnetic_family_covers_late_first_gauge_and_refinement_source(self):
        src = _read(SRC / "SeaStateFusionFilter_OU_III.h")
        for token in (
            # Provisional acquisition can complete after an ungauged timeout.
            "if (usingProxyInit_() && stage_ != Stage::Live)",
            "pending_yaw_abs_rad_ = yaw_abs_rad;",
            "boatQuatWithAbsoluteYaw_(q_bw, yaw_abs_rad)",
            "impl_.mekf().set_quaternion_boat(q_new);",
            "mag_ref_set_ = true;",
            "mag_north_lock_time_sec_ = t_;",
            "syncLinearBlockGate_();",
            # Second-stage refinement is a distinct one-time member.
            "void maybeRefineMagReference_",
            "if (mag_refine_done_) return;",
            "mag_refine_done_    = true;",
            "impl_.setAccBiasHold(false);",
            "maybeApplyContinuousHardIron_();",
        ):
            self.assertIn(token, src)

        hybrid = _flat(_read(DOC / "w3d-hybrid-stability.tex-part"))
        for marker in (
            "First Live gauge acquisition",
            "Second-stage refinement",
            r"\rho_M\in\{0,1\}",
            "first-gauge branch",
            r"\texttt{mag\_ref\_set\_}",
            r"\texttt{mag\_refine\_done\_}",
            "writes no attitude state or covariance",
            r"not to $\mathcal R_M$",
        ):
            self.assertIn(marker, hybrid)

    def test_tilt_relock_contract_matches_wrapper_source(self):
        src = _read(SRC / "SeaStateFusionFilter_OU_III.h")
        for token in (
            "constexpr float TILT_RESET_DEG = 70.0f;",
            "constexpr float TILT_RESET_HOLD_SEC = 0.35f;",
            "constexpr float TILT_RESET_COOLDOWN_SEC = 3.0f;",
            "mekf_->initialize_from_acc_preserve_yaw(acc);",
            "tilt_reset_cooldown_sec_ = TILT_RESET_COOLDOWN_SEC;",
        ):
            self.assertIn(token, src)

        hybrid = _read(DOC / "w3d-hybrid-stability.tex-part")
        self.assertIn(r"T_{\rm cool}=3\ \mathrm{s}", hybrid)
        self.assertIn("does not itself impose", hybrid)

    def test_tilt_relock_contract_matches_mekf_reset_source(self):
        src = _read(SRC / "Kalman3D_Wave_OU_III.h")
        for token in (
            "initialize_from_acc_preserve_yaw",
            "const T yaw_old = std::atan2",
            "initialize_from_acc(acc_body);",
            "q_new_bw = q_yaw * q_pitch * q_roll;",
            "set_quaternion_boat(q_new_bw);",
            "set_accel_only_attitude_covariance_();",
            "zero_AL_cross_cov_once_();",
        ):
            self.assertIn(token, src)

        self.assertIn("Pext.template block<3,3>(0,3).setZero();", src)
        self.assertIn("Pext.template block<3,3>(0, OFF_BA).setZero();", src)

    def test_repeated_reset_claim_is_safety_not_false_convergence(self):
        hybrid = _flat(_read(DOC / "w3d-hybrid-stability.tex-part")).casefold()
        self.assertIn("if the hard-jump sequence is finite", hybrid)
        self.assertIn("if tilt re-locks occur infinitely often", hybrid)
        self.assertIn("prevents zeno behavior", hybrid)
        self.assertIn("recurrent capture claim", hybrid)
        self.assertIn("resets eventually cease", hybrid)
        self.assertIn("two source magnetic jumps are each one-time events", hybrid)


if __name__ == "__main__":
    unittest.main()
