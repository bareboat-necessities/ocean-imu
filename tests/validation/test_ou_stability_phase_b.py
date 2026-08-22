"""Semantic/source contract for OU--III stability hardening Phase B."""

import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
SRC = REPO_ROOT / "src" / "kalman_ou_iii"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


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

        self.assertIn("not asserted numerically here", hybrid)
        self.assertIn("not convergence to zero", hybrid)
        self.assertIn("continuous hard-iron correction", hybrid)

    def test_magnetic_refinement_contract_matches_source(self):
        src = _read(SRC / "SeaStateFusionFilter_OU_III.h")
        for token in (
            "void maybeRefineMagReference_",
            "if (mag_refine_done_) return;",
            "setMagWorldRef_(mag_world_ref_uT);",
            "boatQuatWithAbsoluteYaw_(q_bw, yaw_abs_rad)",
            "impl_.mekf().set_quaternion_boat(q_new);",
            "mag_refine_done_    = true;",
            "impl_.setAccBiasHold(false);",
            "maybeApplyContinuousHardIron_();",
        ):
            self.assertIn(token, src)

        hybrid = _read(DOC / "w3d-hybrid-stability.tex-part")
        self.assertIn("one-shot by source construction", hybrid)
        self.assertIn("writes no attitude state or covariance", hybrid)

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
        self.assertIn("does not itself impose this", hybrid)

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

        # initialize_from_acc() clears the attitude cross-covariances with
        # gyro bias and accelerometer bias before preserve-yaw rewrites yaw.
        self.assertIn("Pext.template block<3,3>(0,3).setZero();", src)
        self.assertIn("Pext.template block<3,3>(0, OFF_BA).setZero();", src)

    def test_repeated_reset_claim_is_safety_not_false_convergence(self):
        hybrid = _read(DOC / "w3d-hybrid-stability.tex-part")
        self.assertIn("If the hard-jump sequence is finite", hybrid)
        self.assertIn("if tilt re-locks occur infinitely often", hybrid)
        self.assertIn("prevents Zeno behavior", hybrid)
        self.assertIn("recurrent capture claim", hybrid)
        self.assertIn("resets eventually cease", hybrid)


if __name__ == "__main__":
    unittest.main()
