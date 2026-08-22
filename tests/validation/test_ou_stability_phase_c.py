"""Semantic/source contract for OU--III stability hardening Phase C."""

import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
SRC = REPO_ROOT / "src"


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _flat(text: str) -> str:
    return " ".join(text.split())


class OUIIIStabilityPhaseCContractTests(unittest.TestCase):
    def test_explicit_mahony_lyapunov_replaces_converse_existence(self):
        proof = _read(DOC / "w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\label{eq:semiglobal-proxy-gains}",
            r"\label{eq:semiglobal-proxy-reduced-dynamics}",
            r"\label{eq:semiglobal-proxy-energy}",
            r"\label{eq:semiglobal-proxy-lyapunov}",
            r"\label{eq:semiglobal-proxy-delta-condition}",
            r"\label{eq:semiglobal-proxy-rate}",
            r"\label{eq:semiglobal-proxy-input-gain}",
            r"\label{eq:semiglobal-proxy-iss-differential}",
            r"\label{eq:semiglobal-proxy-iss}",
            r"\label{eq:semiglobal-proxy-explicit-bounds}",
            r"\label{eq:semiglobal-proxy-beta-gamma}",
            r"\label{thm:semiglobal-proxy-explicit}",
        ):
            self.assertIn(marker, proof)

        self.assertNotIn("A converse Lyapunov theorem", proof)
        self.assertNotIn(r"\alpha_{1,\epsilon}", proof)
        self.assertNotIn(r"\alpha_{3,\epsilon}", proof)
        self.assertIn(r"\delta_P=\SI{5}{s}", proof)
        self.assertIn(r"q_P=0.5", proof)
        self.assertIn(r"a_P=0.025\ \mathrm{s}^{-1}", proof)

    def test_proxy_gain_and_sign_contract_matches_source(self):
        wrapper = _read(SRC / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h")
        vertical = _read(SRC / "tuner" / "VerticalAccelComplementary.h")
        mahony = _read(SRC / "ahrs" / "Mahony_AHRS.h")
        proof = _read(DOC / "w3d-semiglobal-stability.tex-part")

        self.assertIn("constexpr float STARTUP_PROXY_TWO_KP_DEFAULT = 0.2f;", wrapper)
        self.assertIn("constexpr float STARTUP_PROXY_TWO_KI_DEFAULT = 0.02f;", wrapper)
        self.assertIn("STARTUP_PROXY_TWO_KP_DEFAULT,", wrapper)
        self.assertIn("STARTUP_PROXY_TWO_KI_DEFAULT};", wrapper)

        self.assertIn("ahrs_.reset(-1.0f, -1.0f);", vertical)
        self.assertIn("ahrs_.init(two_kp_, two_ki_);", vertical)
        self.assertIn("seed_from_acc_(acc / acc_norm);", vertical)
        self.assertIn("-acc.x(), -acc.y(), -acc.z(),", vertical)

        for token in (
            "halfex = (ay * halfvz - az * halfvy);",
            "halfey = (az * halfvx - ax * halfvz);",
            "halfez = (ax * halfvy - ay * halfvx);",
            "integralFBx += twoKi * halfex * delta_t_sec;",
            "gx += integralFBx;",
            "gx += twoKp * halfex;",
        ):
            self.assertIn(token, mahony)

        self.assertIn(r"k_P=0.1\ \mathrm{s}^{-1}", proof)
        self.assertIn(r"k_I=0.01\ \mathrm{s}^{-2}", proof)
        self.assertIn(
            r"\widetilde{\vct\beta}:=\vct b_g+\vct\beta",
            proof,
        )

    def test_source_seed_and_physical_direction_bound_are_explicit(self):
        proof = _read(DOC / "w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\label{eq:semiglobal-proxy-force-bound}",
            r"\label{eq:semiglobal-proxy-direction-error}",
            r"\label{eq:semiglobal-proxy-source-seed}",
            r"\label{eq:semiglobal-proxy-disturbance-bounds}",
            r"\label{eq:semiglobal-proxy-source-disturbance}",
            r"\overline r_h",
            r"\overline r_b",
        ):
            self.assertIn(marker, proof)
        self.assertIn("first valid accelerometer sample", proof)
        self.assertIn("integral state is zero at reset", proof)

    def test_quantitative_chart_keeps_almost_global_limitation_honest(self):
        proof = _read(DOC / "w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\mathcal N_P:=W^s(\mathcal U_P)",
            r"\label{eq:semiglobal-proxy-chart}",
            r"\theta_\star:=\frac{\pi}{2}-\epsilon",
            r"\label{eq:semiglobal-proxy-chart-invariance}",
            "antipodal",
            "almost-global",
            "not a replacement for the topology",
        ):
            self.assertIn(marker, proof)
        self.assertIn(r"\cite{Mahony2008NonlinearComplementary}", proof)

    def test_quality_gate_has_calculable_entry_time(self):
        proof = _read(DOC / "w3d-semiglobal-stability.tex-part")
        for marker in (
            r"\label{eq:semiglobal-quality-target}",
            r"\label{eq:semiglobal-quality-entry-time}",
            r"V_{P,\infty}<V_Q",
            r"\sin^{-1}(0.075)",
            r"\SI{2}{s}",
            r"\SI{8}{s}",
        ):
            self.assertIn(marker, proof)

    def test_timeout_envelope_is_analytical_and_source_exact(self):
        proof = _read(DOC / "w3d-semiglobal-stability.tex-part")
        wrapper = _flat(
            _read(SRC / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h")
        )

        for marker in (
            r"\label{eq:semiglobal-timeout-source-time}",
            r"\label{eq:semiglobal-timeout-150}",
            r"\label{eq:semiglobal-timeout-energy}",
            r"\label{eq:semiglobal-timeout-tilt}",
            r"\label{eq:semiglobal-timeout-proxy-bias}",
            r"\label{eq:semiglobal-timeout-handoff-set}",
            r"\label{eq:semiglobal-timeout-radius}",
            r"\label{eq:semiglobal-timeout-capture-entry}",
            r"\label{thm:semiglobal-timeout-live}",
            r"\vct R_{H,T}^{\max}\prec\vct R_C",
        ):
            self.assertIn(marker, proof)

        self.assertIn("float proxy_startup_timeout_sec = 150.0f;", wrapper)
        self.assertIn("float proxy_mag_settle_sec = 0.0f;", wrapper)
        self.assertIn("float mag_min_window_sec = 15.0f;", wrapper)
        self.assertIn("float mag_tilt_fallback_sec = 30.0f;", wrapper)
        self.assertIn(
            "const float timeout_sec = std::max(cfg_.proxy_startup_timeout_sec, mag_acquire_deadline);",
            wrapper,
        )
        self.assertIn(
            "const bool ready_by_timeout = proxy_ready && (t_ >= timeout_sec) && mag_gravity_aligned_branch_;",
            wrapper,
        )

    def test_timeout_does_not_invent_heading_or_copy_proxy_bias(self):
        proof = _read(DOC / "w3d-semiglobal-stability.tex-part")
        wrapper = _flat(
            _read(SRC / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h")
        )
        self.assertIn("ready_by_timeout", wrapper)
        self.assertIn("const bool north_ready = !cfg_.with_mag || mag_ref_set_;", wrapper)
        self.assertIn("const bool have_yaw_gauge = std::isfinite(pending_yaw_abs_rad_);", wrapper)
        self.assertIn("impl_.goLive(q_seed,", wrapper)
        self.assertIn("/*allow_acc_bias=*/false", wrapper)

        self.assertIn("does not require", proof)
        self.assertIn(r"\texttt{north\_ready}", proof)
        self.assertIn("Ungauged timeout branch", proof)
        self.assertIn(r"\ref{thm:iss-nomag-quotient}", proof)
        self.assertIn("copies the proxy quaternion, not its private integral state", proof)


if __name__ == "__main__":
    unittest.main()
