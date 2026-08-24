#!/usr/bin/env python3
"""Pin deployed OU-II/OU-III paper claims to the implementation.

Numerical evidence has separate provenance tests.  This source-level publication
contract protects state/model semantics and, deliberately, the complete adaptive
smoothing chain: estimator EWMAs, scheduler EWMAs, clamps, activation/hold
semantics, startup filtering, and magnetic exponential memory/slew.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]


def text(path: str) -> str:
    return (ROOT / path).read_text(encoding="utf-8")


class OUPaperCodeParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ou3_core = text("src/kalman_ou_iii/Kalman3D_Wave_OU_III.h")
        cls.ou3_wrap = text("src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h")
        cls.ou2_core = text("src/kalman_ou_ii/Kalman3D_Wave_OU_II.h")
        cls.ou2_wrap = text("src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h")
        cls.tuner = text("src/tuner/SeaStateAutoTuner.h")
        cls.period = text("src/tuner/WavePeriodEstimator.h")
        cls.limits = text("src/tuner/SeaStateAdaptationLimits.h")
        cls.mag_hi = text("src/tuner/ContinuousMagHardIronEstimator.h")

        cls.ou3_state = text("doc/kalman_ou_iii/w3d-state.tex-part")
        cls.ou3_proc = text("doc/kalman_ou_iii/w3d-proc-cont.tex-part")
        cls.ou3_meas = text("doc/kalman_ou_iii/w3d-meas.tex-part")
        cls.ou3_adapt = text("doc/kalman_ou_iii/w3d-adaptation-motivation.tex-part")
        cls.ou3_obs = text("doc/kalman_ou_iii/w3d-adaptation-observables.tex-part")
        cls.ou3_impl = text("doc/kalman_ou_iii/w3d-fus-methods.tex-part")
        cls.ou3_iss = text("doc/kalman_ou_iii/w3d-iss-stability.tex-part")
        cls.ou3_init = text("doc/kalman_ou_iii/w3d-init.tex-part")
        cls.ou3_mag = text("doc/kalman_ou_iii/w3d-mag-hard-iron.tex-part")
        cls.ou2_paper = text("doc/kalman_ou_ii/ou2-dual-regularization-mse.tex")
        cls.ou2_iss = text("doc/kalman_ou_ii/ou2-iss-stability.tex-part")
        cls.ou2_charts = text("doc/kalman_ou_ii/kalman_ou_ii-charts.tex")

    def test_ou3_state_and_process_match_21_state_core(self):
        self.assertIn("static constexpr int EXT_ADD = 12", self.ou3_core)
        self.assertIn("static constexpr int OFF_S   = BASE_N + 6", self.ou3_core)
        self.assertIn("static constexpr int OFF_AW  = BASE_N + 9", self.ou3_core)
        self.assertIn("BASE_N + 12", self.ou3_core)

        self.assertIn("N_X=21", self.ou3_state)
        for token in (
            r"\delta\vct{\theta}", r"\delta\vct{b}_g", r"\delta\vct{v}",
            r"\delta\vct{p}", r"\delta\vct{S}", r"\delta\vct{a}_w",
            r"\delta\vct{b}_a",
        ):
            self.assertIn(token, self.ou3_state)
        self.assertIn(r"\dot{\vct{v}}=\vct{a}_w", self.ou3_proc)
        self.assertIn(r"\dot{\vct{p}}=\vct{v}", self.ou3_proc)
        self.assertIn(r"\dot{\vct{S}}=\vct{p}", self.ou3_proc)
        self.assertIn(r"\vct{z}_S=0=\vct{S}+\vct{n}_S", self.ou3_meas)

    def test_ou3_spectral_mse_default_matches_paper(self):
        self.assertIn(
            "RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;",
            self.ou3_wrap,
        )
        self.assertRegex(self.ou3_wrap, r"R_S_MSE_COEFF_DEFAULT\s*=\s*0\.0538f")
        self.assertIn("0.0148f * 0.0148f * FREQ_SMOOTHER_DT", self.ou3_wrap)
        self.assertIn("std::pow(u, 6.0f / 7.0f)", self.ou3_wrap)
        self.assertIn("/ std::sqrt(TS)", self.ou3_wrap)
        self.assertIn("if (rs_law_ != RSAdaptationLaw::Cubic) return 1.0f;", self.ou3_wrap)

        self.assertIn(r"C_J=0.0538", self.ou3_impl)
        self.assertIn(r"\widehat\sigma_{a,B,\star}^{\,6/7}", self.ou3_impl)
        self.assertIn(r"\tau_\star^{24/7}T_{S,\star}^{-1/2}", self.ou3_impl)
        self.assertIn(r"\sqrt{R_a}=\SI{0.0148}{m.s^{-2}}", self.ou3_adapt)
        self.assertIn("No additional reference-cadence factor is applied", self.ou3_impl)

    def test_common_variance_ewma_is_single_stage_debiased_and_documented(self):
        self.assertIn("struct DebiasedEMA", self.tuner)
        self.assertIn("value  = (1.0f - alpha) * value + alpha * x;", self.tuner)
        self.assertIn("weight = (1.0f - alpha) * weight + alpha;", self.tuner)
        self.assertIn("A_mean.update(accel, alpha_var);", self.tuner)
        self.assertIn("A_sq.update(accel * accel, alpha_var);", self.tuner)
        self.assertIn("explicit SeaStateAutoTuner(float K_periods_ = 4.0f,", self.tuner)
        self.assertIn("K_periods * T_eff", self.tuner)
        self.assertIn("float tau_var_min_sec = 0.3f;", self.tuner)
        self.assertIn("float tau_var_max_sec = 60.0f;", self.tuner)
        self.assertIn("There is no second EWMA of the already-smoothed variance", self.ou3_obs)
        self.assertIn(r"$K_{\mathrm{periods}}=4$", self.ou3_obs)
        self.assertIn(r"\SI{30}{s} guard therefore binds on the largest", self.ou3_obs)

    def test_wave_period_ewmas_match_log_period_source_exactly(self):
        for token in (
            "float moment_horizon_periods = 4.0f",
            "float log_smoothing_periods = 0.05f",
            "float min_horizon_sec = 20.0f",
            "float max_horizon_sec = 180.0f",
            "const float settle_sec = 6.0f / lambda_;",
            "return (std::isfinite(period) && period > 0.0f) ? period : 6.0f;",
            "log_period_sec_ = log_raw;",
            "const float requested = log_smoothing_periods_ * sea_period;",
            "log_period_sec_ += alpha * (log_raw - log_period_sec_);",
            "return std::isfinite(log_period_sec_) ? std::exp(-log_period_sec_) : NAN;",
        ):
            self.assertIn(token, self.period)
        self.assertIn(r"\tau_{\mathrm{mom}}", self.ou3_obs)
        self.assertIn(r"\clip\!\left(4T_{z,\mathrm{can}},\SI{20}{s},\SI{180}{s}\right)", self.ou3_obs)
        self.assertIn("The canonical period is not an arithmetic EMA of frequency", self.ou3_obs)
        self.assertIn(r"0.05e^{\ell_{k-1}}", self.ou3_obs)
        self.assertIn(r"f_{\mathrm{tune}}=e^{-\ell}", self.ou3_obs)
        self.assertNotIn("bounded EMA of $1/\\widehat T_z$", self.ou3_obs)

    def test_dynamic_ema_guards_match_paper(self):
        self.assertIn("kDynamicEmaTimeScaleMinSec = 0.5f", self.limits)
        self.assertIn("kDynamicEmaTimeScaleMaxSec = 6.0f", self.limits)
        self.assertIn("kDynamicEmaHorizonMinSec = 0.05f", self.limits)
        self.assertIn("kDynamicEmaHorizonMaxSec = 30.0f", self.limits)
        self.assertIn("$[\\SI{0.5}{s},\\SI{6}{s}]$", self.ou3_impl)
        self.assertIn("$[\\SI{0.05}{s},\\SI{30}{s}]$", self.ou3_impl)

    def test_ou3_candidate_emas_and_activation_hold_match_source(self):
        for token in (
            "ADAPT_TAU_SEA_PERIODS          = 0.40f",
            "ADAPT_RS_MULT              = 1.5f",
            "ADAPT_RS_SLEW_LOG          = 0.0f",
            "ADAPT_EVERY_SECS           = 0.1f",
            "if (time_ - last_adapt_time_sec_ > adapt_every_secs_)",
            "online_tune_apply_pending_ = true;",
            "apply_pending_online_tune_();",
        ):
            self.assertIn(token, self.ou3_wrap)
        self.assertIn("candidate tuple", self.ou3_impl)
        self.assertIn("sample-and-hold schedule", self.ou3_impl)
        self.assertIn("parameter activation cadence & approximately $0.1$ s", self.ou3_impl)
        self.assertIn("$1.5\\tau_\\star$", self.ou3_impl)
        self.assertIn("discrepancy-based\nhorizon-shortening term is disabled", self.ou3_impl)

        # The stability source family must model the same two-layer recurrence,
        # not silently apply the candidate EMA directly to the active MEKF.
        self.assertIn(r"\widetilde x_{k+1}", self.ou3_iss)
        self.assertIn(r"\label{eq:iss-source-hold-commit}", self.ou3_iss)
        self.assertIn("active schedule is sample-and-hold", self.ou3_iss)
        self.assertIn("activation timer", self.ou3_iss)
        self.assertIn("on non-commit samples it is exactly\nzero", self.ou3_iss)

    def test_ou3_startup_gravity_ema_matches_source(self):
        self.assertRegex(
            self.ou3_wrap,
            r"float\s+mag_gravity_align_world_tau_sec\s*=\s*12\.0f;",
        )
        self.assertRegex(
            self.ou3_wrap,
            r"float\s+mag_gravity_align_world_warmup_sec\s*=\s*5\.0f;",
        )
        self.assertIn("const float alpha = 1.0f - std::exp(-dt / tau);", self.ou3_wrap)
        self.assertIn("state = x;\n                initialized = true;", self.ou3_wrap)
        self.assertIn("gravity_gate_acc_world_lpf_.reset();", self.ou3_wrap)
        self.assertIn(r"\alpha_{g,k}=1-e^{-\Delta t_k/\SI{12}{s}}", self.ou3_init)
        self.assertIn("first valid sample initializes the EMA state directly", self.ou3_init)
        self.assertIn(r"ignored until this world-frame average has run for\n\SI{5}{s}", self.ou3_init)

    def test_ou3_cadence_bias_and_outer_warmup_match_paper(self):
        for token in (
            "PSEUDO_UPDATE_PERIOD_NOMINAL_S = 0.015f",
            "PSEUDO_UPDATE_TAU_NOMINAL_S = 1.1f",
            "PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT = FREQ_SMOOTHER_DT",
            "PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT = 0.25f",
        ):
            self.assertIn(token, self.ou3_wrap)
        self.assertIn(r"c_T=\frac{\SI{15}{ms}}{\SI{1.1}{s}}", self.ou3_impl)
        self.assertIn(r"T_{S,\min}=\SI{5}{ms}", self.ou3_impl)
        self.assertIn(r"T_{S,\max}=\SI{250}{ms}", self.ou3_impl)
        self.assertIn(r"\tau_b=\SI{5000}{s}", self.ou3_impl)
        self.assertIn("float online_tune_warmup_sec = 10.0f;", self.ou3_wrap)
        self.assertIn("impl_.setOnlineTuneWarmupSec(cfg_.online_tune_warmup_sec);", self.ou3_wrap)
        self.assertIn("online tuning warmup / magnetometer delay & $10/7$ s", self.ou3_impl)

    def test_ou3_magnetic_exponential_memory_and_slew_match_paper(self):
        self.assertIn("float memory_sec = 600.0f;", self.mag_hi)
        self.assertIn("std::exp(-double(dt) / double(cfg_.memory_sec))", self.mag_hi)
        self.assertIn("float solve_period_sec = 1.0f;", self.mag_hi)
        self.assertIn("float mag_hi_slew_tau_sec             = 45.0f;", self.ou3_wrap)
        self.assertIn("1.0f - std::exp(-dt_apply / tau)", self.ou3_wrap)
        self.assertIn(r"\SI{600}{s} memory", self.ou3_mag)
        self.assertIn(r"\SI{45}{s} time", self.ou3_mag)

    def test_ou3_magnetic_refinement_uses_independent_proxy_tilt(self):
        self.assertIn("const Eigen::Quaternionf q_tilt_bw = impl_.startupProxyTiltQuat();", self.ou3_wrap)
        self.assertIn("Magnetic acquisition uses proxy tilt rather than the MEKF", self.ou3_init)
        self.assertIn("later\nindependent acquisition refines that reference", self.ou3_init)

    def test_ou2_state_and_dual_pseudos_match_core(self):
        self.assertIn("static constexpr int EXT_ADD = 9", self.ou2_core)
        self.assertIn("static constexpr int OFF_P   = BASE_N + 3", self.ou2_core)
        self.assertIn("static constexpr int OFF_AW  = BASE_N + 6", self.ou2_core)
        self.assertIn("measurement_update_position_pseudo", self.ou2_core)
        self.assertIn("measurement_update_velocity_pseudo", self.ou2_core)
        self.assertIn(
            "local error coordinates containing attitude, gyroscope bias, velocity,\n"
            "displacement, latent OU acceleration, and accelerometer bias",
            self.ou2_paper,
        )
        self.assertIn("one on world-frame\ndisplacement $p$ and one on world-frame velocity $v$", self.ou2_paper)
        self.assertIn("both applied with period $T_S$", self.ou2_paper)

    def test_ou2_physical_mse_default_matches_paper(self):
        self.assertIn("PseudoAdaptationLaw pseudo_law_ = PseudoAdaptationLaw::PhysicalMSE;", self.ou2_wrap)
        self.assertRegex(self.ou2_wrap, r"R_PSEUDO_MSE_COEFF_DEFAULT\s*=\s*0\.1116f")
        self.assertRegex(self.ou2_wrap, r"R_PSEUDO_MSE_RATIO_DEFAULT\s*=\s*0\.4611f")
        self.assertIn("ACC_NOISE_FLOOR_SIGMA_DEFAULT * ACC_NOISE_FLOOR_SIGMA_DEFAULT * FREQ_SMOOTHER_DT", self.ou2_wrap)
        self.assertIn("std::pow(u, 0.8f) / std::sqrt(TS)", self.ou2_wrap)
        self.assertIn("r_v = r_p / (ratio * tau);", self.ou2_wrap)
        self.assertIn("if (pseudo_law_ != PseudoAdaptationLaw::Empirical) return 1.0f;", self.ou2_wrap)
        for token in (
            r"C_P=0.1116", r"\frac{C_P}{C_V}=0.4611",
            r"\sigma_{\rm floor}=\SI{0.12}{m.s^{-2}}", r"h=\SI{0.005}{s}",
            r"\tau^{12/5}T_S^{-1/2}", r"\tau^{7/5}T_S^{-1/2}",
        ):
            self.assertIn(token, self.ou2_paper)
        self.assertIn("no additional historical\ncadence renormalization is applied", self.ou2_paper)

    def test_ou2_candidate_emas_and_activation_hold_match_source(self):
        for token in (
            "ADAPT_TAU_SEA_PERIODS          = 0.40f",
            "ADAPT_R_p0_MULT            = 3.0f",
            "ADAPT_R_v0_MULT            = 3.0f",
            "ADAPT_R_SLEW_LOG           = 0.0f",
            "ADAPT_EVERY_SECS               = 0.1f",
            "if (time_ - last_adapt_time_sec_ > adapt_every_secs_)",
            "online_tune_apply_pending_ = true;",
            "apply_pending_online_tune_();",
        ):
            self.assertIn(token, self.ou2_wrap)
        self.assertIn("sample-and-hold", self.ou2_iss)
        self.assertIn("approximately \\SI{0.1}{s}", self.ou2_iss)
        self.assertIn("activation timer", self.ou2_charts)
        self.assertNotIn("candidate becomes active\nonly at step $k+1$", self.ou2_charts)

        # OU-II uses the same tau-scaled pseudo cadence as OU-III, not a fixed
        # 15 ms schedule.  The paper must describe the source scheduler.
        self.assertIn("PSEUDO_UPDATE_TAU_RATIO_DEFAULT", self.ou2_wrap)
        self.assertIn("PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT = FREQ_SMOOTHER_DT", self.ou2_wrap)
        self.assertIn("PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT = 0.25f", self.ou2_wrap)
        self.assertIn(r"\label{eq:ou2-iss-pseudo-cadence-source}", self.ou2_iss)
        self.assertIn(r"\SI{5}{ms},\SI{250}{ms}", self.ou2_iss)
        self.assertNotIn("OU--II retains the fixed pseudo-update period", self.ou2_iss)

    def test_ou2_residual_noise_parameter_is_not_a_live_tuner_alias(self):
        self.assertIn("float acc_noise_floor_sigma_ = ACC_NOISE_FLOOR_SIGMA_DEFAULT;", self.ou2_wrap)
        self.assertIn("float pseudo_accel_noise_density_ = R_PSEUDO_ACCEL_NOISE_DENSITY_DEFAULT;", self.ou2_wrap)
        self.assertIn("setPseudoAccelNoiseDensity", self.ou2_wrap)
        self.assertIn("stored as a separate scheduler parameter", self.ou2_paper)
        self.assertIn("does not live-read the tuner's noise-floor setting", self.ou2_paper)


if __name__ == "__main__":
    unittest.main()
