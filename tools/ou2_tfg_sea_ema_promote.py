#!/usr/bin/env python3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OU2 = ROOT / "src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
TFG = ROOT / "src/kalman_tfg/SeaStateFusionFilter_TFG.h"
OU2_TEST = ROOT / "tests/kalman_ou_ii/tuner_schedule-test.cpp"
TFG_TEST = ROOT / "tests/kalman_tfg/tfg_orchestrator-test.cpp"


def replace_once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise RuntimeError(f"{label}: expected one anchor, found {n}")
    return text.replace(old, new, 1)


def patch_ou2() -> None:
    text = OU2.read_text()

    text = replace_once(
        text,
        "Adaptive update: exponential smoothing toward targets over ADAPT_TAU_SEC",
        "Adaptive update: tau/sigma smoothing over 0.40*T_sea (fixed seconds retained for ablation)",
        "OU-II header summary",
    )

    text = replace_once(
        text,
        """constexpr float ADAPT_TAU_SEC              = 1.8f;\nconstexpr float ADAPT_EVERY_SECS           = 0.1f;\n""",
        """// Legacy fixed-second horizon remains for explicit ablations.  The deployed\n// common tau/sigma EMA follows measured sea time T_sea=T_z/2, matching OU-III.\nconstexpr float ADAPT_TAU_SEC                  = 1.8f;\nconstexpr float ADAPT_TAU_SEA_PERIODS          = 0.40f;\nconstexpr float TUNER_FREQ_EMA_SEA_PERIODS     = 0.10f;\nconstexpr float ADAPT_EVERY_SECS               = 0.1f;\n""",
        "OU-II adaptation constants",
    )

    text = replace_once(
        text,
        """        freq_input_lpf_.setCutoff(max_freq_hz_);\n        freq_stillness_.setTargetFreqHz(min_freq_hz_);\n        startup_stage_   = StartupStage::Cold;\n""",
        """        freq_input_lpf_.setCutoff(max_freq_hz_);\n        freq_stillness_.setTargetFreqHz(min_freq_hz_);\n        // The frequency EMA belongs to the measurement-only tuner, so the same\n        // sea-time normalization measured for OU-III applies here unchanged.\n        tuner_.setFrequencySmoothingSeaPeriods(TUNER_FREQ_EMA_SEA_PERIODS);\n        startup_stage_   = StartupStage::Cold;\n""",
        "OU-II constructor tuner default",
    )

    text = replace_once(
        text,
        """    void setAdaptationTimeConstants(float tau_sec) {\n        if (std::isfinite(tau_sec) && tau_sec > 0.0f) adapt_tau_sec_ = tau_sec;\n    }\n\n""",
        """    // Fixed-seconds compatibility mode. Selecting it explicitly disables\n    // sea-period scaling for the common tau/sigma EMA.\n    void setAdaptationTimeConstants(float tau_sec) {\n        if (std::isfinite(tau_sec) && tau_sec > 0.0f) {\n            adapt_tau_sec_ = tau_sec;\n            adapt_tau_sea_periods_ = 0.0f;\n        }\n    }\n\n    void setAdaptationSeaPeriods(float periods) {\n        if (std::isfinite(periods) && periods > 0.0f) {\n            adapt_tau_sea_periods_ = periods;\n        }\n    }\n\n    void setTunerFreqSmoothingSeaPeriods(float periods) {\n        if (std::isfinite(periods) && periods > 0.0f) {\n            tuner_freq_ema_sea_periods_ = periods;\n            tuner_.setFrequencySmoothingSeaPeriods(periods);\n        }\n    }\n\n    void setTunerFreqSmoothingTimeConstant(float tau_sec) {\n        if (std::isfinite(tau_sec) && tau_sec > 0.0f) {\n            tuner_freq_ema_sea_periods_ = 0.0f;\n            tuner_.setTauFreq(tau_sec);\n        }\n    }\n\n    float getAdaptationSeaPeriods() const noexcept { return adapt_tau_sea_periods_; }\n    float getTunerFreqSmoothingSeaPeriods() const noexcept {\n        return tuner_.getFrequencySmoothingSeaPeriods();\n    }\n\n""",
        "OU-II adaptation API",
    )

    text = replace_once(
        text,
        """        adapt_mekf(dt, tau_target_, sigma_target_, R_p0_std_target_, R_v0_std_target_);\n    }\n\n    void adapt_mekf(float dt, float tau_t, float sigma_t, float R_p0_t, float R_v0_t) {\n        const float alpha = 1.0f - std::exp(-dt / adapt_tau_sec_);\n""",
        """        const float sea_time_sec = 0.5f / f_tune;  // T_sea = T_z/2\n        adapt_mekf(dt, tau_target_, sigma_target_, R_p0_std_target_, R_v0_std_target_,\n                   sea_time_sec);\n    }\n\n    void adapt_mekf(float dt, float tau_t, float sigma_t, float R_p0_t, float R_v0_t,\n                    float sea_time_sec) {\n        float adapt_sec = adapt_tau_sec_;\n        if (adapt_tau_sea_periods_ > 0.0f &&\n            std::isfinite(sea_time_sec) && sea_time_sec > 0.0f) {\n            adapt_sec = std::max(dt, adapt_tau_sea_periods_ * sea_time_sec);\n        }\n        const float alpha = 1.0f - std::exp(-dt / adapt_sec);\n""",
        "OU-II sea-time adaptation",
    )

    text = replace_once(
        text,
        """    float adapt_tau_sec_          = ADAPT_TAU_SEC;\n    float adapt_R_p0_mult_        = ADAPT_R_p0_MULT;\n""",
        """    float adapt_tau_sec_              = ADAPT_TAU_SEC;\n    float adapt_tau_sea_periods_      = ADAPT_TAU_SEA_PERIODS;\n    float tuner_freq_ema_sea_periods_ = TUNER_FREQ_EMA_SEA_PERIODS;\n    float adapt_R_p0_mult_            = ADAPT_R_p0_MULT;\n""",
        "OU-II adaptation members",
    )

    OU2.write_text(text)


def patch_tfg() -> None:
    text = TFG.read_text()

    text = replace_once(
        text,
        """constexpr float SIGMA_BAND_MIN_HZ_DEFAULT     = 0.01f;\nconstexpr float SIGMA_BAND_MAX_HZ_DEFAULT     = 6.0f;\n\nstruct TfgTuneState {\n""",
        """constexpr float SIGMA_BAND_MIN_HZ_DEFAULT     = 0.01f;\nconstexpr float SIGMA_BAND_MAX_HZ_DEFAULT     = 6.0f;\n\n// The two EMA horizons that remained fixed in wall-clock seconds now use the\n// same measured sea-time normalization as OU-III. Fixed-second compatibility\n// remains available through the setters below.\nconstexpr float ADAPT_TAU_SEC_DEFAULT              = 1.8f;\nconstexpr float ADAPT_TAU_SEA_PERIODS_DEFAULT      = 0.40f;\nconstexpr float TUNER_FREQ_EMA_SEA_PERIODS_DEFAULT = 0.10f;\n\nstruct TfgTuneState {\n""",
        "TFG sea-EMA constants",
    )

    text = replace_once(
        text,
        """        tuner_ = ::SeaStateAutoTuner(2.0f, 1.0f);\n        wave_period_.reset();\n""",
        """        tuner_ = ::SeaStateAutoTuner(2.0f, tuner_freq_tau_sec_);\n        if (tuner_freq_ema_sea_periods_ > 0.0f) {\n            tuner_.setFrequencySmoothingSeaPeriods(tuner_freq_ema_sea_periods_);\n        }\n        sea_time_sec_ = 0.5f / 0.2f;\n        wave_period_.reset();\n""",
        "TFG begin tuner default",
    )

    text = replace_once(
        text,
        """    void setAdaptationTimeConstants(float tau_sec) {\n        if (tau_sec > 0.0f && std::isfinite(tau_sec)) adapt_tau_sec_ = tau_sec;\n    }\n    void setRSAdaptMult(float m)     { if (m > 0.0f && std::isfinite(m)) adapt_RS_mult_ = m; }\n""",
        """    // Fixed-seconds compatibility mode for the common tau/sigma EMA.\n    void setAdaptationTimeConstants(float tau_sec) {\n        if (tau_sec > 0.0f && std::isfinite(tau_sec)) {\n            adapt_tau_sec_ = tau_sec;\n            adapt_tau_sea_periods_ = 0.0f;\n        }\n    }\n    void setAdaptationSeaPeriods(float periods) {\n        if (periods > 0.0f && std::isfinite(periods)) {\n            adapt_tau_sea_periods_ = periods;\n        }\n    }\n    void setTunerFreqSmoothingSeaPeriods(float periods) {\n        if (periods > 0.0f && std::isfinite(periods)) {\n            tuner_freq_ema_sea_periods_ = periods;\n            tuner_.setFrequencySmoothingSeaPeriods(periods);\n        }\n    }\n    void setTunerFreqSmoothingTimeConstant(float tau_sec) {\n        if (tau_sec > 0.0f && std::isfinite(tau_sec)) {\n            tuner_freq_tau_sec_ = tau_sec;\n            tuner_freq_ema_sea_periods_ = 0.0f;\n            tuner_.setTauFreq(tau_sec);\n        }\n    }\n    [[nodiscard]] float getAdaptationSeaPeriods() const noexcept {\n        return adapt_tau_sea_periods_;\n    }\n    [[nodiscard]] float getTunerFreqSmoothingSeaPeriods() const noexcept {\n        return tuner_.getFrequencySmoothingSeaPeriods();\n    }\n    void setRSAdaptMult(float m)     { if (m > 0.0f && std::isfinite(m)) adapt_RS_mult_ = m; }\n""",
        "TFG adaptation API",
    )

    text = replace_once(
        text,
        """        float f_tune = tuner_.isFreqReady() ? tuner_.getFrequencyHz() : kMinTuneFreqHz;\n        if (!std::isfinite(f_tune) || f_tune < kMinTuneFreqHz) f_tune = kMinTuneFreqHz;\n        f_tune = std::min(f_tune, kMaxFreqHz);\n\n        const float band_noise_sigma = bandNoiseFloorSigma_();\n""",
        """        float f_tune = tuner_.isFreqReady() ? tuner_.getFrequencyHz() : kMinTuneFreqHz;\n        if (!std::isfinite(f_tune) || f_tune < kMinTuneFreqHz) f_tune = kMinTuneFreqHz;\n        f_tune = std::min(f_tune, kMaxFreqHz);\n        sea_time_sec_ = 0.5f / f_tune;  // T_sea = T_z/2\n\n        const float band_noise_sigma = bandNoiseFloorSigma_();\n""",
        "TFG measured sea time",
    )

    text = replace_once(
        text,
        """        // The smoother sees every sample; only the commit is cadenced.\n        if (!freeze_ou_channel_) {\n            const float a = 1.0f - std::exp(-dt / adapt_tau_sec_);\n            tune_.tau_applied   += a * (tau_target_ - tune_.tau_applied);\n            tune_.sigma_applied += a * (sigma_target_ - tune_.sigma_applied);\n        }\n""",
        """        // The smoother sees every sample; only the commit is cadenced.\n        if (!freeze_ou_channel_) {\n            float adapt_sec = adapt_tau_sec_;\n            if (adapt_tau_sea_periods_ > 0.0f &&\n                std::isfinite(sea_time_sec_) && sea_time_sec_ > 0.0f) {\n                adapt_sec = std::max(dt, adapt_tau_sea_periods_ * sea_time_sec_);\n            }\n            const float a = 1.0f - std::exp(-dt / adapt_sec);\n            tune_.tau_applied   += a * (tau_target_ - tune_.tau_applied);\n            tune_.sigma_applied += a * (sigma_target_ - tune_.sigma_applied);\n        }\n""",
        "TFG sea-time adaptation",
    )

    text = replace_once(
        text,
        """    float adapt_tau_sec_    = 1.8f;\n    float adapt_RS_mult_    = 3.0f;\n    float adapt_RS_slew_log_ = 0.0f;\n""",
        """    float adapt_tau_sec_              = ADAPT_TAU_SEC_DEFAULT;\n    float adapt_tau_sea_periods_      = ADAPT_TAU_SEA_PERIODS_DEFAULT;\n    float tuner_freq_tau_sec_          = 1.0f;\n    float tuner_freq_ema_sea_periods_ = TUNER_FREQ_EMA_SEA_PERIODS_DEFAULT;\n    float sea_time_sec_                = 0.5f / 0.2f;\n    float adapt_RS_mult_                = 3.0f;\n    float adapt_RS_slew_log_            = 0.0f;\n""",
        "TFG adaptation members",
    )

    TFG.write_text(text)


def patch_ou2_test() -> None:
    text = OU2_TEST.read_text()

    text = replace_once(
        text,
        """    if (!near(f.getPseudoUpdateTauRatio(), 0.015f / 1.1f) ||\n        !near(f.getPseudoUpdatePeriodSec(), 0.015f)) {\n        std::cerr << \"FAIL: OU-II tau-scaled pseudo cadence did not preserve the nominal point\\n\";\n        return 1;\n    }\n\n""",
        """    if (!near(f.getPseudoUpdateTauRatio(), 0.015f / 1.1f) ||\n        !near(f.getPseudoUpdatePeriodSec(), 0.015f)) {\n        std::cerr << \"FAIL: OU-II tau-scaled pseudo cadence did not preserve the nominal point\\n\";\n        return 1;\n    }\n    if (!near(f.getAdaptationSeaPeriods(), 0.40f) ||\n        !near(f.getTunerFreqSmoothingSeaPeriods(), 0.10f)) {\n        std::cerr << \"FAIL: OU-II sea-scaled EMA defaults are not 0.40/0.10 of T_sea\\n\";\n        return 1;\n    }\n\n    // A sample inside the activation interval must still move the EMA state.\n    {\n        Filter s(false);\n        s.time_ = 2.0;\n        s.last_adapt_time_sec_ = 2.0;\n        s.online_tune_apply_pending_ = false;\n        const float before = s.tune_.tau_applied;\n        constexpr float T_sea = 2.5f;\n        constexpr float dt_ema = 0.1f;\n        const float expected_alpha = 1.0f - std::exp(-dt_ema / (0.40f * T_sea));\n        const float expected = before + expected_alpha * (3.0f - before);\n        s.adapt_mekf(dt_ema, 3.0f, 0.8f, 1.2f, 0.7f, T_sea);\n        if (!near(s.tune_.tau_applied, expected, 1e-5f)) {\n            std::cerr << \"FAIL: OU-II sea-scaled EMA did not consume the physical sample\\n\";\n            return 1;\n        }\n        if (s.online_tune_apply_pending_) {\n            std::cerr << \"FAIL: OU-II EMA sample incorrectly bypassed activation cadence\\n\";\n            return 1;\n        }\n    }\n\n""",
        "OU-II default/cadence test",
    )

    text = replace_once(
        text,
        "f.adapt_mekf(0.1f, 3.0f, 0.8f, 1.2f, 0.7f);",
        "f.adapt_mekf(0.1f, 3.0f, 0.8f, 1.2f, 0.7f, 2.5f);",
        "OU-II adapt_mekf test call",
    )

    OU2_TEST.write_text(text)


def patch_tfg_test() -> None:
    text = TFG_TEST.read_text()

    anchor = """void test_schedule_is_exogenous_to_the_current_sample() {\n"""
    new_test = """void test_sea_scaled_ema_defaults() {\n    Fusion f;\n    f.begin(default_config());\n    check(std::fabs(f.getAdaptationSeaPeriods() - 0.40f) <= 1e-6f,\n          \"TFG common tau/sigma EMA is not 0.40*T_sea by default\");\n    check(std::fabs(f.getTunerFreqSmoothingSeaPeriods() - 0.10f) <= 1e-6f,\n          \"TFG tuner frequency EMA is not 0.10*T_sea by default\");\n\n    // Fixed-second setters remain real ablations rather than being ignored by\n    // the new defaults.\n    f.setAdaptationTimeConstants(1.8f);\n    f.setTunerFreqSmoothingTimeConstant(1.0f);\n    check(f.getAdaptationSeaPeriods() == 0.0f,\n          \"TFG fixed tau/sigma EMA setter did not disable sea scaling\");\n    check(f.getTunerFreqSmoothingSeaPeriods() == 0.0f,\n          \"TFG fixed frequency EMA setter did not disable sea scaling\");\n}\n\n""" + anchor
    text = replace_once(text, anchor, new_test, "TFG sea-EMA test function")

    text = replace_once(
        text,
        """    test_adaptation_is_cadenced();\n    test_schedule_is_exogenous_to_the_current_sample();\n""",
        """    test_adaptation_is_cadenced();\n    test_sea_scaled_ema_defaults();\n    test_schedule_is_exogenous_to_the_current_sample();\n""",
        "TFG main test call",
    )

    TFG_TEST.write_text(text)


def main() -> None:
    patch_ou2()
    patch_tfg()
    patch_ou2_test()
    patch_tfg_test()
    print("PROMOTED OU-II and TFG sea-scaled EMA defaults (0.40, 0.10 of T_sea)")


if __name__ == "__main__":
    main()
