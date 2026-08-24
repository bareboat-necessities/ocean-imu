#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    Sea-state orchestration for the right-invariant two-frame Lie-group filter.

    The translational plant remains the integrated OU chain [v,p,S,a_w], but
    the estimator geometry is TFG-specific.  The measurement-only front end is
    intentionally kept at coefficient parity with deployed OU-III: canonical
    log-period statistics, period-scaled acceleration variance, self-scaled
    parameter EMAs, startup tilt/magnetic acquisition and continuous hard iron.

    r_S defaults to the same reduced spectral-MSE law as current OU-III,

        r_S = C_J q_eff^(1/14) sigma_a,B^(6/7) tau^(24/7) / sqrt(T_S),

    with T_S proportional to tau away from cadence clamps, giving the effective
    tau^(41/14) dependence.  The previously shipped TFG law

        r_S,base = C_R sigma_aw tau^3,
        r_S,filter = r_S,base sqrt(T_0/T_S)

    remains available as LegacyCubic for embedded targets and exact historical
    comparisons.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdint>

#include "kalman_common/SeaStateFusionFilterCommon.h"
#include "kalman_tfg/Kalman3D_Wave_TFG.h"
#include "tuner/AdaptiveWaveBandPass.h"
#include "tuner/ContinuousMagHardIronEstimator.h"
#include "tuner/MagAutoTuner.h"
#include "tuner/SeaStateAutoTuner.h"
#include "tuner/VerticalAccelComplementary.h"
#include "tuner/WavePeriodEstimator.h"

namespace ocean_imu::tfg {

constexpr float MAG_DELAY_SEC = 7.0f;
constexpr float SIGMA_BAND_LOW_RATIO_DEFAULT  = 0.5f;
constexpr float SIGMA_BAND_HIGH_RATIO_DEFAULT = 4.0f;
constexpr float SIGMA_BAND_MIN_HZ_DEFAULT     = 0.01f;
constexpr float SIGMA_BAND_MAX_HZ_DEFAULT     = 6.0f;
constexpr float ADAPT_TAU_SEC_DEFAULT          = 1.8f;
constexpr float ADAPT_TAU_SEA_PERIODS_DEFAULT  = 0.40f;
constexpr float ADAPT_RS_MULT_DEFAULT          = 1.5f;
constexpr float TUNER_SIGMA_VAR_K_PERIODS_DEFAULT = 4.0f;

// Same strong-observation sensor point and analytical spectral coefficient as
// deployed OU-III.  The physical wave RMS is recovered from sigma_aw below, so
// TFG's independently fitted sigma_coeff does not alter the distortion cost.
constexpr float TFG_NOMINAL_DT = 1.0f / 200.0f;
constexpr float R_S_ACCEL_NOISE_DENSITY_DEFAULT =
    0.0148f * 0.0148f * TFG_NOMINAL_DT;
constexpr float R_S_MSE_COEFF_DEFAULT = 0.0538f;

enum class RSAdaptationLaw : uint8_t {
    LegacyCubic = 0,
    SpectralMSE = 1,
};

struct TfgTuneState {
    float tau_applied   = 1.1f;
    float sigma_applied = 1e-2f;
    float RS_applied    = 0.5f;
};

template <typename MekfT = ocean_imu::kalman::Kalman3D_Wave_TFG<float>>
class SeaStateFusionFilter_TFG {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Mekf = MekfT;
    using Vector3f = Eigen::Vector3f;
    enum class StartupStage { Cold, TunerWarm, Live };
    enum class StartupInitPolicy { StagedMekf, MahonyProxy };
    using RSLaw = RSAdaptationLaw;

    struct Config {
        Vector3f sigma_a{Vector3f::Constant(0.5f)};
        float    gyro_noise_density = 0.005f;
        Vector3f sigma_m{Vector3f::Constant(0.1f)};
        float    gyro_bias_rw_var   = 1e-10f;
        float    initial_covariance = 1e-4f;
        bool     with_mag = true;
        float    mag_delay_sec = MAG_DELAY_SEC;
        bool     freeze_acc_bias_until_live = true;
        float    Racc_warmup_std = 0.5f;
        float    gravity_magnitude = 9.80665f;

        // Deployed OU-III wrapper value.  The underlying low-level OU class has
        // a 5 s compatibility default, but the deployed front end sets 10 s.
        float    online_tune_warmup_sec = 10.0f;
        StartupInitPolicy startup_init_policy = StartupInitPolicy::MahonyProxy;

        float proxy_startup_min_sec     = 8.0f;
        float proxy_startup_timeout_sec = 150.0f;
        float proxy_handoff_tilt_sigma_rad = 0.035f;
        float proxy_handoff_yaw_sigma_rad  = 0.087f;
        float proxy_handoff_yaw_sigma_free_rad = 1.5708f;

        // Same private Mahony gains as deployed OU-III.
        float proxy_two_kp = 0.2f;
        float proxy_two_ki = 0.02f;

        // Compatibility field retained for existing TFG studies.  At the
        // reference 25 Hz mag ODR, 10 s corresponds to OU-III's 250 updates.
        float acc_bias_unlock_sec = 10.0f;
        float handoff_acc_bias_std = 0.03f;

        // The actual deployed OU-III magnetic gravity gate is in WORLD frame.
        // Keep the historical TFG field names so source compatibility is not
        // broken, but give them the deployed OU values and interpretation.
        float proxy_gravity_align_sin = 0.075f;
        float proxy_gravity_lpf_sec   = 12.0f;
        float proxy_gravity_hold_sec  = 2.0f;
        float proxy_gravity_warmup_sec = 5.0f;
        float mag_extreme_gyro_dps    = 30.0f;
        float mag_tilt_fallback_sec   = 30.0f;
        float mag_init_min_mag_norm   = 1e-3f;

        int   mag_min_samples    = 128;
        float mag_min_window_sec = 15.0f;
        float mag_max_window_sec = 0.0f;
        float mag_sample_dt_sec  = 1.0f / 200.0f;
        float proxy_mag_settle_sec = 0.0f;

        bool  mag_refine_enabled    = true;
        float mag_refine_start_sec  = 90.0f;
        float mag_refine_window_sec = 30.0f;

        bool  mag_enable_quality_weighting = false;
        float mag_min_effective_weight     = 0.0f;
        float mag_acc_norm_rel_soft        = 0.22f;
        float mag_gyro_soft_dps            = 45.0f;
        bool  mag_estimate_hard_iron       = false;

        bool  mag_continuous_hard_iron        = true;
        float mag_hi_memory_sec               = 600.0f;
        float mag_hi_model_ridge              = 5.0e-4f;
        float mag_hi_model_ridge_relative     = 0.5f;
        float mag_hi_min_information          = 2.0f;
        float mag_hi_min_effective_weight     = 500.0f;
        float mag_hi_max_residual_rms_uT      = 3.0f;
        float mag_hi_max_bias_fraction        = 0.35f;
        float mag_hi_apply_fraction           = 1.0f;
        float mag_hi_slew_tau_sec             = 45.0f;

        bool  wave_band_tuning      = true;
        float sigma_band_low_ratio  = SIGMA_BAND_LOW_RATIO_DEFAULT;
        float sigma_band_high_ratio = SIGMA_BAND_HIGH_RATIO_DEFAULT;
        float sigma_band_min_hz     = SIGMA_BAND_MIN_HZ_DEFAULT;
        float sigma_band_max_hz     = SIGMA_BAND_MAX_HZ_DEFAULT;
    };

    void begin(const Config& cfg) {
        cfg_ = cfg;
        mekf_ = Mekf(cfg.gyro_noise_density, cfg.gravity_magnitude);
        mekf_.initialize_identity(cfg.initial_covariance);
        mekf_.set_Racc_std(cfg.sigma_a);
        mekf_.set_Rmag_std(cfg.sigma_m);
        mekf_.set_Q_bgyro_rw(Vector3f::Constant(cfg.gyro_bias_rw_var));
        Racc_nominal_ = cfg.sigma_a;

        tuner_ = ::SeaStateAutoTuner(TUNER_SIGMA_VAR_K_PERIODS_DEFAULT);
        sea_time_sec_ = 0.5f / kTuneFreqPriorHz;
        wave_period_.reset();
        vertical_complementary_.setGains(cfg.proxy_two_kp, cfg.proxy_two_ki);
        vertical_complementary_.reset();
        sigma_wave_band_.setRatios(cfg.sigma_band_low_ratio, cfg.sigma_band_high_ratio);
        sigma_wave_band_.setLimitsHz(cfg.sigma_band_min_hz, cfg.sigma_band_max_hz);
        sigma_wave_band_.reset();
        bootstrap_tilt_obs_.reset();
        bootstrap_gravity_slow_lpf_.reset();
        bootstrap_gravity_good_sec_ = 0.0f;
        elapsed_sec_ = 0.0f;
        live_sec_ = 0.0f;
        mag_elapsed_sec_ = 0.0f;
        tuner_warm_sec_ = 0.0f;
        pseudo_elapsed_ = 0.0f;
        adapt_elapsed_sec_ = 0.0f;
        stage_ = StartupStage::Cold;
        beginMagAcquisition_();
        enterCold_();
        commitTune_();
        mekf_.reset_aw_covariance_to_stationary();
    }

    void update(float dt, const Vector3f& gyro, const Vector3f& acc, float tempC = 35.0f) {
        if (!(dt > 0.0f) || !std::isfinite(dt) || !gyro.allFinite() || !acc.allFinite()) return;

        applyPendingTune_();
        elapsed_sec_ += dt;

        vertical_complementary_.update(dt, gyro, acc, cfg_.gravity_magnitude);
        last_acc_body_ = acc;
        last_gyro_body_ = gyro;
        have_last_imu_ = true;
        updateProxyGravityQuality_(dt, gyro, acc);

        const float a_up = vertical_complementary_.verticalAccelUpMs2();
        // OU-III feeds the canonical period estimator every sample.  Readiness
        // is a trust gate on the statistic, not permission to update it.
        if (std::isfinite(a_up)) wave_period_.update(dt, a_up);
        updateTuner_(dt, a_up);

        if (stage_ != StartupStage::Live) {
            tuner_warm_sec_ += dt;
        } else {
            live_since_sec_ += dt;
            maybeUnlockAccBias_();
        }

        if (stage_ == StartupStage::Cold) {
            if (cfg_.startup_init_policy == StartupInitPolicy::MahonyProxy) {
                tryProxyHandoff_();
            } else {
                stagedColdStep_(gyro, acc, dt);
            }
            return;
        }

        periodicAwCovSyncTick_(dt);
        mekf_.time_update(gyro, dt);
        mekf_.measurement_update_acc_only(acc, tempC);

        pseudo_elapsed_ += dt;
        if (pseudo_elapsed_ >= pseudo_period_sec_) {
            pseudo_elapsed_ = 0.0f;
            mekf_.applyIntegralZeroPseudoMeas();
        }

        if (stage_ == StartupStage::TunerWarm) {
            live_sec_ += dt;
            if (live_sec_ >= cfg_.online_tune_warmup_sec && wave_period_.isReady()) {
                enterLive_();
            }
        } else {
            adaptMekf_(dt);
        }
    }

    void updateMag(const Vector3f& mag_body) {
        if (!cfg_.with_mag || !mag_body.allFinite()) return;
        if (!(mag_body.norm() > 1e-9f)) return;
        if (elapsed_sec_ < cfg_.mag_delay_sec) return;

        accumulateContinuousHardIron_(mag_body);
        if (usingProxyInit_()) {
            proxyUpdateMag_(mag_body);
            return;
        }
        if (stage_ == StartupStage::Cold) return;

        if (!mekf_.has_magnetic_reference()) {
            const Vector3f b_world = mekf_.R_bw() * mag_body;
            const float horiz = std::hypot(b_world.x(), b_world.y());
            if (horiz > 1e-9f) {
                const float psi = std::atan2(b_world.y(), b_world.x());
                mekf_.apply_world_yaw_gauge(-psi);
            }
            mekf_.set_magnetic_reference_world(Vector3f(horiz, 0.0f, b_world.z()));
            return;
        }
        mekf_.measurement_update_mag_only(mag_body);
    }

    [[nodiscard]] bool magReferenceLearned() const noexcept { return mag_reference_learned_; }
    [[nodiscard]] bool magReferenceRefined() const noexcept { return mag_refine_done_; }
    [[nodiscard]] float magRefineTimeSec() const noexcept { return mag_refine_time_sec_; }
    [[nodiscard]] float magNorthLockTimeSec() const noexcept { return mag_north_lock_time_sec_; }
    [[nodiscard]] const Vector3f& magHardIronBodyUT() const noexcept { return mag_hard_iron_body_uT_; }
    [[nodiscard]] const Vector3f& magContinuousHardIronAppliedUT() const noexcept { return mag_hi_applied_body_uT_; }
    [[nodiscard]] const ContinuousMagHardIronEstimator& magContinuousHardIron() const noexcept { return mag_hi_estimator_; }
    [[nodiscard]] StartupInitPolicy startupInitPolicy() const noexcept { return cfg_.startup_init_policy; }
    [[nodiscard]] bool isTunerReady() const noexcept {
        return wave_period_.isReady() && tuner_warm_sec_ >= cfg_.online_tune_warmup_sec;
    }
    [[nodiscard]] float pseudoUpdatePeriodSec() const noexcept { return pseudo_period_sec_; }
    [[nodiscard]] bool handoffTimedOut() const noexcept { return handoff_timed_out_; }
    [[nodiscard]] float getRSFilterInput() const noexcept { return RS_filter_input_; }

    void setAdaptEverySecs(float s) { if (s >= 0.0f && std::isfinite(s)) adapt_every_secs_ = s; }
    void setTauScaledPseudoCadence(bool on) { tau_scaled_pseudo_cadence_ = on; applyPseudoCadence_(); }
    [[nodiscard]] bool tauScaledPseudoCadence() const noexcept { return tau_scaled_pseudo_cadence_; }

    bool setFixedTuning(float tau_s, float sigma_a, float RS) {
        if (!(tau_s > 0.0f) || !(sigma_a >= 0.0f) || !(RS > 0.0f)) return false;
        if (!std::isfinite(tau_s) || !std::isfinite(sigma_a) || !std::isfinite(RS)) return false;
        fixed_tuning_ = true;
        tune_.tau_applied = tau_s;
        tune_.sigma_applied = sigma_a;
        tune_.RS_applied = RS;
        tau_target_ = tau_s;
        sigma_target_ = sigma_a;
        RS_target_ = RS;
        commitTune_();
        return true;
    }

    bool setChannelFreeze(bool freeze_ou, float tau_s, float sigma_a,
                          bool freeze_RS, float RS) {
        if (freeze_ou && freeze_RS) return false;
        if (freeze_ou) {
            if (!(tau_s > 0.0f) || !(sigma_a >= 0.0f)) return false;
            freeze_ou_channel_ = true;
            tune_.tau_applied = tau_s;
            tune_.sigma_applied = sigma_a;
            tau_target_ = tau_s;
            sigma_target_ = sigma_a;
        }
        if (freeze_RS) {
            if (!(RS > 0.0f)) return false;
            freeze_RS_channel_ = true;
            tune_.RS_applied = RS;
            RS_target_ = RS;
        }
        commitTune_();
        return true;
    }

    void enableTuner(bool on) { enable_tuner_ = on; }
    void enableLinearBlock(bool on) { mekf_.set_linear_block_enabled(on); }
    void setTauCoeff(float c)      { if (c > 0.0f && std::isfinite(c)) tau_coeff_ = c; }
    void setSigmaCoeff(float c)    { if (c > 0.0f && std::isfinite(c)) sigma_coeff_ = c; }
    void setRSCoeff(float c)       { if (c > 0.0f && std::isfinite(c)) R_S_coeff_ = c; }
    void setSFactor(float s)       { if (s > 0.0f && std::isfinite(s)) S_factor_ = s; }
    void setRSXFactor(float k)     { if (k > 0.0f && std::isfinite(k)) R_S_x_factor_ = k; }
    void setRSYFactor(float k)     { if (k > 0.0f && std::isfinite(k)) R_S_y_factor_ = k; }
    [[nodiscard]] float getRSXFactor() const noexcept { return R_S_x_factor_; }
    [[nodiscard]] float getRSYFactor() const noexcept { return R_S_y_factor_; }
    void setAccNoiseFloorSigma(float s) { if (s >= 0.0f && std::isfinite(s)) noise_floor_sigma_ = s; }

    void setRSLaw(RSLaw law) noexcept { rs_law_ = law; }
    [[nodiscard]] RSLaw getRSLaw() const noexcept { return rs_law_; }
    void setEmbeddedFriendlyLegacyRSLaw(bool legacy) noexcept {
        rs_law_ = legacy ? RSLaw::LegacyCubic : RSLaw::SpectralMSE;
    }
    void setRSMseCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) rs_mse_coeff_ = c;
    }
    void setRSAccelNoiseDensity(float r_a) {
        if (std::isfinite(r_a) && r_a > 0.0f) {
            rs_accel_noise_density_ = r_a;
            refreshQeffPow_();
        }
    }
    [[nodiscard]] float getRSMseCoeff() const noexcept { return rs_mse_coeff_; }
    [[nodiscard]] float getRSAccelNoiseDensity() const noexcept { return rs_accel_noise_density_; }

    void setAdaptationTimeConstants(float tau_sec) {
        if (tau_sec > 0.0f && std::isfinite(tau_sec)) {
            adapt_tau_sec_ = tau_sec;
            adapt_tau_sea_periods_ = 0.0f;
        }
    }
    void setAdaptationSeaPeriods(float periods) {
        if (periods > 0.0f && std::isfinite(periods)) adapt_tau_sea_periods_ = periods;
    }
    void setTunerFreqSmoothingSeaPeriods(float /*periods*/) {}
    void setTunerFreqSmoothingTimeConstant(float /*tau_sec*/) {}
    [[nodiscard]] float getAdaptationSeaPeriods() const noexcept { return adapt_tau_sea_periods_; }
    [[nodiscard]] float getTunerFreqSmoothingSeaPeriods() const noexcept {
        return tuner_.getFrequencySmoothingSeaPeriods();
    }
    [[nodiscard]] float getSigmaVarianceHorizonPeriods() const noexcept { return tuner_.getKPeriods(); }
    [[nodiscard]] float getSigmaVarianceHorizonSec() const noexcept { return tuner_.getVarianceHorizonSec(); }

    void setRSAdaptMult(float m)     { if (m > 0.0f && std::isfinite(m)) adapt_RS_mult_ = m; }
    void setRSAdaptSlewLog(float d)  { if (d >= 0.0f && std::isfinite(d)) adapt_RS_slew_log_ = d; }
    [[nodiscard]] float getRSAdaptMult() const noexcept { return adapt_RS_mult_; }
    void setTauBounds(float lo, float hi) { if (lo > 0.0f && hi > lo) { min_tau_ = lo; max_tau_ = hi; } }
    void setRSBounds(float lo, float hi) { if (lo > 0.0f && hi > lo) { min_RS_ = lo; max_RS_ = hi; } }
    void setMaxSigmaA(float m)       { if (m > 0.0f && std::isfinite(m)) max_sigma_a_ = m; }
    void setWaveBandTuning(bool on)  { cfg_.wave_band_tuning = on; }
    [[nodiscard]] bool waveBandTuning() const noexcept { return cfg_.wave_band_tuning; }
    void setPeriodicAwCovSync(bool on) { periodic_aw_cov_sync_ = on; }
    [[nodiscard]] bool periodicAwCovSync() const noexcept { return periodic_aw_cov_sync_; }
    void setSigmaBandRatios(float low, float high) {
        cfg_.sigma_band_low_ratio = low;
        cfg_.sigma_band_high_ratio = high;
        sigma_wave_band_.setRatios(low, high);
    }

    [[nodiscard]] Mekf& mekf() noexcept { return mekf_; }
    [[nodiscard]] const Mekf& mekf() const noexcept { return mekf_; }
    [[nodiscard]] StartupStage stage() const noexcept { return stage_; }
    [[nodiscard]] bool isLive() const noexcept { return stage_ == StartupStage::Live; }
    [[nodiscard]] float getTauApplied()   const noexcept { return tune_.tau_applied; }
    [[nodiscard]] float getSigmaApplied() const noexcept { return tune_.sigma_applied; }
    [[nodiscard]] float getRSApplied()    const noexcept { return tune_.RS_applied; }
    [[nodiscard]] float getTauTarget()    const noexcept { return tau_target_; }
    [[nodiscard]] float getSigmaTarget()  const noexcept { return sigma_target_; }
    [[nodiscard]] float getRSTarget()     const noexcept { return RS_target_; }
    [[nodiscard]] float getWavePeriodSec() const noexcept { return wave_period_.getPeriodSec(); }
    [[nodiscard]] bool  wavePeriodReady() const noexcept { return wave_period_.isReady(); }
    [[nodiscard]] float getAccelVariance() const noexcept { return tuner_.getAccelVariance(); }
    [[nodiscard]] Eigen::Quaternionf quaternion() const { return mekf_.quaternion(); }
    [[nodiscard]] Vector3f get_velocity() const { return mekf_.get_velocity(); }
    [[nodiscard]] Vector3f get_position() const { return mekf_.get_position(); }
    [[nodiscard]] Vector3f get_world_accel() const { return mekf_.get_world_accel(); }

private:
    void beginMagAcquisition_() {
        ::MagAutoTuner::Config mag_cfg;
        mag_cfg.mag_norm_min             = cfg_.mag_init_min_mag_norm;
        mag_cfg.min_samples              = cfg_.mag_min_samples;
        mag_cfg.min_window_sec           = cfg_.mag_min_window_sec;
        mag_cfg.max_window_sec           = cfg_.mag_max_window_sec;
        mag_cfg.sample_dt_sec            = cfg_.mag_sample_dt_sec;
        mag_cfg.gravity_ref              = cfg_.gravity_magnitude;
        mag_cfg.enable_quality_weighting = cfg_.mag_enable_quality_weighting;
        mag_cfg.estimate_hard_iron       = cfg_.mag_estimate_hard_iron;
        mag_cfg.min_effective_weight     = cfg_.mag_min_effective_weight;
        mag_cfg.acc_norm_rel_soft        = cfg_.mag_acc_norm_rel_soft;
        mag_cfg.gyro_soft_dps            = cfg_.mag_gyro_soft_dps;
        mag_auto_tuner_.setConfig(mag_cfg);

        ContinuousMagHardIronEstimator::Config hi_cfg;
        hi_cfg.memory_sec           = cfg_.mag_hi_memory_sec;
        hi_cfg.model_ridge          = cfg_.mag_hi_model_ridge;
        hi_cfg.model_ridge_relative = cfg_.mag_hi_model_ridge_relative;
        hi_cfg.min_information      = cfg_.mag_hi_min_information;
        hi_cfg.min_effective_weight = cfg_.mag_hi_min_effective_weight;
        hi_cfg.max_residual_rms_uT  = cfg_.mag_hi_max_residual_rms_uT;
        hi_cfg.max_bias_fraction    = cfg_.mag_hi_max_bias_fraction;
        hi_cfg.min_mag_norm_uT      = cfg_.mag_init_min_mag_norm;
        mag_hi_estimator_.setConfig(hi_cfg);

        mag_reference_learned_ = false;
        mag_world_ref_valid_ = false;
        mag_world_ref_uT_.setZero();
        mag_yaw_anchor_rad_ = 0.0f;
        have_mag_yaw_anchor_ = false;
        mag_hard_iron_body_uT_.setZero();
        mag_hi_startup_body_uT_.setZero();
        mag_hi_applied_body_uT_.setZero();
        mag_hi_anchor_bias_body_uT_.setZero();
        mag_hi_anchor_world_ref_uT_.setZero();
        mag_hi_anchored_ = false;
        mag_refine_started_ = false;
        mag_refine_done_ = false;
        mag_refine_time_sec_ = NAN;
        mag_north_lock_time_sec_ = NAN;
        mag_init_eligible_t0_ = NAN;
        last_mag_sample_t_ = NAN;
        last_hi_sample_t_ = NAN;
        last_hi_apply_t_ = NAN;
        gravity_gate_acc_world_lpf_.reset();
        gravity_gate_world_elapsed_sec_ = 0.0f;
        proxy_gravity_good_sec_ = 0.0f;
        proxy_gravity_aligned_branch_ = false;
        last_acc_body_.setZero();
        last_gyro_body_.setZero();
        have_last_imu_ = false;
        acc_bias_unlocked_ = false;
        live_since_sec_ = 0.0f;
        acc_bias_hold_ = usingProxyInit_() && cfg_.with_mag && cfg_.mag_refine_enabled;
    }

    void enterCold_() {
        mekf_.set_linear_block_enabled(false);
        if (cfg_.freeze_acc_bias_until_live) mekf_.set_acc_bias_updates_enabled(false);
        mekf_.set_Racc_std(Vector3f::Constant(cfg_.Racc_warmup_std));
    }

    void enterLive_() {
        stage_ = StartupStage::Live;
        live_since_sec_ = 0.0f;
        mekf_.set_linear_block_enabled(true);
        acc_bias_unlocked_ = false;
        mekf_.set_acc_bias_updates_enabled(false);
        maybeUnlockAccBias_();
        mekf_.set_Racc_std(Racc_nominal_);
        commitTune_();
        mekf_.reset_aw_covariance_to_stationary();
    }

    void maybeUnlockAccBias_() {
        if (acc_bias_unlocked_) return;
        if (cfg_.freeze_acc_bias_until_live) {
            if (live_since_sec_ < cfg_.acc_bias_unlock_sec) return;
            if (acc_bias_hold_ && !mag_refine_done_) return;
        }
        acc_bias_unlocked_ = true;
        mekf_.set_acc_bias_updates_enabled(true);
    }

    [[nodiscard]] float bandNoiseFloorSigma_() const noexcept {
        if (!cfg_.wave_band_tuning || !sigma_wave_band_.isReady()) return noise_floor_sigma_;
        const float gain = sigma_wave_band_.whiteNoiseVarianceGain();
        if (!(std::isfinite(gain) && gain >= 0.0f)) return noise_floor_sigma_;
        return noise_floor_sigma_ * std::sqrt(gain);
    }

    [[nodiscard]] float wavePeriodFrequencyOrPrior_() const noexcept {
        const float f = wave_period_.getFrequencyHz();
        return (std::isfinite(f) && f > 0.0f) ? f : kTuneFreqPriorHz;
    }

    [[nodiscard]] float pseudoUpdatePeriodFor_(float tau) const noexcept {
        if (!tau_scaled_pseudo_cadence_) return kPseudoPeriodNominalS;
        if (!(std::isfinite(tau) && tau > 0.0f)) return kPseudoPeriodNominalS;
        return std::min(std::max(kPseudoTauRatio * tau, kPseudoPeriodMinS),
                        kPseudoPeriodMaxS);
    }

    [[nodiscard]] float rsSpectralMseTarget_(float tau, float sigma) const noexcept {
        const float TS = pseudoUpdatePeriodFor_(tau);
        const float c_sigma = (std::isfinite(sigma_coeff_) && sigma_coeff_ > 0.0f)
                                  ? sigma_coeff_ : 1.0f;
        const float sigma_aB = std::max(sigma / c_sigma, 1e-6f);
        if (!(TS > 0.0f) || !(tau > 0.0f)) return rs_mse_coeff_;
        const float tau2 = tau * tau;
        const float u = sigma_aB * tau2 * tau2;
        return rs_mse_coeff_ * rs_qeff_pow_
             * std::pow(u, 6.0f / 7.0f)
             / std::sqrt(TS);
    }

    [[nodiscard]] float rsTargetFromLaw_(float tau, float sigma) const noexcept {
        if (rs_law_ == RSLaw::LegacyCubic) {
            return R_S_coeff_ * sigma * tau * tau * tau;
        }
        return rsSpectralMseTarget_(tau, sigma);
    }

    void refreshQeffPow_() noexcept {
        float r_a = rs_accel_noise_density_;
        if (!(std::isfinite(r_a) && r_a > 0.0f)) r_a = R_S_ACCEL_NOISE_DENSITY_DEFAULT;
        rs_qeff_pow_ = std::pow(2.0f * r_a, 1.0f / 14.0f);
    }

    void updateTuner_(float dt, float a_up) {
        if (!enable_tuner_) return;

        const float f_hint = wavePeriodFrequencyOrPrior_();
        float f_band = tuner_.isFreqReady() ? tuner_.getFrequencyHz() : f_hint;
        if (!std::isfinite(f_band) || f_band < kMinTuneFreqHz) f_band = kMinTuneFreqHz;
        f_band = std::min(f_band, kMaxTuneFreqHz);

        const float a_for_variance = cfg_.wave_band_tuning
            ? sigma_wave_band_.step(a_up, dt, f_band)
            : a_up;
        tuner_.update(dt, a_for_variance, f_hint);
        if (fixed_tuning_) return;

        float f_tune = tuner_.getFrequencyHz();
        if (!std::isfinite(f_tune) || f_tune < kMinTuneFreqHz) f_tune = kMinTuneFreqHz;
        f_tune = std::min(f_tune, kMaxTuneFreqHz);
        sea_time_sec_ = 0.5f / f_tune;

        const float band_noise_sigma = bandNoiseFloorSigma_();
        const float var_noise = band_noise_sigma * band_noise_sigma;
        const float var_total = tuner_.isVarReady()
            ? std::max(0.0f, tuner_.getAccelVariance())
            : var_noise;
        const float sigma_wave = std::sqrt(std::max(1e-6f, var_total - var_noise));

        const float tau_live = std::min(max_tau_, std::max(min_tau_, tau_coeff_ * 0.5f / f_tune));
        float sigma_live = std::min(max_sigma_a_, sigma_coeff_ * sigma_wave);
        if (!tuner_.isVarReady()) {
            sigma_live = std::max(sigma_live, std::max(0.05f, band_noise_sigma));
        }

        // As in OU-III, r_S is derived from the live front-end estimate before
        // an ablation freezes either channel.  This keeps the partial-adaptation
        // arms orthogonal rather than silently freezing r_S with tau/sigma.
        const float rs_live = std::min(max_RS_, std::max(min_RS_,
                                                         rsTargetFromLaw_(tau_live, sigma_live)));

        if (!freeze_ou_channel_) {
            tau_target_ = tau_live;
            sigma_target_ = sigma_live;
        } else {
            tau_target_ = tune_.tau_applied;
            sigma_target_ = tune_.sigma_applied;
        }
        if (!freeze_RS_channel_) RS_target_ = rs_live;
    }

    void stagedColdStep_(const Vector3f& gyro, const Vector3f& acc, float dt) {
        // StagedMekf is a legacy ablation.  Keep its historical TFG bootstrap;
        // the deployed MahonyProxy path below is the parity path.
        const bool levelled = ::seastate::common::runStartupGravityInit(
            gyro, acc, dt, elapsed_sec_, cfg_.gravity_magnitude,
            0.5f, 1.0f, 0.12f, 0.5f, 0.5f, 8.0f, 0.15f,
            bootstrap_tilt_obs_, bootstrap_gravity_slow_lpf_,
            bootstrap_gravity_good_sec_,
            [&](const Eigen::Vector3f& g_body) { mekf_.initialize_from_acc(g_body); });
        if (!levelled) return;
        stage_ = StartupStage::TunerWarm;
        enterCold_();
        commitTune_();
        mekf_.reset_aw_covariance_to_stationary();
    }

    void tryProxyHandoff_() {
        if (elapsed_sec_ < cfg_.proxy_startup_min_sec) return;
        if (!vertical_complementary_.isInitialized()) return;

        const float mag_acquire_deadline = cfg_.with_mag
            ? cfg_.proxy_mag_settle_sec + 2.0f * std::max(cfg_.mag_min_window_sec, 1.0f)
                  + cfg_.mag_tilt_fallback_sec
            : 0.0f;
        const float timeout_sec = std::max(cfg_.proxy_startup_timeout_sec, mag_acquire_deadline);
        const bool timed_out = (elapsed_sec_ >= timeout_sec) && proxy_gravity_aligned_branch_;

        if (!timed_out) {
            if (!vertical_complementary_.isReady()) return;
            if (!proxyGravityTrusted_()) return;
            if (!isTunerReady()) return;
            if (cfg_.with_mag && !mag_reference_learned_) return;
        }

        const Eigen::Quaternionf q_tilt = vertical_complementary_.tiltQuaternion();
        const Eigen::Quaternionf q_bw = have_mag_yaw_anchor_
            ? Eigen::Quaternionf(Eigen::AngleAxisf(mag_yaw_anchor_rad_, Eigen::Vector3f::UnitZ()) * q_tilt)
            : q_tilt;

        mekf_.initialize_from_truth(q_bw.normalized(), Vector3f::Zero(), Vector3f::Zero(),
                                    Vector3f::Zero(), Vector3f::Zero(),
                                    Vector3f::Zero(), Vector3f::Zero());
        seedHandoffAttitudeCovariance_();
        mekf_.set_initial_linear_uncertainty(1.0f, 2.0f, 5.0f, 1.0f);
        mekf_.set_initial_acc_bias_std(cfg_.handoff_acc_bias_std);
        if (mag_world_ref_valid_) mekf_.set_magnetic_reference_world(mag_world_ref_uT_);
        handoff_timed_out_ = timed_out;
        enterLive_();
    }

    [[nodiscard]] bool usingProxyInit_() const noexcept {
        return cfg_.startup_init_policy == StartupInitPolicy::MahonyProxy;
    }

    void proxyUpdateMag_(const Vector3f& mag_body) {
        if (!mag_reference_learned_ && elapsed_sec_ < cfg_.proxy_mag_settle_sec) return;
        if (!vertical_complementary_.isInitialized()) return;
        if (!std::isfinite(mag_init_eligible_t0_)) mag_init_eligible_t0_ = elapsed_sec_;

        if (!mag_reference_learned_) {
            const bool fallback_ok = (elapsed_sec_ - mag_init_eligible_t0_) >= cfg_.mag_tilt_fallback_sec;
            if (!proxyGravityTrusted_() && !fallback_ok) return;
            if (!have_last_imu_) return;
            learnMagReferenceWindowed_(mag_body);
        }

        maybeRefineMagReference_(mag_body);
        maybeApplyContinuousHardIron_();
        if (mag_reference_learned_ && stage_ == StartupStage::Live) {
            mekf_.measurement_update_mag_only(mag_body - mag_hard_iron_body_uT_);
        }
    }

    void learnMagReferenceWindowed_(const Vector3f& mag_body) {
        const float dt_mag = (std::isfinite(last_mag_sample_t_) && elapsed_sec_ > last_mag_sample_t_)
            ? (elapsed_sec_ - last_mag_sample_t_) : cfg_.mag_sample_dt_sec;
        last_mag_sample_t_ = elapsed_sec_;

        if (!mag_auto_tuner_.addSampleWithTiltQuatDt(
                dt_mag, vertical_complementary_.tiltQuaternion(),
                last_acc_body_, last_gyro_body_, mag_body)) return;

        Vector3f ref;
        if (!mag_auto_tuner_.getMagWorldRef(ref) || !ref.allFinite() ||
            !(ref.norm() > cfg_.mag_init_min_mag_norm)) return;
        setMagWorldRef_(ref);

        const float gauge = mag_auto_tuner_.getYawGaugeCorrectionRad();
        if (std::isfinite(gauge)) {
            mag_yaw_anchor_rad_ = wrapPi_(-gauge);
            have_mag_yaw_anchor_ = true;
        }

        Vector3f hard_iron;
        mag_hi_startup_body_uT_ = mag_auto_tuner_.getHardIronBodyUT(hard_iron)
            ? hard_iron : Vector3f::Zero();
        mag_hard_iron_body_uT_ = mag_hi_startup_body_uT_ + mag_hi_applied_body_uT_;
        mag_reference_learned_ = true;
        mag_north_lock_time_sec_ = elapsed_sec_;
    }

    void maybeRefineMagReference_(const Vector3f& mag_body) {
        if (!cfg_.mag_refine_enabled || mag_refine_done_) return;
        if (!mag_reference_learned_ || stage_ != StartupStage::Live) return;
        if (elapsed_sec_ < cfg_.mag_refine_start_sec || !have_last_imu_) return;

        if (!mag_refine_started_) {
            ::MagAutoTuner::Config refine_cfg = mag_auto_tuner_.config();
            refine_cfg.min_window_sec = cfg_.mag_refine_window_sec;
            refine_cfg.min_samples = cfg_.mag_min_samples;
            mag_auto_tuner_.setConfig(refine_cfg);
            mag_auto_tuner_.reset();
            mag_refine_started_ = true;
            last_mag_sample_t_ = NAN;
        }

        const float dt_mag = (std::isfinite(last_mag_sample_t_) && elapsed_sec_ > last_mag_sample_t_)
            ? (elapsed_sec_ - last_mag_sample_t_) : cfg_.mag_sample_dt_sec;
        last_mag_sample_t_ = elapsed_sec_;
        const Vector3f mag_corrected = mag_body - mag_hard_iron_body_uT_;
        if (!mag_auto_tuner_.addSampleWithTiltQuatDt(
                dt_mag, vertical_complementary_.tiltQuaternion(),
                last_acc_body_, last_gyro_body_, mag_corrected)) return;

        Vector3f ref;
        if (!mag_auto_tuner_.getMagWorldRef(ref) || !ref.allFinite() ||
            !(ref.norm() > cfg_.mag_init_min_mag_norm)) return;
        setMagWorldRef_(ref);

        const float gauge = mag_auto_tuner_.getYawGaugeCorrectionRad();
        if (std::isfinite(gauge)) mekf_.set_attitude_yaw_absolute(wrapPi_(-gauge));
        mag_refine_done_ = true;
        mag_refine_time_sec_ = elapsed_sec_;
        maybeUnlockAccBias_();
    }

    void setMagWorldRef_(const Vector3f& ref) {
        mag_world_ref_uT_ = ref;
        mag_world_ref_valid_ = true;
        if (stage_ == StartupStage::Live || !usingProxyInit_()) {
            mekf_.set_magnetic_reference_world(ref);
        }
    }

    void accumulateContinuousHardIron_(const Vector3f& mag_body) {
        if (!cfg_.mag_continuous_hard_iron || !usingProxyInit_()) return;
        if (!vertical_complementary_.isInitialized()) return;
        const float dt_mag = (std::isfinite(last_hi_sample_t_) && elapsed_sec_ > last_hi_sample_t_)
            ? (elapsed_sec_ - last_hi_sample_t_) : cfg_.mag_sample_dt_sec;
        last_hi_sample_t_ = elapsed_sec_;
        mag_hi_estimator_.update(dt_mag, vertical_complementary_.tiltQuaternion(), mag_body);
    }

    void maybeApplyContinuousHardIron_() {
        if (!cfg_.mag_continuous_hard_iron || !usingProxyInit_()) return;
        if (!mag_reference_learned_ || stage_ != StartupStage::Live) return;
        if (cfg_.mag_refine_enabled && !mag_refine_done_) return;
        if (!mag_world_ref_valid_) return;
        const auto& est = mag_hi_estimator_.estimate();
        if (!est.valid) return;

        if (!mag_hi_anchored_) {
            mag_hi_anchor_bias_body_uT_ = mag_hi_applied_body_uT_;
            mag_hi_anchor_world_ref_uT_ = mag_world_ref_uT_;
            mag_hi_anchored_ = true;
        }
        const Vector3f target = cfg_.mag_hi_apply_fraction * est.bias_body_uT;
        if (!target.allFinite()) return;

        const float dt_apply = (std::isfinite(last_hi_apply_t_) && elapsed_sec_ > last_hi_apply_t_)
            ? (elapsed_sec_ - last_hi_apply_t_) : cfg_.mag_sample_dt_sec;
        last_hi_apply_t_ = elapsed_sec_;
        const float tau = cfg_.mag_hi_slew_tau_sec;
        const float alpha = (std::isfinite(tau) && tau > 1.0e-3f)
            ? (1.0f - std::exp(-dt_apply / tau)) : 1.0f;
        const Vector3f applied = mag_hi_applied_body_uT_ + alpha * (target - mag_hi_applied_body_uT_);
        if (!applied.allFinite()) return;

        Vector3f level_new, level_anchor;
        if (!mag_hi_estimator_.levelReferenceForBias(applied, level_new) ||
            !mag_hi_estimator_.levelReferenceForBias(mag_hi_anchor_bias_body_uT_, level_anchor)) return;

        const float h_new = level_new.head<2>().norm();
        const float h_anchor = level_anchor.head<2>().norm();
        if (!std::isfinite(h_new) || !std::isfinite(h_anchor)) return;
        const float h = mag_hi_anchor_world_ref_uT_.x() + (h_new - h_anchor);
        const float z = mag_hi_anchor_world_ref_uT_.z() + (level_new.z() - level_anchor.z());
        if (!(h > cfg_.mag_init_min_mag_norm) || !std::isfinite(h) || !std::isfinite(z)) return;

        mag_hi_applied_body_uT_ = applied;
        mag_hard_iron_body_uT_ = mag_hi_startup_body_uT_ + applied;
        setMagWorldRef_(Vector3f(h, 0.0f, z));
    }

    static float wrapPi_(float a) {
        constexpr float PI_F = 3.14159265358979323846f;
        if (!std::isfinite(a)) return NAN;
        while (a > PI_F) a -= 2.0f * PI_F;
        while (a <= -PI_F) a += 2.0f * PI_F;
        return a;
    }

    void updateProxyGravityQuality_(float dt, const Vector3f& gyro, const Vector3f& acc) {
        if (!vertical_complementary_.isInitialized()) {
            proxy_gravity_good_sec_ = 0.0f;
            proxy_gravity_aligned_branch_ = false;
            return;
        }

        // Match OU-III: rotate into the attitude's world frame first, then
        // average over the wave band.  Averaging body-frame specific force
        // makes hull roll/pitch smear the zero-mean orbital term into the gate.
        gravity_gate_acc_world_lpf_.step(
            ::seastate::common::accWorldFromBody(vertical_complementary_.quaternion(), acc),
            dt, cfg_.proxy_gravity_lpf_sec);
        gravity_gate_world_elapsed_sec_ += dt;

        const Vector3f acc_world_lp = gravity_gate_acc_world_lpf_.state;
        const bool average_warm = gravity_gate_world_elapsed_sec_ >= cfg_.proxy_gravity_warmup_sec;
        const float sin_res = average_warm
            ? ::seastate::common::gravityAlignResidualSinWorld(acc_world_lp)
            : 1.0f;
        const bool aligned_branch =
            ::seastate::common::gravityAlignedBranchWorld(acc_world_lp);
        proxy_gravity_aligned_branch_ = aligned_branch;

        const float gyro_dps = gyro.norm() * 57.295779513f;
        const bool extreme_motion = !std::isfinite(gyro_dps) ||
                                    (gyro_dps > cfg_.mag_extreme_gyro_dps);
        const bool good_now = std::isfinite(sin_res) &&
                              (sin_res <= cfg_.proxy_gravity_align_sin) &&
                              aligned_branch && !extreme_motion;
        if (good_now) proxy_gravity_good_sec_ = std::min(10.0f, proxy_gravity_good_sec_ + dt);
        else proxy_gravity_good_sec_ = std::max(0.0f, proxy_gravity_good_sec_ - 2.0f * dt);
    }

    [[nodiscard]] bool proxyGravityTrusted_() const {
        return proxy_gravity_aligned_branch_ &&
               (proxy_gravity_good_sec_ >= cfg_.proxy_gravity_hold_sec);
    }

    void seedHandoffAttitudeCovariance_() {
        auto& P = mekf_.covariance_full();
        const float st = std::max(1e-6f, cfg_.proxy_handoff_tilt_sigma_rad);
        const float sy = have_mag_yaw_anchor_
            ? std::max(1e-6f, cfg_.proxy_handoff_yaw_sigma_rad)
            : std::max(1e-6f, cfg_.proxy_handoff_yaw_sigma_free_rad);
        constexpr int PHI = Mekf::OFF_PHI;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < static_cast<int>(Mekf::NX); ++j) {
                P(PHI + i, j) = 0.0f;
                P(j, PHI + i) = 0.0f;
            }
        }
        P(PHI + 0, PHI + 0) = st * st;
        P(PHI + 1, PHI + 1) = st * st;
        P(PHI + 2, PHI + 2) = sy * sy;
    }

    void adaptMekf_(float dt) {
        if (!enable_tuner_ || fixed_tuning_) return;

        if (!freeze_ou_channel_) {
            float adapt_sec = adapt_tau_sec_;
            if (adapt_tau_sea_periods_ > 0.0f &&
                std::isfinite(sea_time_sec_) && sea_time_sec_ > 0.0f) {
                const float safe_sea_time =
                    seastate::tuner::limits::clampDynamicEmaTimeScaleSec(sea_time_sec_);
                adapt_sec = seastate::tuner::limits::clampDynamicEmaHorizonSec(
                    adapt_tau_sea_periods_ * safe_sea_time, dt);
            }
            const float a = 1.0f - std::exp(-dt / adapt_sec);
            tune_.tau_applied += a * (tau_target_ - tune_.tau_applied);
            tune_.sigma_applied += a * (sigma_target_ - tune_.sigma_applied);
        }
        if (!freeze_RS_channel_) {
            const float horizon = ::seastate::common::adaptiveSmoothingHorizonSec(
                adapt_RS_mult_, tau_target_, RS_target_, tune_.RS_applied,
                adapt_RS_slew_log_, dt);
            const float a = (horizon > 0.0f) ? (1.0f - std::exp(-dt / horizon)) : 1.0f;
            tune_.RS_applied += a * (RS_target_ - tune_.RS_applied);
        }

        adapt_elapsed_sec_ += dt;
        if (adapt_elapsed_sec_ >= adapt_every_secs_) {
            adapt_elapsed_sec_ = 0.0f;
            tune_apply_pending_ = true;
        }
    }

    void applyPendingTune_() {
        if (!tune_apply_pending_) return;
        tune_apply_pending_ = false;
        commitTune_();
    }

    void periodicAwCovSyncTick_(float dt) {
        if (!periodic_aw_cov_sync_ || stage_ != StartupStage::Live) return;
        aw_sync_elapsed_sec_ += dt;
        if (aw_sync_elapsed_sec_ < adapt_every_secs_) return;
        aw_sync_elapsed_sec_ = 0.0f;
        mekf_.synchronize_aw_covariance_to_stationary();
    }

    void applyPseudoCadence_() {
        pseudo_period_sec_ = pseudoUpdatePeriodFor_(tune_.tau_applied);
    }

    void commitTune_() {
        mekf_.set_aw_time_constant(tune_.tau_applied);
        applyPseudoCadence_();
        const float sZ = tune_.sigma_applied;
        mekf_.set_aw_stationary_std(Vector3f(sZ * S_factor_, sZ * S_factor_, sZ));

        // SpectralMSE already contains the realized target T_S.  The historical
        // cubic target does not, so only LegacyCubic receives the old
        // information-rate normalization here.
        float rs = tune_.RS_applied;
        if (rs_law_ == RSLaw::LegacyCubic && tau_scaled_pseudo_cadence_ && pseudo_period_sec_ > 0.0f) {
            rs *= std::sqrt(kPseudoPeriodNominalS / pseudo_period_sec_);
        }
        RS_filter_input_ = rs;
        mekf_.set_RS_noise(Vector3f(rs * R_S_x_factor_, rs * R_S_y_factor_, rs));
    }

    struct Vec3LPF {
        Eigen::Vector3f state = Eigen::Vector3f::Zero();
        bool initialized = false;
        void reset() { state.setZero(); initialized = false; }
        Eigen::Vector3f step(const Eigen::Vector3f& x, float dt, float tau_sec) {
            if (!x.allFinite()) return state;
            const float tau = std::max(1.0e-3f, tau_sec);
            const float alpha = 1.0f - std::exp(-dt / tau);
            if (!initialized) { state = x; initialized = true; return state; }
            state += alpha * (x - state);
            return state;
        }
    };

    static constexpr float kTuneFreqPriorHz = 0.2f;
    static constexpr float kMinTuneFreqHz = 0.03f;
    static constexpr float kMaxTuneFreqHz = 1.5f;
    static constexpr float kPseudoPeriodNominalS = 0.015f;
    static constexpr float kPseudoTauNominalS = 1.1f;
    static constexpr float kPseudoTauRatio = kPseudoPeriodNominalS / kPseudoTauNominalS;
    static constexpr float kPseudoPeriodMinS = 0.005f;
    static constexpr float kPseudoPeriodMaxS = 0.25f;

    Config cfg_{};
    Mekf mekf_{};
    StartupStage stage_{StartupStage::Cold};
    ::SeaStateAutoTuner tuner_{TUNER_SIGMA_VAR_K_PERIODS_DEFAULT, 1.0f};
    ::WavePeriodEstimator wave_period_{};
    ::VerticalAccelComplementary vertical_complementary_{};
    ::AdaptiveWaveBandPass sigma_wave_band_{SIGMA_BAND_LOW_RATIO_DEFAULT,
                                            SIGMA_BAND_HIGH_RATIO_DEFAULT,
                                            SIGMA_BAND_MIN_HZ_DEFAULT,
                                            SIGMA_BAND_MAX_HZ_DEFAULT};
    ::MagAutoTuner mag_auto_tuner_{};
    ContinuousMagHardIronEstimator mag_hi_estimator_{};
    ::seastate::common::StartupTiltObserver bootstrap_tilt_obs_{};
    Vec3LPF bootstrap_gravity_slow_lpf_{};
    float bootstrap_gravity_good_sec_ = 0.0f;

    TfgTuneState tune_{};
    float tau_target_ = 1.1f;
    float sigma_target_ = 1e-2f;
    float RS_target_ = 0.5f;

    // TFG-specific physical OU prior coefficients remain independently fitted.
    // Preserve the historical 1.15 horizontal regularization operating point,
    // while exposing X/Y as independent experiment/tuning knobs.
    float tau_coeff_ = 1.0f;
    float sigma_coeff_ = 1.0f;
    float R_S_coeff_ = 0.28f;
    float S_factor_ = 1.00f;
    float R_S_x_factor_ = 1.15f;
    float R_S_y_factor_ = 1.15f;
    float noise_floor_sigma_ = 0.12f;

    RSLaw rs_law_ = RSLaw::SpectralMSE;
    float rs_accel_noise_density_ = R_S_ACCEL_NOISE_DENSITY_DEFAULT;
    float rs_mse_coeff_ = R_S_MSE_COEFF_DEFAULT;
    float rs_qeff_pow_ = std::pow(2.0f * R_S_ACCEL_NOISE_DENSITY_DEFAULT, 1.0f / 14.0f);

    float adapt_tau_sec_ = ADAPT_TAU_SEC_DEFAULT;
    float adapt_tau_sea_periods_ = ADAPT_TAU_SEA_PERIODS_DEFAULT;
    float sea_time_sec_ = 0.5f / kTuneFreqPriorHz;
    float adapt_RS_mult_ = ADAPT_RS_MULT_DEFAULT;
    float adapt_RS_slew_log_ = 0.0f;

    float min_tau_ = 0.02f, max_tau_ = 12.0f;
    float min_RS_ = 0.15f, max_RS_ = 400.0f;
    float max_sigma_a_ = 6.0f;
    bool enable_tuner_ = true;
    bool fixed_tuning_ = false;
    bool freeze_ou_channel_ = false;
    bool freeze_RS_channel_ = false;
    Vector3f Racc_nominal_{Vector3f::Constant(0.5f)};

    float elapsed_sec_ = 0.0f;
    float live_sec_ = 0.0f;
    float mag_elapsed_sec_ = 0.0f;
    float pseudo_elapsed_ = 0.0f;
    float pseudo_period_sec_ = kPseudoPeriodNominalS;
    float adapt_every_secs_ = 0.1f;
    float adapt_elapsed_sec_ = 0.0f;
    bool tune_apply_pending_ = false;
    bool tau_scaled_pseudo_cadence_ = true;
    float RS_filter_input_ = 0.5f;

    float tuner_warm_sec_ = 0.0f;
    bool mag_reference_learned_ = false;
    bool handoff_timed_out_ = false;
    bool acc_bias_unlocked_ = false;
    bool acc_bias_hold_ = false;
    float live_since_sec_ = 0.0f;
    float proxy_gravity_good_sec_ = 0.0f;
    bool proxy_gravity_aligned_branch_ = false;
    Vec3LPF gravity_gate_acc_world_lpf_{};
    float gravity_gate_world_elapsed_sec_ = 0.0f;

    bool periodic_aw_cov_sync_ = true;
    float aw_sync_elapsed_sec_ = 0.0f;

    float mag_yaw_anchor_rad_ = 0.0f;
    bool have_mag_yaw_anchor_ = false;
    Vector3f mag_world_ref_uT_{Vector3f::Zero()};
    bool mag_world_ref_valid_ = false;
    float mag_init_eligible_t0_ = NAN;
    float last_mag_sample_t_ = NAN;
    float mag_north_lock_time_sec_ = NAN;
    bool mag_refine_started_ = false;
    bool mag_refine_done_ = false;
    float mag_refine_time_sec_ = NAN;

    Vector3f mag_hard_iron_body_uT_{Vector3f::Zero()};
    Vector3f mag_hi_startup_body_uT_{Vector3f::Zero()};
    Vector3f mag_hi_applied_body_uT_{Vector3f::Zero()};
    Vector3f mag_hi_anchor_bias_body_uT_{Vector3f::Zero()};
    Vector3f mag_hi_anchor_world_ref_uT_{Vector3f::Zero()};
    bool mag_hi_anchored_ = false;
    float last_hi_sample_t_ = NAN;
    float last_hi_apply_t_ = NAN;

    Vector3f last_acc_body_{Vector3f::Zero()};
    Vector3f last_gyro_body_{Vector3f::Zero()};
    bool have_last_imu_ = false;
};

}  // namespace ocean_imu::tfg
