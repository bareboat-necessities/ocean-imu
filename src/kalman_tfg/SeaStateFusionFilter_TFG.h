#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    Sea-state orchestration for the two-frame Lie-group filter.

    SCOPE. This is deliberately the core of what SeaStateFusionFilter_OU_III
    does, not all of it: startup staging, the sea-state tuner, the adaptation
    smoothers, and the update loop. The wave-direction estimator, the
    alternative frequency trackers (KalmANF, Aranovskiy, PLL), the displacement
    detrender and the wind-heel B' frame are NOT carried over. They are
    orthogonal to the estimator under test, and duplicating them would mean
    maintaining a second copy of logic whose behaviour is already established
    and evidenced for OU-III.

    What that costs is worth stating plainly: the TFG simulator is therefore
    not feature-comparable to the OU-III one on direction and frequency
    reporting, and any paired study has to compare the channels both filters
    actually implement.

    WHAT IS SHARED RATHER THAN COPIED
      - ::seastate::common::runStartupGravityInit and StartupTiltObserver, for
        the Cold-stage levelling;
      - ::seastate::common::adaptiveSmoothingHorizonSec, for the r_S channel;
      - SeaStateAutoTuner              (variance at the wave-band operating point)
      - WavePeriodEstimator            (zero-crossing T_z)
      - VerticalAccelComplementary     (its private Mahony observer)

    THE TUNER STAYS OUTSIDE THE GROUP STATE. tau, sigma_aw and r_S are
    hyperparameters here, not manifold states. The article's ablation shows
    nearly all the vertical benefit comes from adapting the integral
    regularization scale rather than the two OU parameters, so promoting them
    to stochastic states would add complexity without evidence it solves the
    problem that matters.

    WHY THE WAVE-PERIOD INPUT IS THE COMPLEMENTARY OBSERVER. The tuner's
    frequency channel must not be a function of the estimator it tunes, or the
    loop closes on itself. VerticalAccelComplementary runs its own small Mahony
    observer on raw gyro and accelerometer, so the period estimate is
    independent of the filter's attitude. That is the same choice OU-III
    settled on, and for the same reason.
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>

#include "kalman_common/SeaStateFusionFilterCommon.h"
#include "kalman_tfg/Kalman3D_Wave_TFG.h"
#include "tuner/SeaStateAutoTuner.h"
#include "tuner/VerticalAccelComplementary.h"
#include "tuner/WavePeriodEstimator.h"

namespace ocean_imu::tfg {

constexpr float MAG_DELAY_SEC = 7.0f;

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

    struct Config {
        Vector3f sigma_a{Vector3f::Constant(0.5f)};
        float    gyro_noise_density = 0.005f;
        Vector3f sigma_m{Vector3f::Constant(0.1f)};
        bool     with_mag = true;
        float    mag_delay_sec = MAG_DELAY_SEC;
        bool     freeze_acc_bias_until_live = true;
        float    Racc_warmup_std = 0.5f;
        float    gravity_magnitude = 9.80665f;
        float    online_tune_warmup_sec = 20.0f;
    };

    void begin(const Config& cfg) {
        cfg_ = cfg;
        mekf_ = Mekf(cfg.gyro_noise_density, cfg.gravity_magnitude);
        mekf_.initialize_identity();
        mekf_.set_Racc_std(cfg.sigma_a);
        mekf_.set_Rmag_std(cfg.sigma_m);
        Racc_nominal_ = cfg.sigma_a;

        tuner_ = ::SeaStateAutoTuner(2.0f, 1.0f);
        wave_period_.reset();
        vertical_complementary_.reset();
        bootstrap_tilt_obs_.reset();
        bootstrap_gravity_good_sec_ = 0.0f;
        elapsed_sec_ = 0.0f;
        live_sec_ = 0.0f;
        mag_elapsed_sec_ = 0.0f;
        stage_ = StartupStage::Cold;

        enterCold_();
        applyTune_();
        mekf_.reset_aw_covariance_to_stationary();
    }

    void update(float dt, const Vector3f& gyro, const Vector3f& acc, float tempC = 35.0f) {
        if (!(dt > 0.0f) || !std::isfinite(dt) || !gyro.allFinite() || !acc.allFinite()) return;
        elapsed_sec_ += dt;

        vertical_complementary_.update(dt, gyro, acc, cfg_.gravity_magnitude);
        const float a_up = vertical_complementary_.verticalAccelUpMs2();
        if (vertical_complementary_.isReady()) wave_period_.update(dt, a_up);

        if (stage_ == StartupStage::Cold) {
            const bool levelled = ::seastate::common::runStartupGravityInit(
                gyro, acc, dt, elapsed_sec_,
                cfg_.gravity_magnitude,
                0.5f,
                1.0f,
                0.12f,
                0.5f,
                0.5f,
                8.0f,
                0.15f,
                bootstrap_tilt_obs_, bootstrap_gravity_slow_lpf_,
                bootstrap_gravity_good_sec_,
                [&](const Eigen::Vector3f& g_body) {
                    mekf_.initialize_from_acc(g_body);
                });
            if (!levelled) return;
            stage_ = StartupStage::TunerWarm;
            enterCold_();
            applyTune_();
            mekf_.reset_aw_covariance_to_stationary();
            return;
        }

        mekf_.time_update(gyro, dt);
        mekf_.measurement_update_acc_only(acc, tempC);

        pseudo_elapsed_ += dt;
        if (pseudo_elapsed_ >= pseudo_period_sec_) {
            pseudo_elapsed_ = 0.0f;
            mekf_.applyIntegralZeroPseudoMeas();
        }

        updateTuner_(dt, a_up);

        if (stage_ == StartupStage::TunerWarm) {
            live_sec_ += dt;
            if (live_sec_ >= cfg_.online_tune_warmup_sec && wave_period_.isReady()) {
                enterLive_();
            }
        } else {
            // Adapt the process model (tau and Sigma_aw), not the posterior
            // covariance. P remains the Riccati covariance of the invariant
            // error and is changed only by propagation, measurements, reset
            // coordinate transport, or explicit initialization.
            adaptMekf_(dt);
        }
    }

    /*
        Magnetometer first lock is a gauge choice, not a measurement update.
        The filter therefore applies the yaw change through
        apply_world_yaw_gauge(), which rotates R, every world column, all
        world-referred tuning covariances, and the corresponding tangent blocks
        of P. Directly changing only R would leave the group state and its
        covariance in different world frames.
    */
    void updateMag(const Vector3f& mag_body) {
        if (!cfg_.with_mag || !mag_body.allFinite()) return;
        if (!(mag_body.norm() > 1e-9f)) return;
        if (stage_ == StartupStage::Cold) return;
        if (elapsed_sec_ < cfg_.mag_delay_sec) return;

        if (!mekf_.has_magnetic_reference()) {
            const Vector3f b_world = mekf_.R_bw() * mag_body;
            const float horiz = std::hypot(b_world.x(), b_world.y());
            if (horiz > 1e-9f) {
                const float psi = std::atan2(b_world.y(), b_world.x());
                mekf_.apply_world_yaw_gauge(-psi);
            }
            // z is unchanged by the yaw gauge and the horizontal magnitude is
            // gauge-invariant, so this canonical reference is consistent with
            // the transformed state.
            mekf_.set_magnetic_reference_world(Vector3f(horiz, 0.0f, b_world.z()));
            return;
        }
        mekf_.measurement_update_mag_only(mag_body);
    }

    bool setFixedTuning(float tau_s, float sigma_a, float RS) {
        if (!(tau_s > 0.0f) || !(sigma_a >= 0.0f) || !(RS > 0.0f)) return false;
        if (!std::isfinite(tau_s) || !std::isfinite(sigma_a) || !std::isfinite(RS)) return false;
        fixed_tuning_ = true;
        tune_.tau_applied = tau_s;
        tune_.sigma_applied = sigma_a;
        tune_.RS_applied = RS;
        tau_target_ = tau_s; sigma_target_ = sigma_a; RS_target_ = RS;
        applyTune_();
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
        }
        if (freeze_RS) {
            if (!(RS > 0.0f)) return false;
            freeze_RS_channel_ = true;
            tune_.RS_applied = RS;
        }
        applyTune_();
        return true;
    }

    void enableTuner(bool on) { enable_tuner_ = on; }
    void enableLinearBlock(bool on) { mekf_.set_linear_block_enabled(on); }

    void setTauCoeff(float c)      { if (c > 0.0f && std::isfinite(c)) tau_coeff_ = c; }
    void setSigmaCoeff(float c)    { if (c > 0.0f && std::isfinite(c)) sigma_coeff_ = c; }
    void setRSCoeff(float c)       { if (c > 0.0f && std::isfinite(c)) R_S_coeff_ = c; }
    void setSFactor(float s)       { if (s > 0.0f && std::isfinite(s)) S_factor_ = s; }
    void setRSXYFactor(float k)    { if (k > 0.0f && std::isfinite(k)) R_S_xy_factor_ = k; }
    void setAccNoiseFloorSigma(float s) { if (s >= 0.0f && std::isfinite(s)) noise_floor_sigma_ = s; }
    void setAdaptationTimeConstants(float tau_sec) {
        if (tau_sec > 0.0f && std::isfinite(tau_sec)) adapt_tau_sec_ = tau_sec;
    }
    void setRSAdaptMult(float m)     { if (m > 0.0f && std::isfinite(m)) adapt_RS_mult_ = m; }
    void setRSAdaptSlewLog(float d)  { if (d >= 0.0f && std::isfinite(d)) adapt_RS_slew_log_ = d; }
    void setTauBounds(float lo, float hi) { if (lo > 0.0f && hi > lo) { min_tau_ = lo; max_tau_ = hi; } }
    void setRSBounds(float lo, float hi)  { if (lo > 0.0f && hi > lo) { min_RS_ = lo; max_RS_ = hi; } }
    void setMaxSigmaA(float m)       { if (m > 0.0f && std::isfinite(m)) max_sigma_a_ = m; }

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
    [[nodiscard]] bool  wavePeriodReady()  const noexcept { return wave_period_.isReady(); }
    [[nodiscard]] float getAccelVariance() const noexcept { return tuner_.getAccelVariance(); }

    [[nodiscard]] Eigen::Quaternionf quaternion() const { return mekf_.quaternion(); }
    [[nodiscard]] Vector3f get_velocity() const { return mekf_.get_velocity(); }
    [[nodiscard]] Vector3f get_position() const { return mekf_.get_position(); }
    [[nodiscard]] Vector3f get_world_accel() const { return mekf_.get_world_accel(); }

private:
    void enterCold_() {
        mekf_.set_linear_block_enabled(false);
        if (cfg_.freeze_acc_bias_until_live) mekf_.set_acc_bias_updates_enabled(false);
        mekf_.set_Racc_std(Vector3f::Constant(cfg_.Racc_warmup_std));
    }

    void enterLive_() {
        stage_ = StartupStage::Live;
        mekf_.set_linear_block_enabled(true);
        mekf_.set_acc_bias_updates_enabled(true);
        mekf_.set_Racc_std(Racc_nominal_);
        applyTune_();
        mekf_.reset_aw_covariance_to_stationary();
    }

    void updateTuner_(float dt, float a_up) {
        if (!enable_tuner_) return;

        const float f_hint = wave_period_.isReady() ? wave_period_.getFrequencyHz() : 0.2f;
        tuner_.update(dt, a_up, f_hint);
        if (fixed_tuning_) return;

        const float f_tune = std::max(kMinTuneFreqHz,
                                      wave_period_.isReady() ? wave_period_.getFrequencyHz()
                                                             : kMinTuneFreqHz);
        const float var = tuner_.getAccelVariance();
        const float wave_var = std::max(0.0f, var - noise_floor_sigma_ * noise_floor_sigma_);
        const float sigma_wave = std::sqrt(wave_var);

        if (!freeze_ou_channel_) {
            tau_target_ = std::min(max_tau_, std::max(min_tau_, tau_coeff_ * 0.5f / f_tune));
            sigma_target_ = std::min(max_sigma_a_, sigma_coeff_ * sigma_wave);
        }
        if (!freeze_RS_channel_) {
            const float tau_for_RS = freeze_ou_channel_ ? tune_.tau_applied : tau_target_;
            const float sig_for_RS = freeze_ou_channel_ ? tune_.sigma_applied : sigma_target_;
            const float raw = R_S_coeff_ * sig_for_RS * tau_for_RS * tau_for_RS * tau_for_RS;
            RS_target_ = std::min(max_RS_, std::max(min_RS_, raw));
        }
    }

    void adaptMekf_(float dt) {
        if (!enable_tuner_ || fixed_tuning_) return;

        if (!freeze_ou_channel_) {
            const float a = 1.0f - std::exp(-dt / adapt_tau_sec_);
            tune_.tau_applied   += a * (tau_target_ - tune_.tau_applied);
            tune_.sigma_applied += a * (sigma_target_ - tune_.sigma_applied);
        }
        if (!freeze_RS_channel_) {
            const float horizon = ::seastate::common::adaptiveSmoothingHorizonSec(
                adapt_RS_mult_, tune_.tau_applied, RS_target_, tune_.RS_applied,
                adapt_RS_slew_log_, dt);
            const float a = (horizon > 0.0f) ? (1.0f - std::exp(-dt / horizon)) : 1.0f;
            tune_.RS_applied += a * (RS_target_ - tune_.RS_applied);
        }
        applyTune_();
    }

    void applyTune_() {
        mekf_.set_aw_time_constant(tune_.tau_applied);
        const float sZ = tune_.sigma_applied;
        mekf_.set_aw_stationary_std(Vector3f(sZ * S_factor_, sZ * S_factor_, sZ));
        const float rs = tune_.RS_applied;
        mekf_.set_RS_noise(Vector3f(rs * R_S_xy_factor_, rs * R_S_xy_factor_, rs));
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

    static constexpr float kMinTuneFreqHz = 0.03f;

    Config cfg_{};
    Mekf   mekf_{};
    StartupStage stage_{StartupStage::Cold};

    ::SeaStateAutoTuner tuner_{2.0f, 1.0f};
    ::WavePeriodEstimator wave_period_{};
    ::VerticalAccelComplementary vertical_complementary_{};

    ::seastate::common::StartupTiltObserver bootstrap_tilt_obs_{};
    Vec3LPF bootstrap_gravity_slow_lpf_{};
    float bootstrap_gravity_good_sec_ = 0.0f;

    TfgTuneState tune_{};
    float tau_target_   = 1.1f;
    float sigma_target_ = 1e-2f;
    float RS_target_    = 0.5f;

    float tau_coeff_    = 1.0f;
    float sigma_coeff_  = 1.0f;
    float R_S_coeff_    = 0.35f;
    float S_factor_     = 1.87f;
    float R_S_xy_factor_ = 1.0f;
    float noise_floor_sigma_ = 0.12f;

    float adapt_tau_sec_    = 1.8f;
    float adapt_RS_mult_    = 3.0f;
    float adapt_RS_slew_log_ = 0.0f;

    float min_tau_ = 0.02f, max_tau_ = 12.0f;
    float min_RS_  = 0.4f,  max_RS_  = 400.0f;
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
    float pseudo_period_sec_ = 0.015f;
};

}  // namespace ocean_imu::tfg
