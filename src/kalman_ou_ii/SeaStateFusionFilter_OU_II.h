#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy

  Released under the MIT License

  SeaStateFusionFilter_OU_II
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>
#include <memory>
#include <algorithm>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "freq/FirstOrderIIRSmoother.h"
#include "freq/FrequencyTrackerPolicy.h"
#include "tuner/SeaStateAutoTuner.h"
#include "tuner/MagAutoTuner.h"
#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
#include "wave_dir/KalmanWaveDirection.h"
#include "wave_dir/WaveDirectionDetector.h"
#include "detrend/AdaptiveWaveDetrender3D.h"
#include "kalman_common/SeaStateFusionFilterCommon.h"

// Shared constants
extern const float g_std;

#ifndef FREQ_GUESS
#define FREQ_GUESS 0.3f
#endif

#ifndef ZERO_CROSSINGS_SCALE
#define ZERO_CROSSINGS_SCALE 1.0f
#endif

#ifndef ZERO_CROSSINGS_DEBOUNCE_TIME
#define ZERO_CROSSINGS_DEBOUNCE_TIME 0.12f
#endif

#ifndef ZERO_CROSSINGS_STEEPNESS_TIME
#define ZERO_CROSSINGS_STEEPNESS_TIME 0.21f
#endif

constexpr float ACC_NOISE_FLOOR_SIGMA_DEFAULT = 0.12f;

constexpr float MIN_FREQ_HZ = 0.2f;
constexpr float MAX_FREQ_HZ = 6.0f;

constexpr float MIN_TAU_S     = 0.02f;
constexpr float MAX_TAU_S     = 3.0f;
constexpr float MAX_SIGMA_A   = 6.0f;
constexpr float MIN_R_p0_std  = 0.05f;
constexpr float MAX_R_p0_std  = 18.0f;
constexpr float MIN_R_v0_std  = 0.01f;
constexpr float MAX_R_v0_std  = 6.0f;

constexpr float ADAPT_TAU_SEC              = 1.8f;
constexpr float ADAPT_EVERY_SECS           = 0.1f;
constexpr float ADAPT_R_p0_MULT            = 5.0f;
constexpr float ADAPT_R_v0_MULT            = 5.0f;
constexpr float ONLINE_TUNE_WARMUP_SEC     = 5.0f;
constexpr float MAG_DELAY_SEC              = 7.0f;

constexpr float FREQ_SMOOTHER_DT = 1.0f / 200.0f;

struct TuneStateOU2 {
    float tau_applied      = 1.1f;
    float sigma_applied    = 1e-2f;
    float R_p0_std_applied = 0.1f;
    float R_v0_std_applied = 0.1f;
};

template<TrackerType trackerT>
class SeaStateFusionFilter_OU_II {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using TrackingPolicy = TrackerPolicy<trackerT>;

    enum class StartupStage {
        Cold,
        TunerWarm,
        Live
    };

    explicit SeaStateFusionFilter_OU_II(bool with_mag = true)
        : with_mag_(with_mag),
          time_(0.0),
          last_adapt_time_sec_(0.0),
          freq_hz_(FREQ_GUESS),
          freq_hz_slow_(FREQ_GUESS)
    {
        freq_input_lpf_.setCutoff(max_freq_hz_);
        freq_stillness_.setTargetFreqHz(min_freq_hz_);
        startup_stage_   = StartupStage::Cold;
        startup_stage_t_ = 0.0f;
    }

    StartupStage getStartupStage() const noexcept { return startup_stage_; }
    bool isAdaptiveLive() const noexcept { return startup_stage_ == StartupStage::Live; }

    void initialize(const Eigen::Vector3f& sigma_a,
                    const Eigen::Vector3f& sigma_g,
                    const Eigen::Vector3f& sigma_m)
    {
        mekf_ = std::make_unique<Kalman3D_Wave_OU_II<float>>(
            sigma_a, sigma_g, sigma_m);

        seastate::common::finalizeInitialization(
            mekf_,
            [this]() { enterCold_(); },
            [this]() { apply_ou_tune_(); });
    }

    void initialize_ext(const Eigen::Vector3f& sigma_a,
                        const Eigen::Vector3f& sigma_g,
                        const Eigen::Vector3f& sigma_m,
                        float Pq0,
                        float Pb0,
                        float b0,
                        float R_p0_var_init,
                        float R_v0_var_init,
                        float gravity_magnitude)
    {
        mekf_ = std::make_unique<Kalman3D_Wave_OU_II<float>>(
            sigma_a,
            sigma_g,
            sigma_m,
            Pq0,
            Pb0,
            b0,
            R_p0_var_init,
            R_v0_var_init,
            gravity_magnitude);

        seastate::common::finalizeInitialization(
            mekf_,
            [this]() { enterCold_(); },
            [this]() { apply_ou_tune_(); });
    }

    void initialize_from_acc(const Eigen::Vector3f& acc_body_ned) {
        if (mekf_) {
            mekf_->initialize_from_acc(acc_body_ned);
        }
    }

    void updateTime(float dt,
                    const Eigen::Vector3f& gyro,
                    const Eigen::Vector3f& acc,
                    float tempC = 35.0f)
    {
        if (!mekf_) return;
        if (!(dt > 0.0f) || !std::isfinite(dt)) return;

        time_ += dt;
        startup_stage_t_ += dt;

        const float a_x_body = acc.x();
        const float a_y_body = acc.y();

        const float a_z_body_proxy = acc.z() + g_std;

        mekf_->time_update(gyro, dt);
        mekf_->measurement_update_acc_only(acc, tempC);

        {
            Eigen::Quaternionf q_bw = mekf_->quaternion_boat();
            q_bw.normalize();

            const Eigen::Vector3f z_body_down_world =
                q_bw * Eigen::Vector3f(0.0f, 0.0f, 1.0f);

            const Eigen::Vector3f z_world_down(0.0f, 0.0f, 1.0f);

            float cos_tilt =
                z_body_down_world.normalized().dot(z_world_down);

            cos_tilt = std::max(-1.0f, std::min(1.0f, cos_tilt));

            const float tilt_deg =
                std::acos(cos_tilt) * 57.295779513f;

            constexpr float TILT_RESET_DEG = 70.0f;
            constexpr float TILT_RESET_HOLD_SEC = 0.35f;
            constexpr float TILT_RESET_COOLDOWN_SEC = 3.0f;

            if (tilt_reset_cooldown_sec_ > 0.0f) {
                tilt_reset_cooldown_sec_ =
                    std::max(0.0f, tilt_reset_cooldown_sec_ - dt);
            }

            if (tilt_deg > TILT_RESET_DEG) {
                tilt_over_limit_sec_ += dt;
            } else {
                tilt_over_limit_sec_ =
                    std::max(0.0f, tilt_over_limit_sec_ - 2.0f * dt);
            }

            if (tilt_over_limit_sec_ >= TILT_RESET_HOLD_SEC &&
                tilt_reset_cooldown_sec_ <= 0.0f)
            {
                if (startup_stage_ == StartupStage::Live) {
                    mekf_->initialize_from_acc_preserve_yaw(acc);
                } else {
                    mekf_->initialize_from_acc(acc);
                    enterCold_();
                    resetTrackingState_();
                }

                tilt_over_limit_sec_ = 0.0f;
                tilt_reset_cooldown_sec_ = TILT_RESET_COOLDOWN_SEC;
            }
        }

        a_body_z_up_proxy_ = -a_z_body_proxy;

        const float a_body_z_up_lp =
            freq_input_lpf_.step(a_body_z_up_proxy_, dt);

        const float f_tracker =
            static_cast<float>(tracker_policy_.run(a_body_z_up_lp, dt));

        f_raw = f_tracker;

        const float f_after_still =
            freq_stillness_.step(a_body_z_up_lp, dt, f_tracker);

        float f_fast = freq_fast_smoother_.update(f_after_still);
        float f_slow = freq_slow_smoother_.update(f_fast);

        f_fast = std::min(std::max(f_fast, min_freq_hz_), max_freq_hz_);
        f_slow = std::min(std::max(f_slow, min_freq_hz_), max_freq_hz_);

        freq_hz_      = f_fast;
        freq_hz_slow_ = f_slow;

        if (enable_tuner_) {
            update_tuner(dt, a_body_z_up_proxy_, f_after_still);
        }

        if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
            apply_R_p0_tune_();
            apply_R_v0_tune_();
        }

        const float omega = 2.0f * static_cast<float>(M_PI) * freq_hz_;

        dir_filter_.update(a_x_body, a_y_body, omega, dt);
        dir_sign_state_ =
            dir_sign_.update(a_x_body, a_y_body, a_body_z_up_proxy_, dt);
    }

    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!with_mag_ || !mekf_) return;
        if (time_ < mag_delay_sec_) return;

        mekf_->measurement_update_mag_only(mag_body_ned);
        ++mag_updates_applied_;

        if (!std::isfinite(first_mag_update_time_)) {
            first_mag_update_time_ = static_cast<float>(time_);
        }

        if (accel_bias_locked_ &&
            startup_stage_ == StartupStage::Live &&
            mag_updates_applied_ >= MAG_UPDATES_TO_UNLOCK &&
            std::isfinite(first_mag_update_time_) &&
            (static_cast<float>(time_) - first_mag_update_time_) > 1.0f)
        {
            accel_bias_locked_ = false;

            if (freeze_acc_bias_until_live_ &&
                startup_stage_ == StartupStage::Live)
            {
                restoreNominalRaccIfNeeded_();
                refreshAccelBiasLearning_();
            }
        }
    }

    void setWithMag(bool with_mag) {
        with_mag_ = with_mag;
    }

    void setAccelBiasLearningEnabled(bool en) {
        accel_bias_learning_enabled_ = en;
        refreshAccelBiasLearning_();
    }

    void setPFactor(float p) {
        if (std::isfinite(p) && p > 0.0f) {
            P_factor_ = p;
        }
    }

    void setR_p0_XYFactor(float k) {
        if (std::isfinite(k)) {
            R_p0_xy_factor_ = std::min(std::max(k, 0.0f), 1.0f);
        }
    }

    void setTauCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) {
            tau_coeff_ = c;
        }
    }

    void setSigmaCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) {
            sigma_coeff_ = c;
        }
    }

    void setR_p0_Coeff(float c) {
        if (std::isfinite(c) && c > 0.0f) {
            const float prev = R_p0_coeff_;
            R_p0_coeff_ = c;

            if (std::isfinite(prev) && prev > 0.0f) {
                const float scale = c / prev;

                if (std::isfinite(tune_.R_p0_std_applied) &&
                    tune_.R_p0_std_applied > 0.0f)
                {
                    tune_.R_p0_std_applied *= scale;
                }

                if (std::isfinite(R_p0_std_target_) &&
                    R_p0_std_target_ > 0.0f)
                {
                    R_p0_std_target_ *= scale;
                }

                if (enable_linear_block_) {
                    apply_R_p0_tune_();
                }
            }
        }
    }

    void setR_v0_Coeff(float c) {
        if (std::isfinite(c) && c > 0.0f) {
            const float prev = R_v0_coeff_;
            R_v0_coeff_ = c;

            if (std::isfinite(prev) && prev > 0.0f) {
                const float scale = c / prev;

                if (std::isfinite(tune_.R_v0_std_applied) &&
                    tune_.R_v0_std_applied > 0.0f)
                {
                    tune_.R_v0_std_applied *= scale;
                }

                if (std::isfinite(R_v0_std_target_) &&
                    R_v0_std_target_ > 0.0f)
                {
                    R_v0_std_target_ *= scale;
                }

                if (enable_linear_block_) {
                    apply_R_v0_tune_();
                }
            }
        }
    }

    void setAccNoiseFloorSigma(float s) {
        if (std::isfinite(s) && s > 0.0f) {
            acc_noise_floor_sigma_ = s;
        }
    }

    float getAccNoiseFloorSigma() const noexcept {
        return acc_noise_floor_sigma_;
    }

    void setFreqInputCutoffHz(float fc) {
        freq_input_lpf_.setCutoff(fc);
    }

    void enableClamp(bool flag = true) {
        enable_clamp_ = flag;
    }

    void enableTuner(bool flag = true) {
        enable_tuner_ = flag;
    }

    void enableLinearBlock(bool flag = true) {
        enable_linear_block_ = flag;

        if (mekf_) {
            const bool on_now =
                flag && (startup_stage_ == StartupStage::Live);

            mekf_->set_linear_block_enabled(on_now);
        }
    }

    void setFreqBounds(float min_hz, float max_hz) {
        if (!std::isfinite(min_hz) || !std::isfinite(max_hz)) return;
        if (min_hz <= 0.0f || max_hz <= min_hz) return;

        min_freq_hz_ = min_hz;
        max_freq_hz_ = max_hz;

        freq_stillness_.setTargetFreqHz(min_freq_hz_);
    }

    void setTauBounds(float min_tau_s, float max_tau_s) {
        if (!std::isfinite(min_tau_s) || !std::isfinite(max_tau_s)) return;
        if (min_tau_s <= 0.0f || max_tau_s <= min_tau_s) return;

        min_tau_s_ = min_tau_s;
        max_tau_s_ = max_tau_s;
    }

    void setMaxSigmaA(float max_sigma_a) {
        if (!std::isfinite(max_sigma_a) || max_sigma_a <= 0.0f) return;
        max_sigma_a_ = max_sigma_a;
    }

    void setR_p0_Bounds(float min_R_p0_std, float max_R_p0_std) {
        if (!std::isfinite(min_R_p0_std) ||
            !std::isfinite(max_R_p0_std)) return;

        if (min_R_p0_std <= 0.0f ||
            max_R_p0_std <= min_R_p0_std) return;

        MIN_R_p0_std_ = min_R_p0_std;
        MAX_R_p0_std_ = max_R_p0_std;
    }

    void setR_v0_Bounds(float min_R_v0_std, float max_R_v0_std) {
        if (!std::isfinite(min_R_v0_std) ||
            !std::isfinite(max_R_v0_std)) return;

        if (min_R_v0_std <= 0.0f ||
            max_R_v0_std <= min_R_v0_std) return;

        MIN_R_v0_std_ = min_R_v0_std;
        MAX_R_v0_std_ = max_R_v0_std;
    }

    void setAdaptationTimeConstants(float tau_sec) {
        if (std::isfinite(tau_sec) && tau_sec > 0.0f) {
            adapt_tau_sec_ = tau_sec;
        }
    }

    void setAdaptationUpdatePeriod(float every_sec) {
        if (std::isfinite(every_sec) && every_sec > 0.0f) {
            adapt_every_secs_ = every_sec;
        }
    }

    void setOnlineTuneWarmupSec(float warmup_sec) {
        if (std::isfinite(warmup_sec) && warmup_sec >= 0.0f) {
            online_tune_warmup_sec_ = warmup_sec;
        }
    }

    void setMagDelaySec(float delay_sec) {
        if (std::isfinite(delay_sec) && delay_sec >= 0.0f) {
            mag_delay_sec_ = delay_sec;
        }
    }

    void setFreezeAccBiasUntilLive(bool en) {
        freeze_acc_bias_until_live_ = en;
        refreshAccelBiasLearning_();
    }

    void setWarmupRaccStd(float r) {
        if (std::isfinite(r) && r > 0.0f) {
            Racc_warmup_std_ = r;
        }
    }

    void setNominalRaccStd(const Eigen::Vector3f& r) {
        Racc_nominal_std_ = r;
    }

    inline float getFreqHz() const noexcept {
        return freq_hz_;
    }

    inline float getFreqSlowHz() const noexcept {
        return freq_hz_slow_;
    }

    inline float getFreqRawHz() const noexcept {
        return f_raw;
    }

    inline float getTauApplied() const noexcept {
        return mekf_ ? mekf_->get_aw_time_constant() : NAN;
    }

    inline float getSigmaApplied() const noexcept {
        return mekf_ ? mekf_->get_aw_stationary_std().z() : NAN;
    }

    inline float getR_p0_std_applied() const noexcept {
        return mekf_ ? mekf_->get_Rp0_noise_std().z() : NAN;
    }

    inline float getR_v0_std_applied() const noexcept {
        return mekf_ ? mekf_->get_Rv0_noise_std().z() : NAN;
    }

    inline float getTauTarget() const noexcept {
        return tau_target_;
    }

    inline float getSigmaTarget() const noexcept {
        return sigma_target_;
    }

    inline float getR_p0_std_target() const noexcept {
        return R_p0_std_target_;
    }

    inline float getR_v0_std_target() const noexcept {
        return R_v0_std_target_;
    }

    inline float getPeriodSec() const noexcept {
        return (freq_hz_slow_ > 1e-6f) ? 1.0f / freq_hz_slow_ : NAN;
    }

    inline float getAccelVariance() const noexcept {
        return tuner_.getAccelVariance();
    }

    inline float getAccelVertical() const noexcept {
        return a_body_z_up_proxy_;
    }

    inline float getHeaveAbs() const noexcept {
        if (!mekf_) return NAN;
        return std::fabs(mekf_->get_position().z());
    }

    inline float getDisplacementScale(bool smoothed = true) const noexcept {
        const float tau   = smoothed ? tune_.tau_applied   : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;

        if (!std::isfinite(sigma) || !std::isfinite(tau)) {
            return NAN;
        }

        constexpr float C_HS = 2.0f * std::sqrt(2.0f) / (M_PI * M_PI);
        return C_HS * sigma * tau * tau / 2.0f;
    }

    float getVerticalSpeedEnvelopeMps(bool smoothed = true) const noexcept {
        const float tau   = smoothed ? tune_.tau_applied   : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;

        if (!(tau > 1e-6f) ||
            !std::isfinite(tau) ||
            !std::isfinite(sigma))
        {
            return NAN;
        }

        constexpr float K = std::sqrt(2.0f) / M_PI;
        const float v_env = K * sigma * tau;

        return std::isfinite(v_env) ? v_env : NAN;
    }

    inline WaveDirection getDirSignState() const noexcept {
        return dir_sign_state_;
    }

    inline float getWaveDirectionDeg() const noexcept {
        return dir_filter_.getDirectionDegrees();
    }

    inline auto& mekf() noexcept {
        return *mekf_;
    }

    inline const auto& mekf() const noexcept {
        return *mekf_;
    }

    inline KalmanWaveDirection& dir() noexcept {
        return dir_filter_;
    }

    inline const KalmanWaveDirection& dir() const noexcept {
        return dir_filter_;
    }

    inline WaveDirectionDetector<float>& dir_sign() noexcept {
        return dir_sign_;
    }

    inline const WaveDirectionDetector<float>& dir_sign() const noexcept {
        return dir_sign_;
    }

private:
    using FreqInputLPF = seastate::common::FreqInputLPF;
    using StillnessAdapter = seastate::common::StillnessAdapter;

    void refreshAccelBiasLearning_() {
        if (!mekf_) return;

        const bool allow =
            accel_bias_learning_enabled_ &&
            startup_stage_ == StartupStage::Live &&
            !accel_bias_locked_;

        mekf_->set_acc_bias_updates_enabled(allow);
    }

    void restoreNominalRaccIfNeeded_() {
        if (!mekf_) return;

        if (warmup_Racc_active_) {
            if (Racc_nominal_std_.allFinite() &&
                Racc_nominal_std_.maxCoeff() > 0.0f)
            {
                mekf_->set_Racc_std(Racc_nominal_std_);
            }

            warmup_Racc_active_ = false;
        }
    }

    void apply_ou_tune_() {
        if (!mekf_) return;

        mekf_->set_aw_time_constant(tune_.tau_applied);

        const float sigma_floor =
            std::max(0.05f, acc_noise_floor_sigma_);

        const float sZ = std::max(sigma_floor, tune_.sigma_applied);
        const float sH = sZ * P_factor_;

        mekf_->set_aw_stationary_std(Eigen::Vector3f(sH, sH, sZ));
    }

    void apply_R_p0_tune_(float rp_scale = 1.0f) {
        if (!mekf_) return;

        const float p =
            (std::isfinite(rp_scale) && rp_scale > 0.0f)
                ? std::min(rp_scale, 1.0f)
                : 1.0f;

        const float R_p0_b =
            std::min(std::max(tune_.R_p0_std_applied,
                              MIN_R_p0_std_),
                     MAX_R_p0_std_);

        const float rp_xy = R_p0_b * p * R_p0_xy_factor_;

        mekf_->set_Rp0_noise_std(
            Eigen::Vector3f(rp_xy, rp_xy, R_p0_b * p));
    }

    void apply_R_v0_tune_(float rv_scale = 1.0f) {
        if (!mekf_) return;

        const float p =
            (std::isfinite(rv_scale) && rv_scale > 0.0f)
                ? std::min(rv_scale, 1.0f)
                : 1.0f;

        const float R_v0_b =
            std::min(std::max(tune_.R_v0_std_applied,
                              MIN_R_v0_std_),
                     MAX_R_v0_std_);

        mekf_->set_Rv0_noise_std(
            Eigen::Vector3f::Constant(R_v0_b * p));
    }

    void update_tuner(float dt,
                      float a_body_z_up_proxy,
                      float freq_hz_for_tuner)
    {
        tuner_.update(dt, a_body_z_up_proxy, freq_hz_for_tuner);

        switch (startup_stage_) {
            case StartupStage::Cold:
                if (startup_stage_t_ >= online_tune_warmup_sec_) {
                    startup_stage_   = StartupStage::TunerWarm;
                    startup_stage_t_ = 0.0f;
                }
                return;

            case StartupStage::TunerWarm:
                if (!tuner_.isFreqReady()) return;

                if (tuner_.isReady()) {
                    enterLive_();
                }
                break;

            case StartupStage::Live:
                break;
        }

        float f_tune = tuner_.getFrequencyHz();

        if (!std::isfinite(f_tune) || f_tune < min_freq_hz_) {
            f_tune = min_freq_hz_;
        }

        if (f_tune > max_freq_hz_) {
            f_tune = max_freq_hz_;
        }

        float var_total =
            acc_noise_floor_sigma_ * acc_noise_floor_sigma_;

        if (tuner_.isVarReady()) {
            var_total = std::max(0.0f, tuner_.getAccelVariance());
        }

        const float var_noise =
            acc_noise_floor_sigma_ * acc_noise_floor_sigma_;

        float var_wave = var_total - var_noise;

        if (var_wave < 0.0f) {
            var_wave = 0.0f;
        }

        if (freq_stillness_.isStill()) {
            const float still_t =
                std::max(0.0f, freq_stillness_.getStillTime());

            constexpr float STILL_VAR_DECAY_SEC = 1.0f;

            float atten = std::exp(-still_t / STILL_VAR_DECAY_SEC);
            atten = std::min(std::max(atten, 0.0f), 1.0f);

            var_wave *= atten;
        }

        var_wave = std::max(var_wave, 1e-6f);

        const float sigma_wave = std::sqrt(var_wave);
        const float tau_raw = tau_coeff_ * 0.5f / f_tune;

        if (enable_clamp_) {
            tau_target_ =
                std::min(std::max(tau_raw, min_tau_s_), max_tau_s_);

            sigma_target_ =
                std::min(sigma_wave * sigma_coeff_, max_sigma_a_);
        } else {
            tau_target_   = tau_raw;
            sigma_target_ = sigma_wave;
        }

        if (!tuner_.isVarReady()) {
            sigma_target_ =
                std::max(sigma_target_,
                         std::max(0.05f, acc_noise_floor_sigma_));
        }

        const float R_p0_raw =
            R_p0_coeff_ * sigma_target_ * tau_target_ * tau_target_;

        const float R_v0_raw =
            R_v0_coeff_ * sigma_target_ * tau_target_;

        if (enable_clamp_) {
            R_p0_std_target_ =
                std::min(std::max(R_p0_raw, MIN_R_p0_std_),
                         MAX_R_p0_std_);

            R_v0_std_target_ =
                std::min(std::max(R_v0_raw, MIN_R_v0_std_),
                         MAX_R_v0_std_);
        } else {
            R_p0_std_target_ = R_p0_raw;
            R_v0_std_target_ = R_v0_raw;
        }

        adapt_mekf(dt,
                   tau_target_,
                   sigma_target_,
                   R_p0_std_target_,
                   R_v0_std_target_);
    }

    void adapt_mekf(float dt,
                    float tau_t,
                    float sigma_t,
                    float R_p0_t,
                    float R_v0_t)
    {
        const float alpha =
            1.0f - std::exp(-dt / adapt_tau_sec_);

        const float R_p0_sec =
            ADAPT_R_p0_MULT * tau_t;

        const float R_v0_sec =
            ADAPT_R_v0_MULT * tau_t;

        const float alpha_R_p0 =
            1.0f - std::exp(-dt / R_p0_sec);

        const float alpha_R_v0 =
            1.0f - std::exp(-dt / R_v0_sec);

        tune_.tau_applied +=
            alpha * (tau_t - tune_.tau_applied);

        tune_.sigma_applied +=
            alpha * (sigma_t - tune_.sigma_applied);

        tune_.R_p0_std_applied +=
            alpha_R_p0 * (R_p0_t - tune_.R_p0_std_applied);

        tune_.R_v0_std_applied +=
            alpha_R_v0 * (R_v0_t - tune_.R_v0_std_applied);

        if (time_ - last_adapt_time_sec_ > adapt_every_secs_) {
            if (tuner_.isFreqReady()) {
                apply_ou_tune_();
            }

            if (startup_stage_ == StartupStage::Live &&
                enable_linear_block_)
            {
                apply_R_p0_tune_();
                apply_R_v0_tune_();
            }

            last_adapt_time_sec_ = time_;
        }
    }

    void resetTrackingState_() {
        tracker_policy_ = TrackingPolicy{};
        freq_input_lpf_ = FreqInputLPF{};
        freq_stillness_ =
            StillnessAdapter(g_std, min_freq_hz_, FREQ_GUESS);

        freq_input_lpf_.setCutoff(max_freq_hz_);
        freq_stillness_.setTargetFreqHz(min_freq_hz_);

        tuner_.reset();

        freq_fast_smoother_ =
            FirstOrderIIRSmoother<float>(FREQ_SMOOTHER_DT, 3.5f);

        freq_slow_smoother_ =
            FirstOrderIIRSmoother<float>(FREQ_SMOOTHER_DT, 10.0f);

        freq_hz_      = FREQ_GUESS;
        freq_hz_slow_ = FREQ_GUESS;
        f_raw         = FREQ_GUESS;

        dir_filter_ =
            KalmanWaveDirection(
                2.0f * static_cast<float>(M_PI) * FREQ_GUESS);

        dir_sign_state_ = UNCERTAIN;

        last_adapt_time_sec_ = time_;
    }

    void enterCold_() {
        startup_stage_   = StartupStage::Cold;
        startup_stage_t_ = 0.0f;

        if (!mekf_) return;

        mekf_->set_linear_block_enabled(false);

        accel_bias_locked_ =
            with_mag_;

        mag_updates_applied_ = 0;
        first_mag_update_time_ = NAN;

        if (freeze_acc_bias_until_live_) {
            mekf_->set_acc_bias_updates_enabled(false);
            mekf_->set_Racc_std(
                Eigen::Vector3f::Constant(Racc_warmup_std_));

            warmup_Racc_active_ = true;
        } else {
            refreshAccelBiasLearning_();
        }
    }

    void enterLive_() {
        startup_stage_   = StartupStage::Live;
        startup_stage_t_ = 0.0f;

        if (!mekf_) return;

        mekf_->set_linear_block_enabled(enable_linear_block_);

        if (freeze_acc_bias_until_live_) {
            restoreNominalRaccIfNeeded_();
        }

        refreshAccelBiasLearning_();

        apply_ou_tune_();

        if (enable_linear_block_) {
            apply_R_p0_tune_();
            apply_R_v0_tune_();
        }
    }

private:
    StartupStage startup_stage_ = StartupStage::Cold;
    float startup_stage_t_ = 0.0f;

    bool freeze_acc_bias_until_live_ = true;
    bool accel_bias_learning_enabled_ = false;

    float Racc_warmup_std_ = 0.6f;
    bool warmup_Racc_active_ = false;
    Eigen::Vector3f Racc_nominal_std_ =
        Eigen::Vector3f::Constant(0.0f);

    bool accel_bias_locked_ = true;
    int mag_updates_applied_ = 0;
    static constexpr int MAG_UPDATES_TO_UNLOCK = 200;

    bool with_mag_;
    double time_;
    double last_adapt_time_sec_;

    float first_mag_update_time_ = NAN;

    float tilt_over_limit_sec_ = 0.0f;
    float tilt_reset_cooldown_sec_ = 0.0f;

    float freq_hz_ = FREQ_GUESS;
    float freq_hz_slow_ = FREQ_GUESS;
    float f_raw = FREQ_GUESS;

    float a_body_z_up_proxy_ = 0.0f;

    bool enable_clamp_ = true;
    bool enable_tuner_ = true;
    bool enable_linear_block_ = true;

    float min_freq_hz_ = MIN_FREQ_HZ;
    float max_freq_hz_ = MAX_FREQ_HZ;
    float min_tau_s_ = MIN_TAU_S;
    float max_tau_s_ = MAX_TAU_S;
    float max_sigma_a_ = MAX_SIGMA_A;

    float MIN_R_p0_std_ = MIN_R_p0_std;
    float MAX_R_p0_std_ = MAX_R_p0_std;
    float MIN_R_v0_std_ = MIN_R_v0_std;
    float MAX_R_v0_std_ = MAX_R_v0_std;

    float adapt_tau_sec_ = ADAPT_TAU_SEC;
    float adapt_every_secs_ = ADAPT_EVERY_SECS;
    float online_tune_warmup_sec_ = ONLINE_TUNE_WARMUP_SEC;
    float mag_delay_sec_ = MAG_DELAY_SEC;

    float R_p0_xy_factor_ = 0.31f;
    float P_factor_ = 1.5f;

    TrackingPolicy tracker_policy_{};

    FirstOrderIIRSmoother<float> freq_fast_smoother_{
        FREQ_SMOOTHER_DT,
        3.5f
    };

    FirstOrderIIRSmoother<float> freq_slow_smoother_{
        FREQ_SMOOTHER_DT,
        10.0f
    };

    SeaStateAutoTuner tuner_;
    TuneStateOU2 tune_;

    float tau_target_ = NAN;
    float sigma_target_ = NAN;
    float R_p0_std_target_ = NAN;
    float R_v0_std_target_ = NAN;

    float acc_noise_floor_sigma_ = ACC_NOISE_FLOOR_SIGMA_DEFAULT;

    float R_p0_coeff_ = 1.6f;
    float R_v0_coeff_ = 1.4f;
    float tau_coeff_ = 1.5f;
    float sigma_coeff_ = 0.85f;

    std::unique_ptr<Kalman3D_Wave_OU_II<float>> mekf_;

    KalmanWaveDirection dir_filter_{
        2.0f * static_cast<float>(M_PI) * FREQ_GUESS
    };

    FreqInputLPF freq_input_lpf_;
    StillnessAdapter freq_stillness_;

    WaveDirectionDetector<float> dir_sign_{0.002f, 0.005f};
    WaveDirection dir_sign_state_ = UNCERTAIN;
};

template<TrackerType trackerT>
class SeaStateFusion_OU_II {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    struct Config {
        bool with_mag = true;

        float mag_delay_sec = MAG_DELAY_SEC;
        float online_tune_warmup_sec = 10.0f;

        bool freeze_acc_bias_until_live = true;
        bool enable_acc_bias_learning = false;

        float Racc_warmup_std = 1.2f;

        Eigen::Vector3f sigma_a =
            Eigen::Vector3f(0.2f, 0.2f, 0.2f);

        Eigen::Vector3f sigma_g =
            Eigen::Vector3f(0.01f, 0.01f, 0.01f);

        Eigen::Vector3f sigma_m =
            Eigen::Vector3f(0.3f, 0.3f, 0.3f);

        // Optional gate. Off by default because accel-LPF gravity can be wrong in waves.
        bool mag_require_gravity_gate = false;
        float mag_gravity_align_max_sin = 0.070f;
        float mag_gravity_align_hold_sec = 2.0f;
        float mag_gravity_align_lpf_tau = 1.0f;
        float mag_tilt_fallback_sec = 30.0f;
        float mag_extreme_gyro_dps = 45.0f;
        float mag_init_min_mag_norm = 1e-3f;

        // Match OU3 default mag acquisition.
        int mag_min_samples = 1500;
        float mag_min_window_sec = 10.0f;
        float mag_max_window_sec = 0.0f;
        float mag_sample_dt_sec = 1.0f / 200.0f;

        // Off by default in waves. Accel/gyro weighting can phase-select samples.
        bool mag_enable_quality_weighting = false;
        float mag_min_effective_weight = 0.0f;
        float mag_acc_norm_rel_soft = 0.22f;
        float mag_gyro_soft_dps = 45.0f;

        // Dedicated yaw-free gravity direction observer for mag leveling.
        float mag_tilt_obs_acc_tau_sec = 2.5f;
        float mag_tilt_obs_norm_frac = 0.22f;

        // Existing startup tilt observer for initializing the main filter.
        float bootstrap_tilt_obs_acc_tau_sec = 2.50f;
        float bootstrap_gravity_slow_tau_sec = 8.0f;
        float bootstrap_gravity_align_max_sin = 0.070f;
        float bootstrap_gravity_hold_sec = 2.0f;
        float bootstrap_gravity_min_sec = 5.0f;
        float bootstrap_gravity_timeout_sec = 15.0f;
        float bootstrap_gravity_norm_frac = 0.22f;

        bool enable_displacement_detrend = false;
        bool use_custom_displacement_detrend_cfg = false;
        AdaptiveWaveDetrender3D::Config displacement_detrend_cfg{};
    };

    void begin(const Config& cfg) {
        cfg_ = cfg;

        begun_ = true;
        stage_ = Stage::Uninitialized;
        t_ = 0.0f;

        gravity_gate_acc_lpf_.reset();
        mag_gravity_good_sec_ = 0.0f;
        mag_init_eligible_t0_ = NAN;
        last_mag_sample_t_ = NAN;

        mag_ref_set_ = false;
        mag_gravity_obs_.reset();

        last_mag_yaw_level_rad_ = NAN;
        last_mag_yaw_rel_rad_ = NAN;
        last_mag_yaw_correction_rad_ = NAN;

        MagAutoTuner::Config mag_cfg;
        mag_cfg.mag_norm_min = cfg_.mag_init_min_mag_norm;
        mag_cfg.min_samples = cfg_.mag_min_samples;
        mag_cfg.min_window_sec = cfg_.mag_min_window_sec;
        mag_cfg.max_window_sec = cfg_.mag_max_window_sec;
        mag_cfg.sample_dt_sec = cfg_.mag_sample_dt_sec;
        mag_cfg.gravity_ref = g_std;
        mag_cfg.enable_quality_weighting =
            cfg_.mag_enable_quality_weighting;
        mag_cfg.min_effective_weight =
            cfg_.mag_min_effective_weight;
        mag_cfg.acc_norm_rel_soft = cfg_.mag_acc_norm_rel_soft;
        mag_cfg.gyro_soft_dps = cfg_.mag_gyro_soft_dps;

        mag_auto_tuner_.setConfig(mag_cfg);

        resetTiltInit_();

        last_acc_body_ned_.setZero();
        last_gyro_body_ned_.setZero();
        have_last_imu_ = false;

        impl_.setWithMag(cfg_.with_mag);
        impl_.setFreezeAccBiasUntilLive(cfg_.freeze_acc_bias_until_live);
        impl_.setAccelBiasLearningEnabled(cfg_.enable_acc_bias_learning);
        impl_.setWarmupRaccStd(cfg_.Racc_warmup_std);
        impl_.setMagDelaySec(0.0f);
        impl_.setOnlineTuneWarmupSec(cfg_.online_tune_warmup_sec);

        impl_.initialize(cfg_.sigma_a, cfg_.sigma_g, cfg_.sigma_m);
        last_impl_startup_stage_ = impl_.getStartupStage();

        impl_.setNominalRaccStd(cfg_.sigma_a);

        displacement_up_m_.setZero();
        displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};

        if (cfg_.enable_displacement_detrend) {
            if (cfg_.use_custom_displacement_detrend_cfg) {
                displacement_detrender_.setConfig(
                    cfg_.displacement_detrend_cfg);
            } else {
                displacement_detrender_.setConfig(
                    seastate::common::defaultDisplacementDetrenderConfig<
                        AdaptiveWaveDetrender3D::Config>(FREQ_GUESS));
            }

            displacement_detrender_.reset(0.0f, 0.0f, 0.0f);
        }
    }

    void update(float dt,
                const Eigen::Vector3f& gyro_body_ned,
                const Eigen::Vector3f& acc_body_ned,
                float tempC = 35.0f)
    {
        if (!begun_) return;
        if (!(dt > 0.0f) || !std::isfinite(dt)) return;

        t_ += dt;

        mag_gravity_obs_.update(
            gyro_body_ned,
            acc_body_ned,
            dt,
            g_std,
            cfg_.mag_tilt_obs_acc_tau_sec,
            cfg_.mag_tilt_obs_norm_frac);

        if (stage_ == Stage::Uninitialized) {
            const bool tilt_ready =
                seastate::common::runStartupGravityInit(
                    gyro_body_ned,
                    acc_body_ned,
                    dt,
                    t_,
                    g_std,
                    cfg_.bootstrap_tilt_obs_acc_tau_sec,
                    cfg_.bootstrap_gravity_slow_tau_sec,
                    cfg_.bootstrap_gravity_align_max_sin,
                    cfg_.bootstrap_gravity_hold_sec,
                    cfg_.bootstrap_gravity_min_sec,
                    cfg_.bootstrap_gravity_timeout_sec,
                    cfg_.bootstrap_gravity_norm_frac,
                    bootstrap_tilt_obs_,
                    bootstrap_gravity_slow_lpf_,
                    bootstrap_gravity_good_sec_,
                    [this](const Eigen::Vector3f& acc_init) {
                        impl_.initialize_from_acc(acc_init);
                        mag_gravity_obs_.resetFromAcc(acc_init);
                    });

            if (tilt_ready) {
                stage_ = Stage::Warming;
            }
        }

        last_acc_body_ned_ = acc_body_ned;
        last_gyro_body_ned_ = gyro_body_ned;
        have_last_imu_ = true;

        if (stage_ != Stage::Uninitialized) {
            impl_.updateTime(dt, gyro_body_ned, acc_body_ned, tempC);

            const Eigen::Vector3f acc_gate_lp =
                gravity_gate_acc_lpf_.step(
                    acc_body_ned,
                    dt,
                    cfg_.mag_gravity_align_lpf_tau);

            const float align_sin =
                seastate::common::gravityAlignResidualSin(
                    impl_.mekf().quaternion_boat(),
                    acc_gate_lp);

            const float gyro_dps =
                gyro_body_ned.norm() * 57.295779513f;

            const bool extreme_motion =
                !std::isfinite(gyro_dps) ||
                (gyro_dps > cfg_.mag_extreme_gyro_dps);

            const bool gravity_good_now =
                std::isfinite(align_sin) &&
                (align_sin <= cfg_.mag_gravity_align_max_sin) &&
                !extreme_motion;

            if (gravity_good_now) {
                mag_gravity_good_sec_ += dt;

                if (mag_gravity_good_sec_ > 10.0f) {
                    mag_gravity_good_sec_ = 10.0f;
                }
            } else {
                mag_gravity_good_sec_ =
                    std::max(0.0f, mag_gravity_good_sec_ - 2.0f * dt);
            }

            const Eigen::Vector3f pos_ned_m =
                impl_.mekf().get_position();

            displacement_up_m_ =
                Eigen::Vector3f(
                    pos_ned_m.x(),
                    pos_ned_m.y(),
                    -pos_ned_m.z());

            if (cfg_.enable_displacement_detrend) {
                const float wave_hz = impl_.getFreqHz();

                const bool ext_freq_valid =
                    isLive() &&
                    std::isfinite(wave_hz) &&
                    (wave_hz >=
                     displacement_detrender_.config().min_wave_freq_hz) &&
                    (wave_hz <=
                     displacement_detrender_.config().max_wave_freq_hz);

                displacement_det_out_ =
                    displacement_detrender_.update(
                        displacement_up_m_,
                        dt,
                        wave_hz,
                        ext_freq_valid);
            } else {
                displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};
                displacement_det_out_.input = displacement_up_m_;
                displacement_det_out_.baseline_slow =
                    Eigen::Vector3f::Zero();
                displacement_det_out_.wave_raw = displacement_up_m_;
                displacement_det_out_.wave_clean = displacement_up_m_;
            }
        }

        const auto cur_stage = impl_.getStartupStage();

        if (cur_stage != last_impl_startup_stage_) {
            if (cur_stage ==
                SeaStateFusionFilter_OU_II<trackerT>::StartupStage::Cold)
            {
                mag_ref_set_ = false;
                mag_auto_tuner_.reset();

                gravity_gate_acc_lpf_.reset();
                mag_gravity_good_sec_ = 0.0f;
                mag_init_eligible_t0_ = NAN;
                last_mag_sample_t_ = NAN;

                last_mag_yaw_level_rad_ = NAN;
                last_mag_yaw_rel_rad_ = NAN;
                last_mag_yaw_correction_rad_ = NAN;

                if (stage_ != Stage::Live) {
                    stage_ = Stage::Warming;

                    displacement_up_m_.setZero();
                    displacement_det_out_ =
                        AdaptiveWaveDetrender3D::Output{};

                    if (cfg_.enable_displacement_detrend) {
                        displacement_detrender_.reset(0.0f, 0.0f, 0.0f);
                    }
                }
            }

            last_impl_startup_stage_ = cur_stage;
        }

        if (stage_ == Stage::Warming && impl_.isAdaptiveLive()) {
            stage_ = Stage::Live;
        }
    }

    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!begun_ || !cfg_.with_mag) return;
        if (stage_ == Stage::Uninitialized) return;
        if (t_ < cfg_.mag_delay_sec) return;
        if (!mag_gravity_obs_.ready()) return;

        if (!std::isfinite(mag_init_eligible_t0_)) {
            mag_init_eligible_t0_ = t_;
        }

        const bool gravity_trusted =
            (mag_gravity_good_sec_ >= cfg_.mag_gravity_align_hold_sec);

        const bool fallback_ok =
            ((t_ - mag_init_eligible_t0_) >= cfg_.mag_tilt_fallback_sec);

        if (cfg_.mag_require_gravity_gate &&
            !gravity_trusted &&
            !fallback_ok)
        {
            return;
        }

        if (!mag_ref_set_) {
            if (!have_last_imu_) return;

            const float dt_mag =
                (std::isfinite(last_mag_sample_t_) &&
                 t_ > last_mag_sample_t_)
                    ? (t_ - last_mag_sample_t_)
                    : cfg_.mag_sample_dt_sec;

            last_mag_sample_t_ = t_;

            const Eigen::Vector3f down_body =
                mag_gravity_obs_.downBody();

            if (mag_auto_tuner_.addSampleWithGravityDirDt(
                    dt_mag,
                    down_body,
                    last_acc_body_ned_,
                    last_gyro_body_ned_,
                    mag_body_ned))
            {
                Eigen::Vector3f mag_world_ref_uT;

                if (mag_auto_tuner_.getMagWorldRef(mag_world_ref_uT) &&
                    mag_world_ref_uT.allFinite() &&
                    mag_world_ref_uT.norm() > cfg_.mag_init_min_mag_norm)
                {
                    impl_.mekf().set_mag_world_ref(mag_world_ref_uT);

                    applyInitialMagYawGaugeCorrection_();

                    mag_ref_set_ = true;
                }
            }
        }

        if (mag_ref_set_) {
            impl_.updateMag(mag_body_ned);
        }
    }

    bool hasMagNorthLock() const noexcept {
        return mag_ref_set_;
    }

    bool isLive() const {
        return stage_ == Stage::Live;
    }

    float freqHz() const {
        return impl_.getFreqHz();
    }

    float waveDirectionDeg() const {
        return impl_.getWaveDirectionDeg();
    }

    const Eigen::Vector3f& displacementUpMeters() const {
        return displacement_up_m_;
    }

    const AdaptiveWaveDetrender3D::Output& displacementDetrend() const {
        return displacement_det_out_;
    }

    SeaStateFusionFilter_OU_II<trackerT>& raw() {
        return impl_;
    }

    const SeaStateFusionFilter_OU_II<trackerT>& raw() const {
        return impl_;
    }

    int magAcceptedCount() const noexcept {
        return mag_auto_tuner_.acceptedCount();
    }

    int magRejectedCount() const noexcept {
        return mag_auto_tuner_.rejectedCount();
    }

    float magAcceptedWindowSec() const noexcept {
        return mag_auto_tuner_.acceptedWindowSec();
    }

    float magEffectiveWeight() const noexcept {
        return mag_auto_tuner_.effectiveWeight();
    }

    float magYawLevelDeg() const noexcept {
        return std::isfinite(last_mag_yaw_level_rad_)
            ? last_mag_yaw_level_rad_ * 57.29577951308232f
            : NAN;
    }

    float magYawRelativeDeg() const noexcept {
        return std::isfinite(last_mag_yaw_rel_rad_)
            ? last_mag_yaw_rel_rad_ * 57.29577951308232f
            : NAN;
    }

    float magYawCorrectionDeg() const noexcept {
        return std::isfinite(last_mag_yaw_correction_rad_)
            ? last_mag_yaw_correction_rad_ * 57.29577951308232f
            : NAN;
    }

private:
    enum class Stage {
        Uninitialized,
        Warming,
        Live
    };

    struct Vec3LPF {
        Eigen::Vector3f state = Eigen::Vector3f::Zero();
        bool initialized = false;

        void reset() {
            state.setZero();
            initialized = false;
        }

        Eigen::Vector3f step(const Eigen::Vector3f& x,
                             float dt,
                             float tau_sec)
        {
            if (!x.allFinite()) return state;

            const float tau = std::max(1.0e-3f, tau_sec);
            const float alpha = 1.0f - std::exp(-dt / tau);

            if (!initialized) {
                state = x;
                initialized = true;
                return state;
            }

            state += alpha * (x - state);
            return state;
        }
    };

    struct MagGravityObserver {
        Eigen::Vector3f down_body =
            Eigen::Vector3f(0.0f, 0.0f, 1.0f);

        bool initialized = false;

        void reset() {
            down_body = Eigen::Vector3f(0.0f, 0.0f, 1.0f);
            initialized = false;
        }

        bool ready() const noexcept {
            return initialized &&
                   down_body.allFinite() &&
                   down_body.norm() > 1.0e-6f;
        }

        Eigen::Vector3f downBody() const {
            return down_body;
        }

        void resetFromAcc(const Eigen::Vector3f& acc_body_ned) {
            if (!acc_body_ned.allFinite()) return;

            const float an = acc_body_ned.norm();

            if (!(an > 1.0e-6f) || !std::isfinite(an)) return;

            Eigen::Vector3f d = -acc_body_ned / an;

            if (!d.allFinite() || !(d.norm() > 1.0e-6f)) return;

            down_body = d.normalized();
            initialized = true;
        }

        void update(const Eigen::Vector3f& gyro_body_rad_s,
                    const Eigen::Vector3f& acc_body_ned,
                    float dt,
                    float g_ref,
                    float acc_tau_sec,
                    float norm_frac)
        {
            if (!(dt > 0.0f) || !std::isfinite(dt)) return;

            if (!initialized) {
                resetFromAcc(acc_body_ned);
                return;
            }

            if (!down_body.allFinite() ||
                !(down_body.norm() > 1.0e-6f))
            {
                resetFromAcc(acc_body_ned);
                return;
            }

            down_body.normalize();

            if (gyro_body_rad_s.allFinite()) {
                down_body += -dt * gyro_body_rad_s.cross(down_body);

                const float pn = down_body.norm();

                if (pn > 1.0e-6f && std::isfinite(pn)) {
                    down_body /= pn;
                }
            }

            if (!acc_body_ned.allFinite()) return;

            const float an = acc_body_ned.norm();

            if (!(an > 1.0e-6f) || !std::isfinite(an)) return;

            Eigen::Vector3f down_meas = -acc_body_ned / an;

            if (!down_meas.allFinite()) return;

            const float g = std::max(g_ref, 1.0e-6f);
            const float rel = std::fabs(an - g) / g;

            const float soft = std::max(norm_frac, 1.0e-3f);

            float w = 1.0f - rel / soft;
            w = std::min(std::max(w, 0.0f), 1.0f);

            const float tau = std::max(acc_tau_sec, 1.0e-3f);
            const float alpha = w * (1.0f - std::exp(-dt / tau));

            down_body += alpha * (down_meas - down_body);

            const float cn = down_body.norm();

            if (cn > 1.0e-6f && std::isfinite(cn)) {
                down_body /= cn;
            }
        }

        Eigen::Quaternionf levelQuatBodyToWorld() const {
            if (!ready()) {
                return Eigen::Quaternionf::Identity();
            }

            Eigen::Vector3f a = down_body.normalized();
            const Eigen::Vector3f b(0.0f, 0.0f, 1.0f);

            float d = a.dot(b);
            d = std::min(std::max(d, -1.0f), 1.0f);

            if (d > 1.0f - 1.0e-6f) {
                return Eigen::Quaternionf::Identity();
            }

            if (d < -1.0f + 1.0e-6f) {
                return Eigen::Quaternionf(
                    Eigen::AngleAxisf(
                        float(M_PI),
                        Eigen::Vector3f(1.0f, 0.0f, 0.0f)));
            }

            Eigen::Vector3f axis = a.cross(b);
            const float axis_n = axis.norm();

            if (!(axis_n > 1.0e-6f) || !axis.allFinite()) {
                return Eigen::Quaternionf::Identity();
            }

            axis /= axis_n;

            Eigen::Quaternionf q(
                Eigen::AngleAxisf(std::acos(d), axis));

            q.normalize();
            return q;
        }
    };

    using StartupTiltObserver = seastate::common::StartupTiltObserver;

    void resetTiltInit_() {
        bootstrap_tilt_obs_.reset();
        bootstrap_gravity_slow_lpf_.reset();
        bootstrap_gravity_good_sec_ = 0.0f;
    }

    static float wrapPi_(float a) {
        constexpr float PI_F = 3.14159265358979323846f;
        constexpr float TWO_PI_F = 2.0f * PI_F;

        if (!std::isfinite(a)) return NAN;

        while (a > PI_F) {
            a -= TWO_PI_F;
        }

        while (a <= -PI_F) {
            a += TWO_PI_F;
        }

        return a;
    }

    static float yawFromQuatRad_(const Eigen::Quaternionf& q_in) {
        if (!q_in.coeffs().allFinite()) return NAN;

        Eigen::Quaternionf q = q_in;
        const float qn = q.norm();

        if (!(qn > 1.0e-6f) || !std::isfinite(qn)) {
            return NAN;
        }

        q.normalize();

        const Eigen::Matrix3f R = q.toRotationMatrix();

        const float c = R(0, 0);
        const float s = R(1, 0);
        const float h2 = c * c + s * s;

        if (!(h2 > 1.0e-8f) ||
            !std::isfinite(h2) ||
            !std::isfinite(c) ||
            !std::isfinite(s))
        {
            return NAN;
        }

        return std::atan2(s, c);
    }

    static float relativeYawLevelToCurrentWorldRad_(
        const Eigen::Quaternionf& q_body_to_current_world,
        const Eigen::Quaternionf& q_body_to_level_world)
    {
        if (!q_body_to_current_world.coeffs().allFinite() ||
            !q_body_to_level_world.coeffs().allFinite())
        {
            return NAN;
        }

        Eigen::Quaternionf q_bw = q_body_to_current_world;
        Eigen::Quaternionf q_level = q_body_to_level_world;

        q_bw.normalize();
        q_level.normalize();

        Eigen::Quaternionf q_rel = q_bw * q_level.conjugate();
        q_rel.normalize();

        return yawFromQuatRad_(q_rel);
    }

    void applyInitialMagYawGaugeCorrection_() {
        const float yaw_level_rad =
            mag_auto_tuner_.getYawGaugeCorrectionRad();

        if (!std::isfinite(yaw_level_rad)) return;

        Eigen::Quaternionf q_bw =
            impl_.mekf().quaternion_boat();

        if (!q_bw.coeffs().allFinite()) return;
        q_bw.normalize();

        const Eigen::Quaternionf q_level =
            mag_gravity_obs_.levelQuatBodyToWorld();

        if (!q_level.coeffs().allFinite()) return;

        const float yaw_rel_rad =
            relativeYawLevelToCurrentWorldRad_(q_bw, q_level);

        if (!std::isfinite(yaw_rel_rad)) return;

        const float yaw_corr_rad =
            wrapPi_(-(yaw_rel_rad + yaw_level_rad));

        const Eigen::Quaternionf q_corr(
            Eigen::AngleAxisf(
                yaw_corr_rad,
                Eigen::Vector3f::UnitZ()));

        Eigen::Quaternionf q_new = q_corr * q_bw;
        q_new.normalize();

        if (!q_new.coeffs().allFinite()) return;

        impl_.mekf().set_quaternion_boat(q_new);

        last_mag_yaw_level_rad_ = yaw_level_rad;
        last_mag_yaw_rel_rad_ = yaw_rel_rad;
        last_mag_yaw_correction_rad_ = yaw_corr_rad;
    }

private:
    Config cfg_{};
    SeaStateFusionFilter_OU_II<trackerT> impl_{false};

    bool begun_ = false;

    Stage stage_ = Stage::Uninitialized;
    float t_ = 0.0f;

    typename SeaStateFusionFilter_OU_II<trackerT>::StartupStage
        last_impl_startup_stage_ =
            SeaStateFusionFilter_OU_II<trackerT>::StartupStage::Cold;

    Eigen::Vector3f last_acc_body_ned_ =
        Eigen::Vector3f::Zero();

    Eigen::Vector3f last_gyro_body_ned_ =
        Eigen::Vector3f::Zero();

    bool have_last_imu_ = false;

    bool mag_ref_set_ = false;
    MagAutoTuner mag_auto_tuner_{};
    MagGravityObserver mag_gravity_obs_{};

    float last_mag_sample_t_ = NAN;

    float last_mag_yaw_level_rad_ = NAN;
    float last_mag_yaw_rel_rad_ = NAN;
    float last_mag_yaw_correction_rad_ = NAN;

    AdaptiveWaveDetrender3D displacement_detrender_{};
    AdaptiveWaveDetrender3D::Output displacement_det_out_{};
    Eigen::Vector3f displacement_up_m_ = Eigen::Vector3f::Zero();

    Vec3LPF gravity_gate_acc_lpf_{};
    float mag_gravity_good_sec_ = 0.0f;
    float mag_init_eligible_t0_ = NAN;

    StartupTiltObserver bootstrap_tilt_obs_{};
    Vec3LPF bootstrap_gravity_slow_lpf_{};
    float bootstrap_gravity_good_sec_ = 0.0f;
};
