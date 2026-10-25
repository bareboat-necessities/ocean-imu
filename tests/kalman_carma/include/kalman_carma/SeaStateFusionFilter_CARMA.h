#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy
  Released under the MIT License

  SeaStateFusionFilter_CARMA -- EXPERIMENTAL research wrapper around
  Kalman3D_Wave_CARMA.  Not a deployed estimator: SeaStateFusionFilter_OU_III
  stays the shipping filter, and this class exists so the two can be compared
  under the same startup, magnetic and adaptation machinery.

  Everything except the vertical motion model is the OU-III wrapper's
  schedule, reused rather than re-tuned:

    * the shared measurement-only front end (Mahony proxy, wave-period
      estimator, sigma band, auto-tuner, wave direction), the vibration
      conditioning and the proxy-startup pipeline (SeaStateFusion_CARMA);
    * the per-sample order of updateCore_(), the one-sample staged commit of
      the operating point, the tilt watchdog, the magnetometer-service gate on
      accelerometer-bias learning, the vibration Racc inflation;
    * the horizontal OU chain with the deployed SpectralMSE r_S law, its x/y
      factors, the tau-scaled S cadence and the periodic a_w covariance floor
      (horizontal axes).

  Vertical operating point.  The same measured (tau, sigma_aw) that tune
  OU-III scale the fixed CARMA spectral shape:

      kappa    = f_tune / f_z,model,   f_tune = c_tau / (2 tau)
                 (the model's displacement zero-crossing frequency is placed on
                 the measured one; WavePeriodEstimator measures T_z),
      sigma_a  = c_carma * max(sigma_floor, sigma_aw / c_sigma)
                 (the CARMA stationary vertical-acceleration std is the
                 band-limited measured acceleration RMS sigma_a,B).

  setFixedVerticalCarma() freezes (kappa, sigma_a) for the fixed-parameter
  evaluation while the horizontal chain keeps adapting.

  With the S_z variant (with_vertical_integral) the S_z = 0 update uses the
  pole-placement law derived in docs/carma-heave-filter.md: near DC the CARMA
  displacement is white with level c0 = S_h(0), so S_z is a random walk and
  the S_z = 0 pseudo-measurement every T_S closes a first-order loop at
  omega_R = sqrt(c0 / (r_Sz^2 T_S)).  Placing omega_R tau = kappa_R (the
  OU-III pole constant) gives r_Sz = sqrt(c0) tau / (kappa_R sqrt(T_S)).
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numbers>

#include "kalman_carma/Kalman3D_Wave_CARMA.h"
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"

namespace ocean_imu::carma {

template<TrackerType trackerT, int NCARMA = 5, bool with_vertical_integral = false>
class SeaStateFusionFilter_CARMA
    : public seastate::common::AccelVibrationConditioning,
      public seastate::common::MarineWaveFrontEnd<trackerT>,
      public seastate::common::SeaStateAdaptationCommon
{
    using FrontEnd   = seastate::common::MarineWaveFrontEnd<trackerT>;
    using Adaptation = seastate::common::SeaStateAdaptationCommon;
    using FrontEnd::levelVertical_;
    using FrontEnd::trackAccelBandFrequency_;
    using FrontEnd::updateWavePeriodAndDirection_;
    using FrontEnd::tuner_frequency_hz_;
    using FrontEnd::frontEndStill_;
    using FrontEnd::frontEndStillTimeSec_;
    using FrontEnd::resetFrontEnd_;
    using FrontEnd::gravity_mps2_;
    using FrontEnd::low_wave_noise_;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Mekf = kalman::Kalman3D_Wave_CARMA<float, NCARMA, with_vertical_integral>;
    using FrontEnd::startupProxyQuat;
    using FrontEnd::wavePeriodUsable;

    explicit SeaStateFusionFilter_CARMA(bool with_mag = true)
        : Adaptation(Adaptation::EstimatorBounds{
              ou3::MAX_TUNE_FREQ_HZ, ou3::MIN_TAU_S, ou3::MAX_TAU_S, ou3::MAX_SIGMA_A,
              ou3::PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT}),
          with_mag_(with_mag)
    {
        setAccelVibrationGuard(ACC_VIBRATION_GUARD_HZ_DEFAULT,
                               ACC_VIBRATION_GUARD_POLES_DEFAULT);
        setAccelVibrationRaccGain(ACC_VIBRATION_RACC_GAIN_DEFAULT);
    }

    void initialize(const Eigen::Vector3f& sigma_a, const Eigen::Vector3f& sigma_g,
                    const Eigen::Vector3f& sigma_m) {
        mekf_ = std::make_unique<Mekf>(sigma_a, sigma_g, sigma_m);
        seastate::common::finalizeInitialization(
            mekf_, [this]() { enterCold_(); }, [this]() { apply_ou_tune_(true); });
    }

    void initialize_ext(const Eigen::Vector3f& sigma_a, const Eigen::Vector3f& sigma_g,
                        const Eigen::Vector3f& sigma_m, float Pq0, float Pb0, float b0,
                        float R_S_noise, float gravity_magnitude) {
        gravity_mps2_ = gravity_magnitude;
        mekf_ = std::make_unique<Mekf>(sigma_a, sigma_g, sigma_m, Pq0, Pb0, b0, R_S_noise,
                                       gravity_magnitude);
        seastate::common::finalizeInitialization(
            mekf_, [this]() { enterCold_(); }, [this]() { apply_ou_tune_(true); });
    }

    bool initialize_from_acc(const Eigen::Vector3f& acc_body_ned) {
        return mekf_ ? mekf_->initialize_from_acc(acc_body_ned) : false;
    }

    void updateTime(float dt, const Eigen::Vector3f& gyro, const Eigen::Vector3f& acc,
                    float tempC = 35.0f) {
        updateCore_(dt, gyro, acc, tempC, /*drive_mekf=*/true);
    }

    void updateFrontEnd(float dt, const Eigen::Vector3f& gyro, const Eigen::Vector3f& acc) {
        updateCore_(dt, gyro, acc, /*tempC=*/35.0f, /*drive_mekf=*/false);
    }

    void goLive(const Eigen::Quaternionf& q_bw, float tilt_sigma_rad, float yaw_sigma_rad,
                bool allow_acc_bias = false) {
        if (!mekf_) return;
        if (!q_bw.coeffs().allFinite() || !(q_bw.norm() > 1e-8f)) return;
        mekf_->initialize_from_attitude(q_bw, tilt_sigma_rad, yaw_sigma_rad);
        if (allow_acc_bias) accel_bias_locked_ = false;
        enterLive_();
    }

private:
    // OU-III's per-sample schedule, unchanged.
    void updateCore_(float dt, const Eigen::Vector3f& gyro, const Eigen::Vector3f& acc,
                     float tempC, bool drive_mekf) {
        if (!mekf_) return;
        if (!(dt > 0.0f) || !std::isfinite(dt)) return;
        apply_pending_online_tune_();
        time_ += dt;
        startup_stage_t_ += dt;

        const Eigen::Vector3f acc_in = conditionAccel_(acc, dt);
        levelVertical_(dt, gyro, acc_in);
        if (drive_mekf) apply_racc_vibration_inflation_();
        if (drive_mekf) {
            mekf_->time_update(gyro, dt);
            mekf_->measurement_update_acc_only(acc_in, tempC);
        }
        if (drive_mekf && tilt_watchdog_.step(mekf_->quaternion_boat(), dt)) {
            if (startup_stage_ == StartupStage::Live) {
                mekf_->initialize_from_acc_preserve_yaw(acc_in);
            } else {
                mekf_->initialize_from_acc(acc_in);
                enterCold_();
                resetTrackingState_();
            }
        }
        const float a_vert_measurement = trackAccelBandFrequency_(dt);
        if (enable_tuner_) update_tuner(dt, a_vert_measurement, tuner_frequency_hz_());
        if (drive_mekf) periodic_aw_cov_sync_tick_();
        if (drive_mekf) {
            const Eigen::Vector3f acc_bias = mekf_->get_acc_bias_body_at_temperature(tempC);
            updateWavePeriodAndDirection_(dt, mekf_->quaternion_boat(), acc_in, &acc_bias);
        } else {
            updateWavePeriodAndDirection_(dt, startupProxyQuat(), acc_in, nullptr);
        }
    }

public:
    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!with_mag_ || !mekf_) return;
        if (time_ < mag_delay_sec_) return;
        mekf_->measurement_update_mag_only(mag_body_ned);
        const float t = static_cast<float>(time_);
        mag_service_.record(mekf_->lastMagDiag().accepted, t);
        if (accel_bias_locked_ && startup_stage_ == StartupStage::Live &&
            mag_service_.sufficient(mag_updates_to_unlock_, t)) {
            accel_bias_locked_ = false;
            if (!acc_bias_hold_) mekf_->set_acc_bias_updates_enabled(true);
        }
    }

    void setWithMag(bool with_mag) { with_mag_ = with_mag; }
    void setSFactor(float s) { if (std::isfinite(s) && s > 0.0f) S_factor_ = s; }
    void setRSXFactor(float k) { if (std::isfinite(k)) R_S_x_factor_ = std::min(std::max(k, 0.0f), 4.0f); }
    void setRSYFactor(float k) { if (std::isfinite(k)) R_S_y_factor_ = std::min(std::max(k, 0.0f), 4.0f); }
    void setTauCoeff(float c) { if (std::isfinite(c) && c > 0.0f) tau_coeff_ = c; }
    void setSigmaCoeff(float c) { if (std::isfinite(c) && c > 0.0f) sigma_coeff_ = c; }
    // Multiplier c_carma on the measured band acceleration RMS.
    void setCarmaSigmaGain(float c) { if (std::isfinite(c) && c > 0.0f) carma_sigma_gain_ = c; }
    // Absolute floor on the CARMA acceleration std, combined (max) with the
    // band-referred accelerometer noise floor.  The default 0.05 m/s^2 is the
    // floor OU-III applies to its a_w prior; a smaller value is the mechanism
    // ablation of docs/carma-heave-filter.md, not a deployed setting.
    void setCarmaSigmaFloor(float f) { if (std::isfinite(f) && f >= 0.0f) carma_sigma_floor_ = f; }
    void setMagDelaySec(float d) { if (std::isfinite(d) && d >= 0.0f) mag_delay_sec_ = d; }
    void setMagUpdatesToUnlockAccBias(int n) { if (n >= 0) mag_updates_to_unlock_ = n; }
    void setNominalRaccStd(const Eigen::Vector3f& r) { Racc_nominal_ = r; }

    void setAccBiasHold(bool hold) {
        if (acc_bias_hold_ == hold) return;
        acc_bias_hold_ = hold;
        if (!mekf_) return;
        if (hold) {
            mekf_->set_acc_bias_updates_enabled(false);
            return;
        }
        if (!accel_bias_locked_ && startup_stage_ == StartupStage::Live)
            mekf_->set_acc_bias_updates_enabled(true);
    }
    bool accBiasHeld() const noexcept { return acc_bias_hold_; }

    // Freeze the whole operating point (as OU-III setFixedTuning): tau and
    // sigma_aw drive the horizontal chain and, through the same mapping, the
    // vertical CARMA scale and amplitude.
    bool setFixedTuning(float tau_s, float sigma_a, float RS) {
        if (!(std::isfinite(tau_s) && tau_s > 0.0f && std::isfinite(sigma_a) && sigma_a > 0.0f &&
              std::isfinite(RS) && RS > 0.0f)) return false;
        enable_tuner_ = false;
        online_tune_apply_pending_ = false;
        tau_target_ = std::min(std::max(tau_s, min_tau_s_), max_tau_s_);
        sigma_target_ = std::min(sigma_a, max_sigma_a_);
        RS_target_ = std::min(std::max(RS, min_R_S_), max_R_S_);
        tune_.tau_applied = tau_target_;
        tune_.sigma_applied = sigma_target_;
        tune_.RS_applied = RS_target_;
        apply_ou_tune_(true);
        if (startup_stage_ == StartupStage::Live) apply_RS_tune_();
        return true;
    }

    // Hold the vertical frequency scale at kappa (kappa = 1 is the 0.29 Hz
    // design shape) and, when accel_std > 0, the amplitude too; with
    // accel_std <= 0 the amplitude keeps following the measured sea.  The
    // horizontal chain and the tuner keep running.  kappa <= 0 releases the
    // hold.
    bool setFixedVerticalCarma(float kappa, float accel_std) {
        if (!(std::isfinite(kappa) && kappa > 0.0f)) {
            vertical_fixed_ = false;
            return true;
        }
        vertical_fixed_ = true;
        fixed_kappa_ = kappa;
        fixed_accel_std_ = (std::isfinite(accel_std) && accel_std > 0.0f) ? accel_std : 0.0f;
        if (mekf_) apply_vertical_tune_();
        return true;
    }

    inline float getTauApplied()   const noexcept { return mekf_ ? mekf_->get_aw_time_constant() : NAN; }
    inline float getSigmaApplied() const noexcept { return tune_.sigma_applied; }
    inline float getRSApplied()    const noexcept { return mekf_ ? mekf_->get_RS_noise_std().x() : NAN; }
    inline float getRSTarget()     const noexcept { return RS_target_; }
    inline float getVerticalKappa() const noexcept { return mekf_ ? mekf_->vertical_model().kappa() : NAN; }
    inline float getVerticalAccelStd() const noexcept {
        return mekf_ ? mekf_->vertical_model().accel_std() : NAN;
    }
    inline float getVerticalDisplacementStd() const noexcept {
        return mekf_ ? mekf_->vertical_model().displacement_std() : NAN;
    }

    inline float getDisplacementScale(bool smoothed = true) const noexcept {
        const float tau = smoothed ? tune_.tau_applied : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;
        if (!std::isfinite(sigma) || !std::isfinite(tau)) return NAN;
        constexpr float C_HS = 2.0f * std::numbers::sqrt2_v<float> /
                               (std::numbers::pi_v<float> * std::numbers::pi_v<float>);
        return C_HS * sigma * tau * tau / 2.0f;
    }
    float getVerticalSpeedEnvelopeMps(bool smoothed = true) const noexcept {
        const float tau = smoothed ? tune_.tau_applied : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;
        if (!(tau > 1e-6f) || !std::isfinite(tau) || !std::isfinite(sigma)) return NAN;
        return std::numbers::sqrt2_v<float> / std::numbers::pi_v<float> * sigma * tau;
    }

    inline Mekf& mekf() noexcept { return *mekf_; }
    inline const Mekf& mekf() const noexcept { return *mekf_; }

private:
    void apply_pending_online_tune_() {
        if (!online_tune_apply_pending_ || !mekf_) return;
        apply_ou_tune_(false);
        if (startup_stage_ == StartupStage::Live) apply_RS_tune_();
        online_tune_apply_pending_ = false;
    }

    void apply_ou_tune_(bool sync_covariance) {
        if (!mekf_) return;
        mekf_->set_aw_time_constant(tune_.tau_applied);
        apply_pseudo_update_cadence_(mekf_.get(), tune_.tau_applied);
        const float sigma_floor = std::max(0.05f, band_noise_floor_sigma_());
        const float sZ = std::max(sigma_floor, tune_.sigma_applied);
        const float sH = sZ * S_factor_;
        mekf_->set_aw_stationary_std(Eigen::Vector3f(sH, sH, sZ));
        apply_vertical_tune_();
        if (sync_covariance) {
            mekf_->synchronize_aw_covariance_to_stationary();
            last_aw_cov_sync_sec_ = time_;
        }
    }

    // Vertical CARMA scale and amplitude from the measured operating point.
    void apply_vertical_tune_() {
        if (!mekf_) return;
        const float c_sigma = (std::isfinite(sigma_coeff_) && sigma_coeff_ > 0.0f) ? sigma_coeff_ : 1.0f;
        const float sigma_floor = std::max({carma_sigma_floor_, band_noise_floor_sigma_(), 1e-4f});
        const float sigma_aB = std::max(sigma_floor, tune_.sigma_applied / c_sigma);
        const float accel_std = (vertical_fixed_ && fixed_accel_std_ > 0.0f)
                                ? fixed_accel_std_ : carma_sigma_gain_ * sigma_aB;
        if (vertical_fixed_) {
            mekf_->set_vertical_carma(fixed_kappa_, accel_std);
            return;
        }
        const float tau = tune_.tau_applied;
        if (!(std::isfinite(tau) && tau > 0.0f)) return;
        const double fz_base = mekf_->vertical_model().base_zero_crossing_hz();
        const float f_tune = tau_coeff_ * 0.5f / tau;
        const float kappa = static_cast<float>(f_tune / fz_base);
        mekf_->set_vertical_carma(kappa, accel_std);
    }

    void periodic_aw_cov_sync_tick_() {
        if (!mekf_ || !periodicAwCovSyncDue_()) return;
        mekf_->synchronize_aw_covariance_to_stationary();
        last_aw_cov_sync_sec_ = time_;
    }

    // OU-III SpectralMSE law for the horizontal S updates (its constants).
    float rs_spectral_mse_target_(float tau, float sigma) const noexcept {
        const float TS = pseudo_update_period_for_(tau);
        const float c_sigma = (std::isfinite(sigma_coeff_) && sigma_coeff_ > 0.0f) ? sigma_coeff_ : 1.0f;
        const float sigma_aB = std::max(sigma / c_sigma, 1e-6f);
        if (!(TS > 0.0f) || !(tau > 0.0f)) return ou3::R_S_MSE_COEFF_DEFAULT;
        const float tau2 = tau * tau;
        return ou3::R_S_MSE_COEFF_DEFAULT * rs_qeff_pow_ *
               std::pow(sigma_aB * tau2 * tau2, 6.0f / 7.0f) / std::sqrt(TS);
    }

    // S_z law of the integral variant: r_Sz = sqrt(c0) tau / (kappa_R sqrt(T_S)).
    float rs_vertical_target_() const noexcept {
        if (!mekf_) return NAN;
        const auto& vm = mekf_->vertical_model();
        const double c0 = static_cast<double>(vm.noise_gain2()) *
                          kalman::carma::carma_unit_psd<NCARMA>(vm.shape(), vm.kappa(), 0.0);
        const float tau = tune_.tau_applied;
        const float TS = mekf_->get_pseudo_update_period_s();
        if (!(c0 > 0.0) || !(tau > 0.0f) || !(TS > 0.0f)) return NAN;
        return static_cast<float>(std::sqrt(c0)) * tau /
               (ou3::R_S_POLE_KAPPA_DEFAULT * std::sqrt(TS));
    }

    void apply_RS_tune_() {
        if (!mekf_) return;
        const float rs = std::min(std::max(tune_.RS_applied, min_R_S_), max_R_S_);
        float rs_z = rs;
        if constexpr (with_vertical_integral) {
            const float v = rs_vertical_target_();
            if (std::isfinite(v) && v > 0.0f) rs_z = v;
        }
        mekf_->set_RS_noise(Eigen::Vector3f(rs * R_S_x_factor_, rs * R_S_y_factor_, rs_z));
    }

    void update_tuner(float dt, float a_vertical_measurement, float freq_hz_for_tuner) {
        stepVarianceChannel_(dt, a_vertical_measurement, freq_hz_for_tuner);
        if (!advanceStartupStage_(wavePeriodUsable())) return;
        constexpr float STILL_VAR_DECAY_SEC = 1.0f;
        const float sea_time_sec = measureOperatingPoint_(
            frontEndStill_(), frontEndStillTimeSec_(), STILL_VAR_DECAY_SEC, tau_coeff_, sigma_coeff_);
        const float RS_raw = rs_spectral_mse_target_(tau_target_, sigma_target_);
        RS_target_ = enable_clamp_ ? std::min(std::max(RS_raw, min_R_S_), max_R_S_) : RS_raw;
        adapt_mekf(dt, tau_target_, sigma_target_, RS_target_, sea_time_sec);
    }

    void adapt_mekf(float dt, float tau_t, float sigma_t, float RS_t, float sea_time_sec) {
        const float alpha = tauSigmaSmoothingAlpha_(dt, sea_time_sec);
        const float RS_sec = seastate::common::adaptiveSmoothingHorizonSec(
            ou3::ADAPT_RS_MULT, tau_t, RS_t, tune_.RS_applied, ou3::ADAPT_RS_SLEW_LOG, dt);
        const float alpha_RS = 1.0f - std::exp(-dt / RS_sec);
        tune_.tau_applied   += alpha    * (tau_t   - tune_.tau_applied);
        tune_.sigma_applied += alpha    * (sigma_t - tune_.sigma_applied);
        tune_.RS_applied    += alpha_RS * (RS_t    - tune_.RS_applied);
        stageOnlineTuneCommit_();
    }

    void resetTrackingState_() {
        resetFrontEnd_();
        resetAdaptation_();
    }

    void enterCold_() {
        startup_stage_ = StartupStage::Cold;
        startup_stage_t_ = 0.0f;
        if (!mekf_) return;
        accel_bias_locked_ = with_mag_;
        mag_service_.reset();
        mekf_->set_acc_bias_updates_enabled(false);
    }

    Eigen::Vector3f racc_base_std_() const {
        if (Racc_nominal_.allFinite() && Racc_nominal_.minCoeff() > 0.0f) return Racc_nominal_;
        return Eigen::Vector3f::Zero();
    }

    void apply_racc_vibration_inflation_() {
        if (!mekf_ || (!(racc_vibration_gain_ > 0.0f) &&
                       !(low_wave_noise_.max_std_scale > 1.0f) && !racc_inflated_)) return;
        const Eigen::Vector3f base = racc_base_std_();
        if (!(base.minCoeff() > 0.0f)) return;
        const Eigen::Vector3f scales = isAdaptiveLive()
            ? low_wave_noise_.scales(tune_.sigma_applied / sigma_coeff_, tuner_frequency_hz_(), base)
            : Eigen::Vector3f::Ones();
        commandRaccStd_(*mekf_, base, scales);
    }

    void enterLive_() {
        startup_stage_ = StartupStage::Live;
        startup_stage_t_ = 0.0f;
        if (!mekf_) return;
        apply_ou_tune_(true);
        // Seat the horizontal a_w marginal and the vertical CARMA block on the
        // committed operating point before the first prediction.
        mekf_->reset_aw_covariance_to_stationary();
        mekf_->set_acc_bias_updates_enabled(!accel_bias_locked_ && !acc_bias_hold_);
        apply_RS_tune_();
    }

    Eigen::Vector3f Racc_nominal_ = Eigen::Vector3f::Constant(0.0f);
    bool accel_bias_locked_ = true;
    seastate::common::MagServiceCount mag_service_{};
    int  mag_updates_to_unlock_ = 250;
    bool acc_bias_hold_ = false;
    bool with_mag_;
    float mag_delay_sec_ = MAG_DELAY_SEC;
    seastate::common::TiltResetWatchdog tilt_watchdog_{};

    float rs_qeff_pow_ = std::pow(2.0f * ou3::R_S_ACCEL_NOISE_DENSITY_DEFAULT, 1.0f / 14.0f);
    float min_R_S_ = ou3::MIN_R_S;
    float max_R_S_ = ou3::MAX_R_S;
    float R_S_x_factor_ = 0.72f;
    float R_S_y_factor_ = 0.50f;
    float S_factor_ = 1.0f;

    ou3::TuneState tune_;
    float RS_target_ = NAN;
    float tau_coeff_ = 1.0f;
    float sigma_coeff_ = 0.9f;
    float carma_sigma_gain_ = 1.0f;
    float carma_sigma_floor_ = 0.05f;

    bool vertical_fixed_ = false;
    float fixed_kappa_ = 1.0f;
    float fixed_accel_std_ = 1.0f;

    std::unique_ptr<Mekf> mekf_;
};

struct SeaStateFusionConfig_CARMA : seastate::common::ProxyStartupFusionConfig {
    float R_S_noise = 1.5f;
};

template<TrackerType trackerT, int NCARMA = 5, bool with_vertical_integral = false>
class SeaStateFusion_CARMA
    : public seastate::common::ProxyStartupFusion<
          SeaStateFusionFilter_CARMA<trackerT, NCARMA, with_vertical_integral>,
          SeaStateFusionConfig_CARMA>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using Config = SeaStateFusionConfig_CARMA;
    using Impl = SeaStateFusionFilter_CARMA<trackerT, NCARMA, with_vertical_integral>;

    void begin(const Config& cfg) {
        this->beginStartup_(cfg, [](Impl& impl, const Config& c) {
            impl.initialize_ext(c.sigma_a, c.sigma_g, c.sigma_m, c.Pq0, c.Pb0, c.b0, c.R_S_noise,
                                c.gravity_magnitude);
        });
    }
};

}  // namespace ocean_imu::carma
