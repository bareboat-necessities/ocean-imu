#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy
  Released under the MIT License

  SeaStateFusionFilter_CARMA3 -- EXPERIMENTAL research wrapper around
  Kalman3D_Wave_CARMA3 (three-axis CARMA translational prior).  Not a deployed
  estimator; SeaStateFusionFilter_OU_III stays the shipping filter.  It keeps
  OU-III's per-sample schedule (shared front end, proxy startup, staged tune
  commit, tilt watchdog, magnetometer-service bias gate, vibration Racc
  inflation) and OU-III's tuning envelope and r_S law
  (kalman_common/SeaStateTuning.h), and maps the measured operating point
  (tau, sigma_aw) onto the CARMA priors:

    kappa_z  = f_tune / f_z,z-model,      f_tune = c_tau / (2 tau)
               (the model's displacement zero-crossing frequency is placed
               on the measured one)
    sigma_z  = c_carma * max(sigma_floor, sigma_aw / c_sigma)
               (the measured band acceleration RMS)
    kappa_h  = f_tune / f_z,h-model       (the same measured wave frequency
                                           that scales the vertical shape)
    sigma_s  = S_x,y * c_carma * max(sigma_floor, sigma_aw / c_sigma)
               (the vertical band acceleration RMS times a keel-frame factor:
               x along the keel (surge), y across it (sway); S_x = S_y = 1 is
               OU-III's isotropic horizontal prior sigma_aw,H = sigma_aw,Z)
    r        = 0 (independent surge and sway) unless setHorizontalCorrelation().

  The CARMA horizontal states are world-frame (NED).  With S_x = S_y and r = 0
  the horizontal prior is isotropic and yaw-invariant.  Otherwise the
  keel-frame covariance sigma_z^2 [S_x^2, r S_x S_y; r S_x S_y, S_y^2] is
  rotated into NED by the inner filter's own heading every time the tune is
  committed, so the prior turns with the hull and the filter stays
  yaw-equivariant (docs/tfg-carma3-research.md section 2.10): a prior fixed in
  world coordinates would carry heading information.

  No horizontal acceleration statistic is measured: the mapping is fixed and
  uses only quantities the shipping tuner already produces.

  Wave excitation (setCarmaExcitationMode, kalman_carma/CarmaWaveEvidence.h,
  docs/carma-wave-excitation.md).  The formula above is the Legacy mode.  Its
  0.05 m/s^2 floor, and the legacy sigma channel under it, are not wave
  evidence: the channel subtracts the configured 0.12 m/s^2 pre-band floor
  (referred through the band, 35 times the band variance the simulated
  1.5 mg sensor produces in still water) and
  attenuates the variance whenever the stillness detector's fixed
  0.28 m/s^2 RMS threshold reports still water, so in still water and on
  every small sea the CARMA driving noise is the floor itself.  The Evidence
  mode replaces sigma_aw / c_sigma by sqrt(p max(E, 0)), E the band variance
  in excess of the sensor-noise model and p its detection weight, and moves
  the CARMA frequency only with p times the confidence of the measured
  period (the acceleration band confirms where the waves are), holding it
  otherwise.  Debiased and FloorOff are its ablations.
  The only floor left, setCarmaNumericalSigmaFloor (1e-5 m/s^2), is
  numerical: it keeps the configured models valid and the covariance out of
  float subnormals.

  MekfT is the inner filter.  It defaults to Kalman3D_Wave_CARMA3; the
  two-frame Lie-group variant of SeaStateFusionFilter_TFG_CARMA3.h runs through
  this same schedule so the two differ only in the estimator.

  IMU lever arm (setImuLeverArmBody, off by default).  An IMU at body offset r
  from the centre of gravity measures alpha x r + omega x (omega x r) on top of
  the CoG specific force.  The wrapper removes it at its input, from the
  measured rate (kalman_common/ImuLeverArm.h), so the vibration guard, the
  startup proxy and the magnetic reference it learns, the tuner, the direction
  estimator and the inner filter all see CoG specific force.  The inner
  filter's own set_imu_lever_arm_body() (OU-III's placement) corrects only its
  accelerometer update; docs/carma-imu-lever-arm.md measures the difference.
  Use one or the other, never both.

  Horizontal S law.  For a horizontal CARMA with S (SH):
    * if the shape has a nonzero DC level c0 = S_p(0) (e.g. CARMA(5,2) as
      fitted), S is a random walk under the prior and S = 0 every T_S closes a
      first-order loop at omega_R = sqrt(c0/(r^2 T_S)); omega_R tau = kappa_R
      (OU-III's pole constant) gives r_S = sqrt(c0) tau / (kappa_R sqrt(T_S))
      per axis, and the same law sets r_Sz for S_z (SZ);
    * if the shape has a zero at DC (c0 = 0), the S mode is uncontrollable
      and no pole-placement law exists; the OU-III SpectralMSE r_S and its
      x/y factors are then carried over unchanged.  This is the negative
      control of docs/carma-three-axis.md (not a derived law).
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

#include "kalman_carma/CarmaWaveEvidence.h"
#include "kalman_carma/Kalman3D_Wave_CARMA3.h"
#include "kalman_common/AccelVibrationConditioning.h"
#include "kalman_common/ImuLeverArm.h"
#include "kalman_common/MarineWaveFrontEnd.h"
#include "kalman_common/ProxyStartupFusion.h"
#include "kalman_common/SeaStateAdaptationCommon.h"
#include "kalman_common/SeaStateFusionDefaults.h"
#include "kalman_common/SeaStateFusionFilterCommon.h"
#include "kalman_common/SeaStateTuning.h"

namespace ocean_imu::carma {

namespace defaults = ::seastate::common::defaults;
namespace tuning = ::seastate::common::tuning;

template<TrackerType trackerT, int NHC, bool SH, int NZC, bool SZ,
         class MekfT = kalman::Kalman3D_Wave_CARMA3<float, NHC, SH, NZC, SZ>>
class SeaStateFusionFilter_CARMA3
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

    using Mekf = MekfT;
    using FrontEnd::startupProxyQuat;
    using FrontEnd::wavePeriodUsable;

    explicit SeaStateFusionFilter_CARMA3(bool with_mag = true)
        : Adaptation(Adaptation::EstimatorBounds{
              tuning::MAX_TUNE_FREQ_HZ, tuning::MIN_TAU_S, tuning::MAX_TAU_S, tuning::MAX_SIGMA_A,
              tuning::PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT}),
          with_mag_(with_mag)
    {
        setAccelVibrationGuard(defaults::ACC_VIBRATION_GUARD_HZ, defaults::ACC_VIBRATION_GUARD_POLES);
        setAccelVibrationRaccGain(defaults::ACC_VIBRATION_RACC_GAIN);
    }

    void initialize(const Eigen::Vector3f& sigma_a, const Eigen::Vector3f& sigma_g,
                    const Eigen::Vector3f& sigma_m) {
        lever_.reset_kinematics();
        this->resetFrontendGyroHistory();
        mekf_ = std::make_unique<Mekf>(sigma_a, sigma_g, sigma_m);
        applyGyroSampling_();
        seastate::common::finalizeInitialization(
            mekf_, [this]() { enterCold_(); }, [this]() { apply_tune_(); });
    }

    void initialize_ext(const Eigen::Vector3f& sigma_a, const Eigen::Vector3f& sigma_g,
                        const Eigen::Vector3f& sigma_m, float Pq0, float Pb0, float b0,
                        float R_S_noise, float gravity_magnitude) {
        gravity_mps2_ = gravity_magnitude;
        lever_.reset_kinematics();  // a restart may follow a gap in the samples
        this->resetFrontendGyroHistory();
        mekf_ = std::make_unique<Mekf>(sigma_a, sigma_g, sigma_m, Pq0, Pb0, b0, R_S_noise,
                                       gravity_magnitude);
        applyGyroSampling_();
        seastate::common::finalizeInitialization(
            mekf_, [this]() { enterCold_(); }, [this]() { apply_tune_(); });
    }

    bool initialize_from_acc(const Eigen::Vector3f& acc_body_ned) {
        return mekf_ ? mekf_->initialize_from_acc(acc_body_ned) : false;
    }

    // Independent core/frontend policies survive raw-core replacement. The
    // live wrapper chooses the validated core-only default.
    void setGyroSampleMode(kalman::filter_detail::GyroSampleMode core,
                           kalman::filter_detail::GyroSampleMode frontend) {
        core_gyro_sample_mode_ = core;
        if (mekf_) mekf_->set_gyro_sample_mode(core);
        this->setFrontendGyroSampleMode(frontend);
    }
    void setGyroSamplePeriod(float period) {
        core_gyro_sample_period_ = std::isfinite(period) && period > 0.0f ? period : 0.0f;
        if (mekf_) mekf_->set_gyro_sample_period(core_gyro_sample_period_);
    }
    // A diagnostic freeze preserves live state, covariance and all learned
    // statistics. The already staged commit is consumed once, then held.
    void setSeaStateAdaptationEnabled(bool enabled) { enable_tuner_ = enabled; }

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
    // The measured specific force moved from the IMU to the CoG.
    Eigen::Vector3f cogSpecificForce_(const Eigen::Vector3f& acc, const Eigen::Vector3f& gyro, float dt) {
        lever_.update_rate(gyro, dt);
        return acc - lever_.acceleration();
    }

    void updateCore_(float dt, const Eigen::Vector3f& gyro, const Eigen::Vector3f& acc,
                     float tempC, bool drive_mekf) {
        if (!mekf_) return;
        // A non-finite sample is dropped whole: the frontend, tuner and
        // trackers share it with the core.
        if (!(dt > 0.0f) || !std::isfinite(dt) || !gyro.allFinite() || !acc.allFinite()) return;
        apply_pending_online_tune_();
        time_ += dt;
        startup_stage_t_ += dt;

        const Eigen::Vector3f acc_in = lever_.enabled() ? conditionAccel_(cogSpecificForce_(acc, gyro, dt), dt)
                                                        : conditionAccel_(acc, dt);
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

    // IMU installation lever arm: r_b is the IMU position relative to the
    // vessel centre of gravity in the body frame (x forward, y starboard,
    // z down), metres.  Zero, the default, switches the correction off.  See
    // the file header; do not also set the inner filter's lever arm.
    void setImuLeverArmBody(const Eigen::Vector3f& r_b) {
        const bool was_enabled = lever_.enabled();
        lever_.set_offset_body(r_b);
        if (lever_.enabled() && !was_enabled) lever_.reset_kinematics();
    }
    void clearImuLeverArm() { lever_.clear(); }
    const Eigen::Vector3f& imuLeverArmBody() const noexcept { return lever_.offset_body(); }
    // Corner of the low-pass ahead of the rate derivative (default 15 Hz).
    void setImuLeverArmCutoffHz(float hz) { lever_.set_derivative_cutoff_hz(hz); }
    float imuLeverArmCutoffHz() const noexcept { return lever_.derivative_cutoff_hz(); }
    // Horizontal amplitude factors in the keel frame (x along the keel, y
    // across it), relative to the vertical; setSFactor sets both.
    void setSFactor(float s) { setSFactorXY(s, s); }
    void setSFactorXY(float sx, float sy) {
        if (!(std::isfinite(sx) && sx > 0.0f && std::isfinite(sy) && sy > 0.0f)) return;
        Sx_factor_ = sx;
        Sy_factor_ = sy;
    }
    float getSxFactor() const noexcept { return Sx_factor_; }
    float getSyFactor() const noexcept { return Sy_factor_; }
    void setRSXFactor(float k) { if (std::isfinite(k)) R_S_x_factor_ = std::min(std::max(k, 0.0f), 4.0f); }
    void setRSYFactor(float k) { if (std::isfinite(k)) R_S_y_factor_ = std::min(std::max(k, 0.0f), 4.0f); }
    void setTauCoeff(float c) { if (std::isfinite(c) && c > 0.0f) tau_coeff_ = c; }
    void setSigmaCoeff(float c) { if (std::isfinite(c) && c > 0.0f) sigma_coeff_ = c; }
    void setCarmaSigmaGain(float c) { if (std::isfinite(c) && c > 0.0f) carma_sigma_gain_ = c; }
    void setCarmaSigmaFloor(float f) { if (std::isfinite(f) && f >= 0.0f) carma_sigma_floor_ = f; }

    // Wave-excitation policy (file header).  Takes effect at the next commit.
    void setCarmaExcitationMode(CarmaExcitationMode m) { excitation_mode_ = m; }
    CarmaExcitationMode carmaExcitationMode() const noexcept { return excitation_mode_; }
    // Evidence-mode ablations: drop the evidence weight p from the amplitude,
    // or let the frequency follow the tuner instead of being held.
    void setCarmaEvidenceGate(bool on) { evidence_gate_ = on; }
    void setCarmaFrequencyHold(bool on) { frequency_hold_ = on; }
    // Accelerometer noise model of the wave-evidence statistic (vertical,
    // after levelling): white per-sample std at the IMU rate, bias random
    // walk, and the relative uncertainty of the modelled noise variance.
    void setCarmaSensorNoise(float white_sigma_mps2, float bias_rw_mps2_sqrts, float rel_var_uncertainty) {
        if (!(std::isfinite(white_sigma_mps2) && white_sigma_mps2 >= 0.0f)) return;
        if (!(std::isfinite(bias_rw_mps2_sqrts) && bias_rw_mps2_sqrts >= 0.0f)) return;
        if (!(std::isfinite(rel_var_uncertainty) && rel_var_uncertainty >= 0.0f)) return;
        sensor_noise_ = CarmaSensorNoise{white_sigma_mps2, bias_rw_mps2_sqrts, rel_var_uncertainty};
    }
    const CarmaSensorNoise& carmaSensorNoise() const noexcept { return sensor_noise_; }
    void setCarmaEvidenceThresholds(const CarmaEvidenceConfig& c) {
        if (!(std::isfinite(c.z0) && std::isfinite(c.z1) && c.z1 > c.z0)) return;
        if (!(std::isfinite(c.zf0) && std::isfinite(c.zf1) && c.zf1 > c.zf0)) return;
        if (!(std::isfinite(c.k0) && std::isfinite(c.k1) && c.k1 > c.k0)) return;
        evidence_cfg_ = c;
    }
    const CarmaEvidenceConfig& carmaEvidenceThresholds() const noexcept { return evidence_cfg_; }
    // Numerical lower bound of the CARMA acceleration std in the non-Legacy
    // modes; independent of the sensor and of the wave evidence.
    void setCarmaNumericalSigmaFloor(float s) { if (std::isfinite(s) && s > 0.0f) carma_numerical_floor_ = s; }
    float carmaNumericalSigmaFloor() const noexcept { return carma_numerical_floor_; }

    // Wave evidence of the last tuner sample (CarmaWaveEvidence.h).
    const CarmaWaveEvidence& waveEvidence() const noexcept { return evidence_; }
    // CARMA vertical acceleration std of the last committed prior.
    float getCarmaAccelStd() const noexcept { return carma_accel_std_applied_; }
    void setMagDelaySec(float d) { if (std::isfinite(d) && d >= 0.0f) mag_delay_sec_ = d; }
    void setMagUpdatesToUnlockAccBias(int n) { if (n >= 0) mag_updates_to_unlock_ = n; }
    void setNominalRaccStd(const Eigen::Vector3f& r) { Racc_nominal_ = r; }
    // Correlation coefficient of the surge and sway CARMA driving noise in the
    // keel frame (coupled prior); 0 = independent axes.
    void setHorizontalCorrelation(float r) { if (std::isfinite(r) && std::abs(r) < 1.0f) h_corr_ = r; }
    bool setHorizontalShape(const kalman::carma::CarmaShape& shape) {
        if (!mekf_ || !mekf_->set_horizontal_shape(shape)) return false;
        apply_tune_();
        if (startup_stage_ != StartupStage::Live) mekf_->reset_aw_covariance_to_stationary();
        return true;
    }
    bool setVerticalShape(const kalman::carma::CarmaShape& shape) {
        if (!mekf_ || !mekf_->set_vertical_shape(shape)) return false;
        apply_tune_();
        if (startup_stage_ != StartupStage::Live) mekf_->reset_aw_covariance_to_stationary();
        return true;
    }

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

    inline float getTauApplied()   const noexcept { return mekf_ ? mekf_->get_aw_time_constant() : NAN; }
    inline float getSigmaApplied() const noexcept { return tune_.sigma_applied; }
    inline float getRSApplied()    const noexcept { return mekf_ ? mekf_->get_RS_noise_std().x() : NAN; }
    inline float getRSTarget()     const noexcept { return RS_target_; }
    inline float getVerticalKappa() const noexcept { return mekf_ ? mekf_->vertical_model().kappa() : NAN; }
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
        apply_tune_();
        if (startup_stage_ == StartupStage::Live) apply_RS_tune_();
        online_tune_apply_pending_ = false;
    }

    float carma_accel_std_() const {
        const float c_sigma = (std::isfinite(sigma_coeff_) && sigma_coeff_ > 0.0f) ? sigma_coeff_ : 1.0f;
        if (excitation_mode_ != CarmaExcitationMode::Legacy) {
            // tune_.sigma_applied carries the mode's own target (update_tuner).
            return carma_sigma_gain_ * std::max(carma_numerical_floor_, tune_.sigma_applied / c_sigma);
        }
        const float sigma_floor = std::max({carma_sigma_floor_, band_noise_floor_sigma_(), 1e-4f});
        const float sigma_aB = std::max(sigma_floor, tune_.sigma_applied / c_sigma);
        return carma_sigma_gain_ * sigma_aB;
    }

    // Wave evidence from the levelled vertical acceleration at the band
    // corners and horizon the variance channel is in force with
    // (CarmaWaveEvidence.h).  Measurement-only, like the channel.
    void update_wave_evidence_(float dt, float a_vertical_measurement) {
        CarmaWaveEvidence::Inputs in;
        in.dt = dt;
        in.accel = a_vertical_measurement;
        in.band_low_hz = sigma_wave_band_.lowHz();
        in.band_high_hz = sigma_wave_band_.highHz();
        in.band_white_gain = sigma_wave_band_.whiteNoiseVarianceGain();
        in.var_horizon_s = tuner_.getVarianceHorizonSec();
        evidence_.update(in, sensor_noise_, evidence_cfg_);
    }

    // sigma_aw target of the non-Legacy modes, in the c_sigma units of the
    // legacy channel (carma_accel_std_ divides it out again).
    float excitation_sigma_target_() const {
        switch (excitation_mode_) {
            case CarmaExcitationMode::FloorOff:
                return sigma_target_;
            case CarmaExcitationMode::Debiased:
            case CarmaExcitationMode::Evidence: {
                const bool gated = excitation_mode_ == CarmaExcitationMode::Evidence && evidence_gate_;
                const float s = sigma_coeff_ * std::sqrt(evidence_.excitationVariance(gated));
                return enable_clamp_ ? std::min(s, max_sigma_a_) : s;
            }
            case CarmaExcitationMode::Legacy:
            default:
                return sigma_target_;
        }
    }

    // Weight of the tau (CARMA frequency) update.  The Evidence mode moves the
    // frequency only as far as the measured period is wave evidence and holds
    // the last committed value otherwise; before the period estimator is
    // usable the tuner runs on the fixed prior, which it follows.
    float frequency_update_weight_() const {
        if (excitation_mode_ != CarmaExcitationMode::Evidence || !frequency_hold_) return 1.0f;
        if (!wavePeriodUsable()) return 1.0f;
        return evidence_.waveProbability() * evidence_.frequencyConfidence();
    }

    void apply_tune_() {
        if (!mekf_) return;
        mekf_->set_aw_time_constant(tune_.tau_applied);
        apply_pseudo_update_cadence_(mekf_.get(), tune_.tau_applied);
        apply_carma_tune_();
    }

    // Vertical and horizontal CARMA scales and amplitudes from the measured
    // operating point (fixed mappings, see the file header).
    void apply_carma_tune_() {
        if (!mekf_) return;
        const float tau = tune_.tau_applied;
        if (!(std::isfinite(tau) && tau > 0.0f)) return;
        const float f_tune = tau_coeff_ * 0.5f / tau;
        const float accel_std = carma_accel_std_();
        carma_accel_std_applied_ = accel_std;
        const double fz_v = mekf_->vertical_model().base_zero_crossing_hz();
        mekf_->set_vertical_carma(static_cast<float>(f_tune / fz_v), accel_std);
        const double fz_h = mekf_->horizontal_model().base_zero_crossing_hz();
        const float kappa_h = static_cast<float>(f_tune / fz_h);
        const float sx = accel_std * Sx_factor_, sy = accel_std * Sy_factor_;
        if (sx == sy && h_corr_ == 0.0f) {
            // Isotropic: yaw-invariant, no heading enters.
            mekf_->set_horizontal_carma(kappa_h, sx, sx, 0.0f);
            return;
        }
        // Keel-frame covariance rotated into NED by the filter's heading:
        // the keel's horizontal direction is the first column of R(body->NED).
        const Eigen::Matrix3f R = mekf_->quaternion_boat().toRotationMatrix();
        float c = R(0, 0), s = R(1, 0);
        const float n = std::sqrt(c * c + s * s);
        if (n > 1e-6f && std::isfinite(n)) { c /= n; s /= n; } else { c = 1.0f; s = 0.0f; }
        const float a = sx * sx, b = sy * sy, m = h_corr_ * sx * sy;
        const float wxx = c * c * a - 2.0f * c * s * m + s * s * b;
        const float wyy = s * s * a + 2.0f * c * s * m + c * c * b;
        const float wxy = c * s * (a - b) + (c * c - s * s) * m;
        const float ox = std::sqrt(std::max(wxx, 0.0f)), oy = std::sqrt(std::max(wyy, 0.0f));
        const float r = (ox > 0.0f && oy > 0.0f) ? std::min(std::max(wxy / (ox * oy), -1.0f), 1.0f) : 0.0f;
        mekf_->set_horizontal_carma(kappa_h, ox, oy, r);
    }

    float rs_spectral_mse_target_(float tau, float sigma) const noexcept {
        const float TS = pseudo_update_period_for_(tau);
        const float c_sigma = (std::isfinite(sigma_coeff_) && sigma_coeff_ > 0.0f) ? sigma_coeff_ : 1.0f;
        const float sigma_aB = std::max(sigma / c_sigma, 1e-6f);
        if (!(TS > 0.0f) || !(tau > 0.0f)) return tuning::R_S_MSE_COEFF_DEFAULT;
        const float tau2 = tau * tau;
        return tuning::R_S_MSE_COEFF_DEFAULT * rs_qeff_pow_ *
               std::pow(sigma_aB * tau2 * tau2, 6.0f / 7.0f) / std::sqrt(TS);
    }

    // Pole-placement r_S for a CARMA axis with DC level c0 (unit gain c0_unit
    // times the axis noise intensity): r = sqrt(c0) tau / (kappa_R sqrt(T_S)).
    template<int N>
    float rs_pole_law_(const kalman::carma::CarmaVerticalModel<float,N>& m, float gain2) const noexcept {
        const double c0 = static_cast<double>(gain2) *
                          kalman::carma::carma_unit_psd<N>(m.shape(), m.kappa(), 0.0);
        const float tau = tune_.tau_applied;
        const float TS = mekf_->get_pseudo_update_period_s();
        if (!(c0 > 1e-14) || !(tau > 0.0f) || !(TS > 0.0f)) return NAN;
        return static_cast<float>(std::sqrt(c0)) * tau / (tuning::R_S_POLE_KAPPA_DEFAULT * std::sqrt(TS));
    }

    void apply_RS_tune_() {
        if (!mekf_) return;
        const float rs = std::min(std::max(tune_.RS_applied, min_R_S_), max_R_S_);
        float rx = rs * R_S_x_factor_, ry = rs * R_S_y_factor_, rz = rs;
        if constexpr (SH) {
            const auto& hm = mekf_->horizontal_model();
            const auto Sg = mekf_->horizontal_sigma();
            const float g2 = hm.noise_gain2();            // unit-variance model gain
            const float vx = rs_pole_law_<NHC>(hm, g2 * Sg(0,0));
            const float vy = rs_pole_law_<NHC>(hm, g2 * Sg(1,1));
            if (std::isfinite(vx) && vx > 0.0f) rx = vx;
            if (std::isfinite(vy) && vy > 0.0f) ry = vy;
        }
        if constexpr (SZ) {
            const auto& vm = mekf_->vertical_model();
            const float v = rs_pole_law_<NZC>(vm, vm.noise_gain2());
            if (std::isfinite(v) && v > 0.0f) rz = v;
        }
        mekf_->set_RS_noise(Eigen::Vector3f(rx, ry, rz));
    }

    void update_tuner(float dt, float a_vertical_measurement, float freq_hz_for_tuner) {
        stepVarianceChannel_(dt, a_vertical_measurement, freq_hz_for_tuner);
        // Legacy keeps its cost and its outputs: it reads no evidence.
        if (excitation_mode_ != CarmaExcitationMode::Legacy) update_wave_evidence_(dt, a_vertical_measurement);
        if (!advanceStartupStage_(wavePeriodUsable())) return;
        constexpr float STILL_VAR_DECAY_SEC = 1.0f;
        const float sea_time_sec = measureOperatingPoint_(
            frontEndStill_(), frontEndStillTimeSec_(), STILL_VAR_DECAY_SEC, tau_coeff_, sigma_coeff_);
        const float sigma_t = excitation_sigma_target_();
        const float RS_raw = rs_spectral_mse_target_(tau_target_, sigma_t);
        RS_target_ = enable_clamp_ ? std::min(std::max(RS_raw, min_R_S_), max_R_S_) : RS_raw;
        adapt_mekf(dt, tau_target_, sigma_t, RS_target_, sea_time_sec);
    }

    void adapt_mekf(float dt, float tau_t, float sigma_t, float RS_t, float sea_time_sec) {
        const float alpha = tauSigmaSmoothingAlpha_(dt, sea_time_sec);
        const float RS_sec = seastate::common::adaptiveSmoothingHorizonSec(
            tuning::ADAPT_RS_MULT, tau_t, RS_t, tune_.RS_applied, tuning::ADAPT_RS_SLEW_LOG, dt);
        const float alpha_RS = 1.0f - std::exp(-dt / RS_sec);
        if (excitation_mode_ == CarmaExcitationMode::Legacy) {
            tune_.tau_applied += alpha * (tau_t - tune_.tau_applied);
        } else {
            // The tau EMA scaled by the frequency weight: zero holds the
            // committed frequency exactly, one is the legacy update.
            const float w = frequency_update_weight_();
            tune_.tau_applied += (alpha * w) * (tau_t - tune_.tau_applied);
        }
        tune_.sigma_applied += alpha    * (sigma_t - tune_.sigma_applied);
        tune_.RS_applied    += alpha_RS * (RS_t    - tune_.RS_applied);
        stageOnlineTuneCommit_();
    }

    void resetTrackingState_() {
        resetFrontEnd_();
        resetAdaptation_();
        evidence_.reset();
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
        apply_tune_();
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
    float mag_delay_sec_ = defaults::MAG_DELAY_SEC;
    seastate::common::TiltResetWatchdog tilt_watchdog_{};
    kalman::filter_detail::ImuLeverArm<float> lever_{};

    float rs_qeff_pow_ = std::pow(2.0f * tuning::R_S_ACCEL_NOISE_DENSITY_DEFAULT, 1.0f / 14.0f);
    float min_R_S_ = tuning::MIN_R_S;
    float max_R_S_ = tuning::MAX_R_S;
    float R_S_x_factor_ = 0.72f;
    float R_S_y_factor_ = 0.50f;
    float Sx_factor_ = 1.0f;
    float Sy_factor_ = 1.0f;
    float h_corr_ = 0.0f;

    tuning::TuneState tune_;
    float RS_target_ = NAN;
    float tau_coeff_ = 1.0f;
    float sigma_coeff_ = 0.9f;
    float carma_sigma_gain_ = 1.0f;
    float carma_sigma_floor_ = 0.05f;

    CarmaExcitationMode excitation_mode_ = CarmaExcitationMode::Legacy;
    bool evidence_gate_ = true;
    bool frequency_hold_ = true;
    float carma_numerical_floor_ = 1.0e-5f;
    CarmaSensorNoise sensor_noise_{};
    CarmaEvidenceConfig evidence_cfg_{};
    CarmaWaveEvidence evidence_{};
    float carma_accel_std_applied_ = NAN;

    kalman::filter_detail::GyroSampleMode core_gyro_sample_mode_ = kalman::filter_detail::GyroSampleMode::IntervalRate;
    float core_gyro_sample_period_ = 0.0f;
    void applyGyroSampling_() {
        mekf_->set_gyro_sample_mode(core_gyro_sample_mode_);
        mekf_->set_gyro_sample_period(core_gyro_sample_period_);
    }
    std::unique_ptr<Mekf> mekf_;
};

struct SeaStateFusionConfig_CARMA3 : seastate::common::ProxyStartupFusionConfig {
    float R_S_noise = 1.5f;
    bool schmidt_acc_bias_hold = false;
    float initial_acc_bias_std = NAN; // calibrated residual uncertainty, m/s^2; NaN keeps native prior
};

template<TrackerType trackerT, int NHC, bool SH, int NZC, bool SZ,
         class MekfT = kalman::Kalman3D_Wave_CARMA3<float, NHC, SH, NZC, SZ>>
class SeaStateFusion_CARMA3
    : public seastate::common::ProxyStartupFusion<
          SeaStateFusionFilter_CARMA3<trackerT, NHC, SH, NZC, SZ, MekfT>,
          SeaStateFusionConfig_CARMA3>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    using Config = SeaStateFusionConfig_CARMA3;
    using Impl = SeaStateFusionFilter_CARMA3<trackerT, NHC, SH, NZC, SZ, MekfT>;

    void begin(const Config& cfg) {
        this->beginStartup_(cfg, [](Impl& impl, const Config& c) {
            impl.initialize_ext(c.sigma_a, c.sigma_g, c.sigma_m, c.Pq0, c.Pb0, c.b0, c.R_S_noise,
                                c.gravity_magnitude);
            if (std::isfinite(c.initial_acc_bias_std) && c.initial_acc_bias_std > 0.0f)
                impl.mekf().set_initial_acc_bias_std(c.initial_acc_bias_std);
        });
    }
};

}  // namespace ocean_imu::carma
