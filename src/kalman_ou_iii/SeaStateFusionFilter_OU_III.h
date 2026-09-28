#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy
  Released under the MIT License

  SeaStateFusionFilter_OU_III

  Marine Inertial Navigational System (INS) Filter for IMU

  Combines multiple real-time estimators into a cohesive ocean-state tracker:

    • Quaternion-based attitude and linear motion estimation via
      Kalman3D_Wave_OU_III: the integrated OU chain [v, p, S, a_w] with the
      integral displacement S regularized by the zero pseudo-measurement
      S = 0 + n_S.

    • The shared measurement-only marine front end
      (seastate::common::MarineWaveFrontEnd): private Mahony observer,
      acceleration-band frequency tracker (wave-direction carrier only),
      WavePeriodEstimator and the two-stage wave-direction stage.  The OU
      operating point is a property of the *wave* band, which the
      acceleration band does not measure: an ocean acceleration spectrum is the
      elevation spectrum weighted by (2πf)⁴, so its apparent frequency sits
      well above the spectral peak and barely moves as the sea grows.  Tuning
      therefore reads WavePeriodEstimator, and reads it at every instant of the
      run -- before that estimator has a value, a fixed wave-band prior stands
      in rather than the tracker.  The tracker has no path into the adaptation.

    • Online auto-tuning of Kalman filter parameters (τ, σ_aw, r_S).  The
      operating point (τ, σ_aw) is measured by the shared adaptation mechanics
      (seastate::common::SeaStateAdaptationCommon): a period-scaled
      measurement-only wave band ahead of SeaStateAutoTuner, whose variance is
      averaged over K wave periods.  r_S follows from the deployed SpectralMSE
      law (RSAdaptationLaw below),
          r_S = C_J q_eff^(1/14) σ_a,B^(6/7) τ^(24/7) / sqrt(T_S),
      with T_S = c_T τ the self-similar pseudo-measurement cadence.  The
      cubic base r_S,base = C_R sqrt(R_a) τ³ is selectable
      as RSAdaptationLaw::Cubic.

  Where
  – τ (tau):  OU process time constant ≈ ½ · T  (half the wave zero-crossing period)
  – σ_aw:     Stationary acceleration scale from period-scaled wave-band RMS
  – r_S:      Pseudo-measurement noise controlling integral drift suppression
  – r_S,xy:   Per-axis horizontal factors (surge/sway anisotropy; see below)

  Adaptive update: τ/σ_aw smoothed over 0.40·T_sea (T_sea = T_z/2; fixed
  seconds retained for ablation), r_S over ADAPT_RS_MULT·τ; the smoothed
  candidate is committed one IMU sample later.

  This header is the OU-III estimator layer: the MEKF, its S=0 regularizer
  law and covariance commits, the magnetometer-gated accelerometer-bias
  release, and the per-sample schedule in updateCore_().  The front end, the
  adaptation mechanics, the vibration conditioning and the proxy-startup
  wrapper (SeaStateFusion_OU_III) are shared with OU-II; see
  src/kalman_common/.

  Features
  • Modular tracker selection via TrackerPolicy template
  • Quaternion-consistent Euler conversion (aerospace → nautical, ENU frame)
  • Magnetometer yaw correction with configurable startup delay
  • Fully compatible with Arduino or native Eigen builds
*/

#ifdef EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#else
#include <ArduinoEigenDense.h>
#endif

#include <cmath>
#include <numbers>
#include <memory>
#include <algorithm>
#include <limits>

#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#include "kalman_common/SeaStateOUFamilyDefaults.h"
#include "kalman_common/AccelVibrationConditioning.h"
#include "kalman_common/MarineWaveFrontEnd.h"
#include "kalman_common/SeaStateAdaptationCommon.h"
#include "kalman_common/ProxyStartupFusion.h"
#include "kalman_common/SeaStateFusionFilterCommon.h"

// Ceiling for the wave-band tuning frequency.  MAX_FREQ_HZ is the tracker's
// bound and does not belong in the adaptation path: with the tuning frequency
// read from the wave band, setFreqBounds() is a wave-direction knob and must
// not move the OU operating point.  Neither bound binds on any reference
// record -- the wave band spans 0.12 to 0.40 Hz -- so this is a safety limit,
// not a tuning surface.  1.2 Hz is a 0.83 s zero-crossing period, shorter than
// any sea a hull responds to while reducing the theorem-only high-frequency tail.
// OU-II and TFG keep 1.5 Hz.
constexpr float MAX_TUNE_FREQ_HZ = 1.2f;

constexpr float MIN_TAU_S   = 0.02f;
// The tau ceiling admits developed seas and long swell when tuning from
// the zero-crossing wave period.
constexpr float MAX_TAU_S   = 12.0f;
constexpr float MAX_SIGMA_A = 4.0f;
// Lower guard against r_S collapsing toward zero in low-motion seas.
constexpr float MIN_R_S     = 0.15f;
// Upper saturation guard for the integral-displacement measurement sigma.
// The Cubic base grows as tau^3; SpectralMSE grows as tau^(24/7) before
// cadence normalization, at fixed acceleration scale.
constexpr float MAX_R_S     = 100.0f;

// Smoothing horizon of the r_S channel, in units of tau_target.  Measured on
// the versioned records against synthesized sea-state transitions: the error
// during a transition falls monotonically as this shortens while the worst
// single realization degrades monotonically the other way, and the value is
// where the mean gain and the worst-record loss cross.  r_S grows at least as
// fast as tau^3 (tau^3 under Cubic, tau^(24/7) under SpectralMSE, before the
// cadence), so it amplifies tau noise by at least the third power, which is
// why this stays above the horizon a tau^1 channel would want.  The common
// tau/sigma_aw EMA (ADAPT_TAU_SEA_PERIODS) and the activation cadence
// (ADAPT_EVERY_SECS) are shared; see SeaStateFusionDefaults.h.
constexpr float ADAPT_RS_MULT              = 1.5f;   // dimensionless
// Discrepancy, in natural-log units of the r_S target-to-applied ratio, above
// which the smoothing horizon shortens.  Zero keeps the plain proportional
// horizon.
constexpr float ADAPT_RS_SLEW_LOG          = 0.0f;   // ln units

// Self-similar S=0 pseudo-measurement cadence; the ratio and the lower clamp
// are shared (SeaStateOUFamilyDefaults.h).  The 150 ms upper guard tightens
// the Live theorem's guaranteed S-update recurrence; OU-II and TFG keep 250 ms.
constexpr float PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT = 0.15f;


// Integral-regularizer adaptation laws.  The first three place the drift-band
// regularization pole of the reduced Riccati model; they differ in which
// asymptotic branch of the posterior acceleration-noise intensity
//     q_eff = 2 r_a (1 - 1/sqrt(1+zeta)),   zeta = 2 sigma_aw^2 tau / r_a
// they assume, where r_a ~ R_a*dt is the acceleration measurement noise
// spectral density seen by the reduced scalar model.
//
//   Cubic            r_S,base = C_R sqrt(R_a) tau^3, then renormalized by
//                    sqrt(T_0/T_S).  With the self-similar cadence
//                    T_S = c_T tau this is effectively r_S ~ tau^(5/2).
//   StrongRiccati    q_eff = 2 r_a: r_S = sqrt(2 r_a) tau^3 / (sqrt(T_S) k^3),
//                    with no leading-order sigma_aw dependence.  At the
//                    analytical C_R this is the *same schedule* as Cubic, not
//                    merely the same shape; see below.
//   PosteriorRiccati the full transition law, reducing to
//                    r_S ~ sigma_aw tau^(7/2) / sqrt(T_S) as zeta -> 0
//                    and to StrongRiccati as
//                    zeta -> infinity.
//
// The deployed envelope has zeta ~ 1e5..1e7 against the bench sensor floor, so
// the strong-observation branch q_eff ~ 2 r_a is the applicable one and the
// drift-driving error is set by the *sensor*, not by the sea.  The base
// schedule takes its acceleration scale accordingly: R_a is the accelerometer
// measurement-noise variance, and sqrt(R_a) tau^3 carries the required units of
// m*s.  The wave amplitude still enters the filter through the OU prior,
// sigma_aw = c_sigma sigma_a,B.  SpectralMSE also uses the physical wave
// amplitude in its displacement-distortion cost.
//
// Writing the base this way makes C_R a pole placement rather than a gain.
// With r_a = R_a h, the cadence-normalized base is
//     C_R sqrt(R_a) tau^3 sqrt(T_0/T_S) = sqrt(2 r_a) tau^3 / (kappa^3 sqrt(T_S))
// exactly when C_R = sqrt(2h/T_0)/kappa^3, so C_R and the normalized corner
// kappa = omega_R tau are two spellings of one number and the Cubic and
// StrongRiccati laws coincide there.  R_S_COEFF_ANALYTICAL_REFERENCE is that
// value for kappa = 0.3627.  The applied C_R comes from a complete-MEKF sweep
// against it, because the scalar reduction omits attitude/gravity leakage,
// residual bias, three-axis covariance coupling and the cadence clamps.
//
// The amplitude tilt of the Riccati laws is a one-parameter ablation:
// p = 0 leaves the Riccati schedule unchanged; p = 1 adds a normalized
// sigma_aw factor.  This tilt does not apply to SpectralMSE.
//   SpectralMSE      the bias-variance law.  The three above all answer "what
//                    r_S holds a chosen normalized pole?"; none answers "what
//                    pole minimizes displacement error for a sea of amplitude
//                    sigma_a?".  The pseudo-measurement both suppresses
//                    integration drift and distorts the real wave, and the
//                    second cost scales with physical wave energy, so
//                        J(omega_R) = 3 q_eff / (2 omega_R^3) + J_wave,
//                        J_wave = (1/2pi) int |G(jw)-1|^2 S_eta^(2) dw,
//                    with |G(jw)-1|^2 = (1+4x^2)/(1+x^6), x = w/omega_R.  Well
//                    below the wave band |G-1|^2 -> 4/x^4, so J_wave -> 4 m_-4
//                    omega_R^4 and the balance gives
//                        omega_R*^7 = (9/32) q_eff / m_-4.
//                    For a self-similar sea m_-4 ~ sigma_a^2 tau^8, hence
//                        r_S* ~ C_J q_eff^(1/14) sigma_a^(6/7) tau^(24/7)
//                               / sqrt(T_S),
//                    i.e. tau^(41/14) away from cadence clamps.  This is where
//                    sigma_a belongs analytically: not in the noise floor, but
//                    in the penalty for suppressing genuine displacement.
//                    DEPLOYED DEFAULT.
//
// Cubic and SpectralMSE are the two laws intended for deployment, and the
// choice between them is a cost/accuracy trade rather than a right/wrong one.
// Over the calibrated envelope the two schedules differ only by
//     r_S,MSE / r_S,Cubic ~ (tau^3 / sigma_aw)^(1/7),
// a seventh root that compresses the exponent difference to about 1.3x across
// the whole H_s = 0.27..8.5 m range, and along the calibrated sea-state
// trajectory (sigma_aw ~ tau^1.1) to a ratio varying only as tau^0.27.  The
// measured difference is correspondingly small but real: SpectralMSE improves
// the paired ten-seed vertical endpoint by 0.0263 +/- 0.0137 %Hs.
//
// SpectralMSE costs one powf per tuner update where Cubic costs none:
// 19.8 ns versus 3.2 ns per evaluation on x86-64 (-O3 -march=native), a factor
// of 6.  On a microcontroller without hardware transcendentals the ratio is
// far larger, and the tuner update is one of the few places the filter calls
// one at all.  Cubic is therefore the supported low-cost configuration for
// embedded targets, and it is a documented operating point rather than a
// deprecated one: on this evidence it gives up well under 1 % of the vertical
// endpoint.  Select it with setRSLaw(RSAdaptationLaw::Cubic), which also needs
// R_S_coeff (C_R) rather than C_J.
enum class RSAdaptationLaw : uint8_t {
    Cubic = 0,
    StrongRiccati = 1,
    PosteriorRiccati = 2,
    SpectralMSE = 3,
};

// Normalized regularization pole omega_R*tau targeted by the Riccati laws.
//
// The ratio between the tilted (p = 1) and untilted (p = 0) members of the
// Riccati family is independent of tau,
//     r_S,p=1 / r_S,p=0 = sigma_aw / sigma_ref,
//     sigma_ref = sqrt(2 r_a) / (kappa^3 C_R sqrt(R_a T_0)),
// because the tau^3/sqrt(T_S) factor is common to both.  sigma_ref is the
// amplitude at which the two members cross, so the ablation measures only the
// sigma_aw dependence of the schedule rather than an overall change of
// regularizer gain.  The anchor is set by kappa and c_R; the ablation varies
// the amplitude dependence at that anchor.
constexpr float R_S_POLE_KAPPA_DEFAULT = 0.3627f;
// r_a = R_a * dt for the reduced scalar acceleration observation.  The default
// is the bench accelerometer noise of the validated configuration
// (0.0148 m/s^2)^2 at the nominal 200 Hz sample interval.
constexpr float R_S_ACCEL_NOISE_DENSITY_DEFAULT = 0.0148f * 0.0148f * FREQ_SMOOTHER_DT;
// Analytical pole-placement value of C_R, Eq. (adapt-cR-kappa) of the OU-III
// paper: C_R = sqrt(2h/T_S,0) / kappa^3.  With h = 5 ms (200 Hz), T_S,0 = 15 ms
// and kappa = 0.3627 this is 17.112.  It is the value at which the deployed
// Cubic base, after the sqrt(T_S,0/T_S) cadence renormalization, is *identical*
// to the StrongRiccati law -- the cubic base and the pole-placement law are two
// spellings of the same schedule, and C_R is the spelling that names the corner
// directly.  Quoted here as the reference the calibration sweep is measured
// against, not as the applied value.
constexpr float R_S_COEFF_ANALYTICAL_REFERENCE = 17.112f;

// Coefficient C_J of the SpectralMSE law,
//     r_S = C_J q_eff^(1/14) sigma_a,B^(6/7) tau^(24/7) / sqrt(T_S).
// It absorbs the dimensionless spectral moment int x^-8 Phi_a(x) dx of the
// self-similar acceleration spectrum, raised to 3/7.  Evaluating that moment on
// the eight reference spectra and solving the exact balance record by record
// gives C_J ~ 0.054, so the default is an analytical prediction rather than a
// fitted number; the sweep is what decides whether it is the complete-MEKF
// optimum.
constexpr float R_S_MSE_COEFF_DEFAULT = 0.0538f;

struct TuneState {
    float tau_applied   = 1.1f;    // s
    float sigma_applied = 1e-2f;   // m/s²
    float RS_applied    = 0.5f;    // m*s
};

//  OU-III sea-state fusion filter: shared front end and adaptation mechanics
//  around the OU-III MEKF and its integral regularizer.
template<TrackerType trackerT>
class SeaStateFusionFilter_OU_III
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

    using FrontEnd::startupProxyQuat;
    using FrontEnd::wavePeriodUsable;

    explicit SeaStateFusionFilter_OU_III(bool with_mag = true)
        : Adaptation(Adaptation::EstimatorBounds{
              MAX_TUNE_FREQ_HZ, MIN_TAU_S, MAX_TAU_S, MAX_SIGMA_A,
              PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT}),
          with_mag_(with_mag)
    {
        setAccelVibrationGuard(ACC_VIBRATION_GUARD_HZ_DEFAULT,
                               ACC_VIBRATION_GUARD_POLES_DEFAULT);
        setAccelVibrationRaccGain(ACC_VIBRATION_RACC_GAIN_DEFAULT);
    }

    void initialize(const Eigen::Vector3f& sigma_a,
                    const Eigen::Vector3f& sigma_g,
                    const Eigen::Vector3f& sigma_m)
    {
        mekf_ = std::make_unique<Kalman3D_Wave_OU_III<float>>(sigma_a, sigma_g, sigma_m);
        seastate::common::finalizeInitialization(
            mekf_,
            [this]() { enterCold_(); },
            [this]() { apply_ou_tune_(true); });
    }

    void initialize_ext(const Eigen::Vector3f& sigma_a,
                        const Eigen::Vector3f& sigma_g,
                        const Eigen::Vector3f& sigma_m,
                        float Pq0, float Pb0,
                        float b0, float R_S_noise,
                        float gravity_magnitude)
    {
        gravity_mps2_ = gravity_magnitude;
        mekf_ = std::make_unique<Kalman3D_Wave_OU_III<float>>(sigma_a, sigma_g, sigma_m, Pq0, Pb0, b0, R_S_noise, gravity_magnitude);
        seastate::common::finalizeInitialization(
            mekf_,
            [this]() { enterCold_(); },
            [this]() { apply_ou_tune_(true); });
    }

    void initialize_from_acc(const Eigen::Vector3f& acc_body_ned) {
        if (mekf_) {
            mekf_->initialize_from_acc(acc_body_ned);
        }
    }

    // Time update (IMU integration + frequency tracking).
    //
    // This drives the MEKF, so it belongs after the handoff.  The filter has no
    // degraded warmup configuration any more -- the linear block, the
    // accelerometer covariance and the bias gate are the live ones from the
    // first sample it sees -- so a caller that starts here instead of at
    // updateFrontEnd() is asking a fully configured wave filter to converge
    // from an unknown attitude on an untuned operating point.  Run
    // updateFrontEnd() until isTunerReady(), hand the proxy attitude over with
    // goLive(), then call this.  SeaStateFusion_OU_III does exactly that.
    void updateTime(float dt, const Eigen::Vector3f& gyro, const Eigen::Vector3f& acc,
                    float tempC = 35.0f)
    {
        updateCore_(dt, gyro, acc, tempC, /*drive_mekf=*/true);
    }

    // Measurement-only front end: the Mahony proxy, the frequency tracker and
    // its stillness detector, the wave-period estimator, the sigma band, the
    // auto-tuner and the wave-direction stage all advance, and the MEKF is
    // left untouched.
    //
    // Every one of those consumers is already exogenous by design -- none of
    // them reads a filter state -- so running them without the MEKF changes
    // none of their outputs.  The one apparent exception is the wave-direction
    // stage, which needs a levelled heading frame: it is given the proxy
    // attitude instead, and heading_frame_acceleration() resolves into the
    // projected bow axis, which is invariant under q -> Rz(psi) q.  The proxy's
    // drifting yaw therefore cannot reach it, and the handoff to the MEKF
    // quaternion later is continuous in everything the stage can see.
    //
    // Model parameters staged by the tuner are still written to the MEKF while
    // this runs.  That is the point: they are parameters, not state, so the
    // filter reaches its first real sample already on the right operating
    // point instead of adapting toward it from FREQ_GUESS.
    void updateFrontEnd(float dt, const Eigen::Vector3f& gyro,
                        const Eigen::Vector3f& acc)
    {
        updateCore_(dt, gyro, acc, /*tempC=*/35.0f, /*drive_mekf=*/false);
    }

    // Hand the attitude over and start the MEKF live.
    //
    // q_bw is the bootstrap solution: proxy tilt, carrying the magnetometer's
    // yaw gauge if one has been acquired.  The sigmas describe how well each
    // part is known, and they are genuinely different -- tilt has been
    // integrated through the wave band, yaw is either gauged or arbitrary --
    // which is why the covariance seed is anisotropic rather than the
    // accel-only default.
    //
    // allow_acc_bias only unlocks the accelerometer-bias gate early.  Left
    // false, the filter keeps waiting for its usual count of magnetometer
    // updates after going live, which is the deployed behaviour.
    void goLive(const Eigen::Quaternionf& q_bw,
                float tilt_sigma_rad,
                float yaw_sigma_rad,
                bool allow_acc_bias = false)
    {
        if (!mekf_) return;
        if (!q_bw.coeffs().allFinite()) return;
        if (!(q_bw.norm() > 1e-8f)) return;

        mekf_->initialize_from_attitude(q_bw, tilt_sigma_rad, yaw_sigma_rad);

        if (allow_acc_bias) {
            accel_bias_locked_ = false;
        }

        enterLive_();
    }

private:
    // The per-sample schedule.  Its order is part of the estimator: the
    // pending schedule is committed before this sample's innovation, the
    // accelerometer is conditioned once for every consumer, the MEKF runs
    // before the tuner sees the sample, and the wave-period estimator runs
    // after it so that the tuner reads the previous sample's period.
    void updateCore_(float dt, const Eigen::Vector3f& gyro,
                     const Eigen::Vector3f& acc, float tempC,
                     bool drive_mekf)
    {
        if (!mekf_) return;
        if (!(dt > 0.0f) || !std::isfinite(dt)) return;

        // Predictable online schedule: commit only coefficients staged
        // after the previous IMU sample. The current sample can update
        // the tuner later in this function, but it cannot choose the
        // model/gain schedule used by its own Kalman innovation.
        apply_pending_online_tune_();
        time_ += dt;
        startup_stage_t_ += dt;

        // Strip out-of-band accelerometer vibration before anything reads it.
        // This is the one place raw measurements enter, so filtering here is
        // what keeps the guard's effect describable: the proxy, the MEKF, and
        // the tilt watchdog below all see the same conditioned signal, and no
        // consumer can be left on a different version of the accelerometer.
        //
        // Armed by default, and transparent below its detector's lower rail,
        // in which case acc_in is acc unchanged.
        const Eigen::Vector3f acc_in = conditionAccel_(acc, dt);

        // Private Mahony observer: levelled vertical channel and startup
        // attitude, fed before the MEKF sees the sample.
        levelVertical_(dt, gyro, acc_in);

        // Tell the MEKF how much it should trust that sample before it uses
        // it, so the covariance and the measurement describe the same
        // conditions.  A no-op unless a gain is set and the guard sees
        // machinery.
        if (drive_mekf) apply_racc_vibration_inflation_();

        // MEKF updates first (attitude + latent a_w)
        if (drive_mekf) {
            mekf_->time_update(gyro, dt);
            mekf_->measurement_update_acc_only(acc_in, tempC);
        }

        // The tilt watchdog reads and rewrites MEKF attitude, so it only has
        // meaning while the MEKF is the one propagating it.  A bootstrap that
        // has not handed over yet has no attitude here to run away.
        if (drive_mekf && tilt_watchdog_.step(mekf_->quaternion_boat(), dt)) {
            if (startup_stage_ == StartupStage::Live) {
                // In Live, re-lock only tilt while preserving yaw/north frame.
                mekf_->initialize_from_acc_preserve_yaw(acc_in);
            } else {
                // During startup stages, accel-only re-lock is acceptable.
                mekf_->initialize_from_acc(acc_in);
                enterCold_();
                resetTrackingState_();
            }
        }

        const float a_vert_measurement = trackAccelBandFrequency_(dt);

        // The sigma channel sees the same exogenous levelled acceleration
        // used by the wave-period estimator, but only after a period-scaled
        // band-pass.  The band is inside update_tuner because its corners use
        // the tuner's own lagged/smoothed wave frequency.
        if (enable_tuner_) {
            update_tuner(dt, a_vert_measurement, tuner_frequency_hz_());
        }

        // R_S is committed with the rest of the staged online schedule at
        // the beginning of the next IMU sample.

        // Bounded covariance inflation of the a_w marginal.
        if (drive_mekf) {
            periodic_aw_cov_sync_tick_();
        }

        // Direction may use the MEKF frame; the default tuner channels do not.
        // Before handoff the MEKF has no attitude to offer, so the proxy's is
        // used.
        if (drive_mekf) {
            const Eigen::Vector3f acc_bias = mekf_->get_acc_bias_body_at_temperature(tempC);
            updateWavePeriodAndDirection_(dt, mekf_->quaternion_boat(), acc_in, &acc_bias);
        } else {
            updateWavePeriodAndDirection_(dt, startupProxyQuat(), acc_in, nullptr);
        }
    }

public:
    //  Magnetometer correction
    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!with_mag_ || !mekf_) return;
        if (time_ < mag_delay_sec_) return;

        mekf_->measurement_update_mag_only(mag_body_ned);
        // Keep counting attempts up to the largest configurable threshold.
        // Measurements and release checks continue after saturation.
        if (mag_updates_applied_ < std::numeric_limits<int>::max()) {
            ++mag_updates_applied_;
        }

        if (!std::isfinite(first_mag_update_time_)) {
            first_mag_update_time_ = static_cast<float>(time_);
        }

        // We can "unlock" once mag has had a few updates, but we DO NOT
        // enable accel-bias learning unless we're already Live.
        if (accel_bias_locked_ &&
            startup_stage_ == StartupStage::Live &&
            mag_updates_applied_ >= mag_updates_to_unlock_ &&
            std::isfinite(first_mag_update_time_) &&
            (static_cast<float>(time_) - first_mag_update_time_) > 1.0f) // 1s guard
        {
            accel_bias_locked_ = false;

            // Only allow accel bias to start learning once no external hold is
            // in force.
            if (!acc_bias_hold_) {
                mekf_->set_acc_bias_updates_enabled(true);
            }
        }
    }

    void setWithMag(bool with_mag) {
        with_mag_ = with_mag;
    }

    // Anisotropy configuration (runtime)
    // S-factor scales horizontal vs vertical stationary std of a_w.
    // The R_S x/y factors scale pseudo-measurement noise per horizontal axis
    // against Z.
    void setSFactor(float s) {
        if (std::isfinite(s) && s > 0.0f) {
            S_factor_ = s;
        }
    }

    // Horizontal integral-noise factors are clamped to [0, 4].  Values above
    // one permit a looser horizontal anchor.  With sigma_aw,H = S_factor *
    // sigma_aw,Z, the similarity law sigma_S ~ sigma_aw tau^3 gives equal
    // dimensionless regularization on all axes when both factors equal S_factor.
    void setRSXFactor(float k) {
        if (std::isfinite(k)) {
            R_S_x_factor_ = std::min(std::max(k, 0.0f), 4.0f);
        }
    }
    void setRSYFactor(float k) {
        if (std::isfinite(k)) {
            R_S_y_factor_ = std::min(std::max(k, 0.0f), 4.0f);
        }
    }

    // tau = tau_coeff * T_z/2 and sigma_aw = sigma_coeff * sigma_a,B.
    void setTauCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) tau_coeff_ = c;
    }
    void setSigmaCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) sigma_coeff_ = c;
    }

    // C_R of the Cubic law.  The in-place rescale keeps the staged and applied
    // values on the same schedule when C_R moves.
    void setRSCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) {
            const float prev = R_S_coeff_;
            R_S_coeff_ = c;

            if (std::isfinite(prev) && prev > 0.0f) {
                const float scale = c / prev;

                if (std::isfinite(tune_.RS_applied) && tune_.RS_applied > 0.0f) {
                    tune_.RS_applied *= scale;
                }
                if (std::isfinite(RS_target_) && RS_target_ > 0.0f) {
                    RS_target_ *= scale;
                }

                apply_RS_tune_();
            }
        }
    }

    // Select how the a_w marginal is re-aligned when a sync happens.
    //
    // false (default, deployed): overwrite the marginal and keep the raw
    // cross-covariances, which rescales the implied correlations by the square
    // root of the marginal change.
    // true: congruence re-alignment, which reaches the same marginal while
    // leaving the whitened cross-covariance untouched and staying PSD by
    // construction.  This is the consistent operation, so running it isolates
    // "is the re-alignment inconsistent?" from "is re-aligning at all a good
    // idea?".
    void setAwCovarianceSyncCongruent(bool flag) { congruent_aw_cov_sync_ = flag; }
    bool awCovarianceSyncCongruent() const noexcept { return congruent_aw_cov_sync_; }

    // Self-similar integral pseudo-measurement cadence T_S = c_T * tau_applied.
    // Enabled by default; disabling selects a fixed 15 ms cadence.
    // Whenever cadence changes while Live,
    // reapply R_S so its per-update covariance stays information-rate matched.
    void setTauScaledPseudoUpdateCadence(bool flag) {
        tau_scaled_pseudo_cadence_ = flag;
        if (!mekf_) return;
        if (flag) apply_pseudo_update_cadence_(mekf_.get(), tune_.tau_applied);
        else mekf_->set_pseudo_update_period_s(pseudo_update_fixed_period_s_);
        if (startup_stage_ == StartupStage::Live) {
            apply_RS_tune_();
        }
    }
    void setPseudoUpdateTauRatio(float ratio) {
        if (!(std::isfinite(ratio) && ratio > 0.0f)) return;
        pseudo_update_tau_ratio_ = ratio;
        if (tau_scaled_pseudo_cadence_) {
            apply_pseudo_update_cadence_(mekf_.get(), tune_.tau_applied);
            if (startup_stage_ == StartupStage::Live) {
                apply_RS_tune_();
            }
        }
    }
    void setPseudoUpdatePeriodBounds(float min_s, float max_s) {
        if (!(std::isfinite(min_s) && std::isfinite(max_s) &&
              min_s > 0.0f && max_s >= min_s)) return;
        pseudo_update_period_min_s_ = min_s;
        pseudo_update_period_max_s_ = max_s;
        if (tau_scaled_pseudo_cadence_) {
            apply_pseudo_update_cadence_(mekf_.get(), tune_.tau_applied);
            if (startup_stage_ == StartupStage::Live) {
                apply_RS_tune_();
            }
        }
    }

    // Select which drift-band regularization law generates the r_S target.
    // Changing it while Live restages the schedule on the next adapt tick; the
    // active covariance is refreshed here so ablations switch cleanly.
    void setRSLaw(RSAdaptationLaw law) {
        rs_law_ = law;
        if (mekf_ && startup_stage_ == StartupStage::Live) {
            apply_RS_tune_();
        }
    }
    RSAdaptationLaw getRSLaw() const noexcept { return rs_law_; }
    void setRSPoleKappa(float kappa) {
        if (!(std::isfinite(kappa) && kappa > 0.0f)) return;
        rs_pole_kappa_ = kappa;
    }
    float getRSPoleKappa() const noexcept { return rs_pole_kappa_; }
    // r_a = R_a * dt of the reduced scalar acceleration observation.
    void setRSAccelNoiseDensity(float r_a) {
        if (!(std::isfinite(r_a) && r_a > 0.0f)) return;
        rs_accel_noise_density_ = r_a;
        refresh_qeff_pow_();
    }
    float getRSAccelNoiseDensity() const noexcept { return rs_accel_noise_density_; }
    // C_J of the SpectralMSE law.
    void setRSMseCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) rs_mse_coeff_ = c;
    }
    float getRSMseCoeff() const noexcept { return rs_mse_coeff_; }
    // Amplitude tilt exponent p of the Riccati laws: r_S ~ sigma^p tau^(5/2).
    // p = 0 is the strong-observation asymptote; p = 1 reproduces Cubic.
    void setRSSigmaExponent(float p) {
        if (std::isfinite(p)) rs_sigma_exponent_ = p;
    }
    float getRSSigmaExponent() const noexcept { return rs_sigma_exponent_; }
    // Diagnostic: the acceleration-information ratio zeta = 2 sigma^2 tau / r_a
    // at the current applied operating point.  zeta >> 1 is the strong
    // acceleration-observation branch.
    float getAccelInformationRatio() const noexcept {
        const float r_a = rs_accel_noise_density_;
        if (!(std::isfinite(r_a) && r_a > 0.0f)) return NAN;
        return 2.0f * tune_.sigma_applied * tune_.sigma_applied
               * tune_.tau_applied / r_a;
    }

    float getPseudoUpdatePeriodSec() const noexcept {
        return mekf_ ? mekf_->get_pseudo_update_period_s() : NAN;
    }

    // Freeze the online tuner at an externally supplied operating point. This
    // is primarily useful for controlled ablations (fixed-nominal and
    // fixed-oracle) after the normal startup sequence has reached Live.
    bool setFixedTuning(float tau_s, float sigma_a, float RS)
    {
        if (!(std::isfinite(tau_s) && tau_s > 0.0f &&
              std::isfinite(sigma_a) && sigma_a > 0.0f &&
              std::isfinite(RS) && RS > 0.0f))
        {
            return false;
        }

        enable_tuner_ = false;
        online_tune_apply_pending_ = false;
        tau_target_ = enable_clamp_
            ? std::min(std::max(tau_s, min_tau_s_), max_tau_s_)
            : tau_s;
        sigma_target_ = enable_clamp_
            ? std::min(sigma_a, max_sigma_a_)
            : sigma_a;
        RS_target_ = enable_clamp_
            ? std::min(std::max(RS, min_R_S_), max_R_S_)
            : RS;

        tune_.tau_applied = tau_target_;
        tune_.sigma_applied = sigma_target_;
        tune_.RS_applied = RS_target_;
        freeze_ou_channel_ = false;
        freeze_RS_channel_ = false;
        apply_ou_tune_(true);
        if (startup_stage_ == StartupStage::Live) {
            apply_RS_tune_();
        }
        return true;
    }

    // Freeze one adaptation channel while the other keeps tracking the sea.
    //
    // The deployed law derives r_S from the live tau and sigma_aw (SpectralMSE;
    // the Cubic base depends on tau alone), so simply freezing tau and sigma_a
    // with setFixedTuning() freezes r_S as well, and an ablation built that way
    // cannot say which channel carries the benefit.  Here the tuner keeps
    // running and keeps deriving r_S from its *live* tau and sigma_a
    // estimates; only the channels named below are held at the supplied
    // operating point on the way to the filter.
    //
    // freeze_ou   holds the OU process parameters (tau, sigma_a) at
    //             (tau_s, sigma_a) while r_S continues to adapt.
    // freeze_RS   holds the integral pseudo-measurement scale at RS while
    //             tau and sigma_a continue to adapt.
    //
    // Freezing both is equivalent to setFixedTuning() and is rejected here so
    // that the two entry points do not silently overlap.
    bool setChannelFreeze(bool freeze_ou,
                          float tau_s,
                          float sigma_a,
                          bool freeze_RS,
                          float RS)
    {
        if (freeze_ou == freeze_RS) return false;
        if (freeze_ou && !(std::isfinite(tau_s) && tau_s > 0.0f &&
                           std::isfinite(sigma_a) && sigma_a > 0.0f))
        {
            return false;
        }
        if (freeze_RS && !(std::isfinite(RS) && RS > 0.0f)) return false;

        enable_tuner_ = true;
        freeze_ou_channel_ = freeze_ou;
        freeze_RS_channel_ = freeze_RS;

        if (freeze_ou) {
            frozen_tau_s_ = enable_clamp_
                ? std::min(std::max(tau_s, min_tau_s_), max_tau_s_)
                : tau_s;
            frozen_sigma_a_ = enable_clamp_
                ? std::min(sigma_a, max_sigma_a_)
                : sigma_a;
            tau_target_ = frozen_tau_s_;
            sigma_target_ = frozen_sigma_a_;
            tune_.tau_applied = tau_target_;
            tune_.sigma_applied = sigma_target_;
            apply_ou_tune_(true);
        }
        if (freeze_RS) {
            frozen_RS_ = enable_clamp_
                ? std::min(std::max(RS, min_R_S_), max_R_S_)
                : RS;
            RS_target_ = frozen_RS_;
            tune_.RS_applied = RS_target_;
            if (startup_stage_ == StartupStage::Live) {
                apply_RS_tune_();
            }
        }
        return true;
    }

    bool frozenOUChannel() const noexcept { return freeze_ou_channel_; }
    bool frozenRSChannel() const noexcept { return freeze_RS_channel_; }

    void setRSBounds(float min_RS, float max_RS) {
        if (!std::isfinite(min_RS) || !std::isfinite(max_RS)) return;
        if (min_RS <= 0.0f || max_RS <= min_RS) return;
        min_R_S_ = min_RS;
        max_R_S_ = max_RS;
    }

    // Smoothing-horizon multiplier for the r_S channel.  The EMA time
    // constant is mult * tau_target, so the horizon follows the sea state
    // instead of being pinned to one second count.  r_S grows at least as
    // fast as tau^3 and so amplifies a tau error by at least the third power,
    // which is what sets the multiplier apart from the tau/sigma one.
    void setRSAdaptMult(float m) {
        if (std::isfinite(m) && m > 0.0f) adapt_RS_mult_ = m;
    }

    // Size, in natural-log units of the r_S target-to-applied ratio, of a
    // discrepancy the smoother should treat as a real sea-state move rather
    // than tuner jitter.  Zero or negative leaves the plain proportional
    // horizon.  See seastate::common::adaptiveSmoothingHorizonSec.
    void setRSAdaptSlewLog(float d) {
        if (std::isfinite(d)) adapt_RS_slew_log_ = d;
    }

    float getRSAdaptMult() const noexcept { return adapt_RS_mult_; }
    float getRSAdaptSlewLog() const noexcept { return adapt_RS_slew_log_; }

    void setMagDelaySec(float delay_sec) {
        if (std::isfinite(delay_sec) && delay_sec >= 0.0f) {
            mag_delay_sec_ = delay_sec;
        }
    }

    // Magnetometer updates that must land after going live before the
    // accelerometer-bias gate opens.  The bias is only weakly observable in
    // waves, so this is a real tuning knob rather than a formality; exposed so
    // it can be swept without editing the filter.
    void setMagUpdatesToUnlockAccBias(int n) {
        if (n >= 0) mag_updates_to_unlock_ = n;
    }
    int magUpdatesToUnlockAccBias() const noexcept { return mag_updates_to_unlock_; }

    // External hold on accelerometer-bias learning, over and above the
    // magnetometer-update count.
    //
    // Accelerometer bias and a tilt error are barely separable in waves, and
    // the bias state has a 5000 s correlation time, so a wrong value learned
    // against a provisional magnetic reference is not a transient -- it
    // outlives the record.  A caller that intends to replace that reference
    // later therefore has to keep this shut until it has, or the bias absorbs
    // an error the reference correction can no longer undo.
    void setAccBiasHold(bool hold) {
        if (acc_bias_hold_ == hold) return;
        acc_bias_hold_ = hold;

        if (!mekf_) return;

        if (hold) {
            mekf_->set_acc_bias_updates_enabled(false);
            return;
        }

        // Releasing does not itself grant learning; the normal gate still
        // decides, and updateMag() re-evaluates it on the next sample.
        if (!accel_bias_locked_ && startup_stage_ == StartupStage::Live) {
            mekf_->set_acc_bias_updates_enabled(true);
        }
    }
    bool accBiasHeld() const noexcept { return acc_bias_hold_; }

    // For SeaStateFusion_OU_III to restore Racc automatically
    void setNominalRaccStd(const Eigen::Vector3f& r) { Racc_nominal_ = r; }

    //  Exposed getters
    inline float getTauApplied()    const noexcept { return mekf_ ? mekf_->get_aw_time_constant() : NAN; }
    inline float getSigmaApplied()  const noexcept { return mekf_ ? mekf_->get_aw_stationary_std().z() : NAN; }
    inline float getRSApplied()     const noexcept { return mekf_ ? mekf_->get_RS_noise_std().z() : NAN; }
    inline float getRSTarget()      const noexcept { return RS_target_;    }

    inline float getHeaveAbs() const noexcept { if (!mekf_) return NAN; return std::fabs(mekf_->get_position().z()); }

    inline float getDisplacementScale(bool smoothed = true) const noexcept {
        const float tau = smoothed ? tune_.tau_applied : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;
        if (!std::isfinite(sigma) || !std::isfinite(tau)) return NAN;
        constexpr float C_HS  = 2.0f * std::numbers::sqrt2_v<float> / (std::numbers::pi_v<float> * std::numbers::pi_v<float>);
        return C_HS * sigma * tau * tau / 2.0f;
    }

    float getVerticalSpeedEnvelopeMps(bool smoothed = true) const noexcept {
        const float tau   = smoothed ? tune_.tau_applied   : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;
        if (!(tau > 1e-6f) || !std::isfinite(tau) || !std::isfinite(sigma)) return NAN;
        constexpr float K = std::numbers::sqrt2_v<float> / std::numbers::pi_v<float>;
        const float v_env = K * sigma * tau;
        return std::isfinite(v_env) ? v_env : NAN;
    }

    inline auto& mekf() noexcept { return *mekf_; }
    inline const auto& mekf() const noexcept { return *mekf_; }

private:
    // Apply the online tuner output only at the next IMU-sample boundary.
    // adapt_mekf() may consume y_k and smooth its candidate during step k,
    // but this function runs before y_{k+1} reaches the MEKF. Therefore the
    // active schedule at step k+1 is measurable with respect to data through k.
    void apply_pending_online_tune_() {
        if (!online_tune_apply_pending_ || !mekf_) return;
        apply_ou_tune_(false);
        if (startup_stage_ == StartupStage::Live) {
            apply_RS_tune_();
        }
        online_tune_apply_pending_ = false;
    }

    // sync_covariance is set only by discrete reconfiguration events. The
    // periodic adaptation path leaves the posterior a_w marginal alone; the
    // new stationary scale reaches the filter through the OU process
    // covariance instead.
    void apply_ou_tune_(bool sync_covariance) {
        if (!mekf_) return;
        mekf_->set_aw_time_constant(tune_.tau_applied);
        // Commit the S=0 cadence with the same applied tau so T_S/tau remains
        // constant apart from explicit safety clamps.
        apply_pseudo_update_cadence_(mekf_.get(), tune_.tau_applied);

        const float sigma_floor = std::max(0.05f, band_noise_floor_sigma_());
        const float sZ = std::max(sigma_floor, tune_.sigma_applied);
        const float sH = sZ * S_factor_;
        const Eigen::Vector3f aw_std(sH, sH, sZ);
        mekf_->set_aw_stationary_std(aw_std);
        if (sync_covariance) {
            apply_aw_cov_sync_();
            last_aw_cov_sync_sec_ = time_;
        }
    }

    // Re-align the posterior a_w marginal with the stationary prior at the
    // adaptation cadence. Runs independently of the tuner so that
    // fixed-tuning modes apply the same policy and remain matched controls.
    void periodic_aw_cov_sync_tick_() {
        if (!mekf_ || !periodicAwCovSyncDue_()) return;
        apply_aw_cov_sync_();
        last_aw_cov_sync_sec_ = time_;
    }

    void apply_aw_cov_sync_() {
        if (congruent_aw_cov_sync_) {
            mekf_->synchronize_aw_covariance_to_stationary_congruent();
        } else {
            mekf_->synchronize_aw_covariance_to_stationary();
        }
    }

    // set_RS_noise() accepts a standard deviation, so one S=0 update has
    // covariance r_S^2. With updates every T_S seconds, the continuous-equivalent
    // information rate is proportional to 1/(r_S^2 T_S).  The Cubic base
    // preserves the nominal 15 ms information rate by normalizing the
    // filter-input standard deviation:
    //     r_S,filter = r_S,base * sqrt(T_0/T_S).
    // The base tuner value remains clamped to [min_R_S_, max_R_S_].  Do not
    // clamp again after this normalization: the smallest-sea operating point
    // may already sit on the base floor and must move below it when
    // T_S > 15 ms to preserve r_S^2 T_S.  With unclamped T_S proportional to
    // tau this turns the Cubic base C_R*sqrt(R_a)*tau^3 schedule into an
    // effective tau^(5/2) schedule at the filter input.
    float pseudo_update_information_rate_scale_() const noexcept {
        // SpectralMSE and the Riccati laws already contain the realized T_S,
        // so renormalizing them again would double-count the cadence.
        if (rs_law_ != RSAdaptationLaw::Cubic) return 1.0f;
        if (!tau_scaled_pseudo_cadence_ || !mekf_) return 1.0f;
        return cadenceRenormalization_(mekf_->get_pseudo_update_period_s());
    }

    // sqrt(R_a), the accelerometer measurement-noise standard deviation the
    // Cubic base uses as its acceleration scale.  The configured quantity is
    // the reduced scalar density r_a = R_a * h, so this divides the sample
    // period back out.  Falls back to the compiled default rather than to
    // zero, because a zero here would silently disable the drift-band
    // regularizer.
    float rs_accel_noise_sigma_() const noexcept {
        const float h = FREQ_SMOOTHER_DT;
        float r_a = rs_accel_noise_density_;
        if (!(std::isfinite(r_a) && r_a > 0.0f)) r_a = R_S_ACCEL_NOISE_DENSITY_DEFAULT;
        return std::sqrt(r_a / h);
    }

    // Bias-variance optimal r_S, Eq. (spectral-mse-rs):
    //     r_S = C_J q_eff^(1/14) sigma_a,B^(6/7) tau^(24/7) / sqrt(T_S).
    // Returned as the *filter input*, like the Riccati laws, because the
    // cadence is already in it; pseudo_update_information_rate_scale_() must
    // therefore not renormalize it again.
    //
    // sigma_a,B is the physical band-limited acceleration RMS, which is what
    // carries the wave energy the pseudo-measurement distorts.  The tuner
    // passes the OU prior sigma_aw = c_sigma sigma_a,B, so divide c_sigma back
    // out: the distortion penalty is a property of the sea, not of our choice
    // of prior, and C_J should not move when c_sigma does.
    float rs_spectral_mse_target_(float tau, float sigma) const noexcept {
        const float TS = pseudo_update_period_for_(tau);
        const float c_sigma = (std::isfinite(sigma_coeff_) && sigma_coeff_ > 0.0f)
                              ? sigma_coeff_ : 1.0f;
        const float sigma_aB = std::max(sigma / c_sigma, 1e-6f);
        if (!(TS > 0.0f) || !(tau > 0.0f)) return rs_mse_coeff_;
        // sigma^(6/7) tau^(24/7) == (sigma tau^4)^(6/7) exactly, so the whole
        // schedule needs one transcendental rather than three: q_eff^(1/14) is
        // constant for a given sensor and is cached in rs_qeff_pow_.  This
        // matters on a microcontroller without hardware transcendentals, where
        // powf dominates the tuner update; see the RSAdaptationLaw comment on
        // the Cubic law as the low-cost alternative.
        const float tau2 = tau * tau;
        const float u = sigma_aB * tau2 * tau2;
        return rs_mse_coeff_ * rs_qeff_pow_
             * std::pow(u, 6.0f / 7.0f)
             / std::sqrt(TS);
    }

    // q_eff^(1/14) for the current sensor noise density, cached because it is
    // constant between setRSAccelNoiseDensity() calls.
    void refresh_qeff_pow_() noexcept {
        float r_a = rs_accel_noise_density_;
        if (!(std::isfinite(r_a) && r_a > 0.0f)) r_a = R_S_ACCEL_NOISE_DENSITY_DEFAULT;
        rs_qeff_pow_ = std::pow(2.0f * r_a, 1.0f / 14.0f);
    }

    // r_S target in m*s for the selected adaptation law.  For the Riccati laws
    // this is Eq. (pole-target) of the OU-III paper,
    //     r_S = sqrt(q_eff) * tau^3 / (sqrt(T_S) * kappa^3),
    // evaluated at the cadence the scheduler will actually use.
    //
    // The Cubic base is Eq. (adapt-rs-base) of the OU-III paper,
    //     r_S,base = C_R sqrt(R_a) tau^3,
    // where R_a is the per-sample accelerometer measurement-noise variance.
    // The acceleration scale is the *sensor* noise, not the measured wave
    // amplitude; see the RSAdaptationLaw comment.
    float rs_target_from_law_(float tau, float sigma) const noexcept {
        const float tau3 = tau * tau * tau;
        if (rs_law_ == RSAdaptationLaw::Cubic) {
            return R_S_coeff_ * rs_accel_noise_sigma_() * tau3;
        }
        if (rs_law_ == RSAdaptationLaw::SpectralMSE) {
            return rs_spectral_mse_target_(tau, sigma);
        }
        const float r_a = rs_accel_noise_density_;
        const float kappa = rs_pole_kappa_;
        if (!(std::isfinite(r_a) && r_a > 0.0f) ||
            !(std::isfinite(kappa) && kappa > 0.0f)) {
            // degenerate config: stay on Cubic
            return R_S_coeff_ * rs_accel_noise_sigma_() * tau3;
        }
        float q_eff = 2.0f * r_a;
        if (rs_law_ == RSAdaptationLaw::PosteriorRiccati) {
            const float zeta = 2.0f * sigma * sigma * tau / r_a;
            // 1 - 1/sqrt(1+zeta); use the numerically stable small-zeta form so
            // the weak branch does not vanish into rounding.
            const float w = (zeta < 1e-3f)
                            ? 0.5f * zeta * (1.0f - 0.75f * zeta)
                            : 1.0f - 1.0f / std::sqrt(1.0f + zeta);
            q_eff *= std::max(0.0f, w);
        }
        const float TS = pseudo_update_period_for_(tau);
        if (!(q_eff > 0.0f) || !(TS > 0.0f))
            return R_S_coeff_ * rs_accel_noise_sigma_() * tau3;
        const float k3 = kappa * kappa * kappa;
        float rs = std::sqrt(q_eff) * tau3 / (std::sqrt(TS) * k3);
        // Optional amplitude tilt (sigma/sigma_ref)^p about the anchor
        // amplitude sigma_ref = sqrt(2 r_a) / (kappa^3 C_R sqrt(R_a T_0)), the
        // sigma at which the p = 0 and p = 1 members cross for every tau.
        // p = 0 leaves the selected Riccati law unchanged; p = 1 adds the
        // normalized sigma_aw multiplier.  In the strong-observation case, at the
        // analytical C_R of Eq. (adapt-cR-kappa) the two coincide exactly,
        // sigma_ref = 1, because that C_R is by construction the one that makes
        // the cadence-normalized cubic base equal to the strong-observation
        // pole-placement law.  Intermediate p models a low-frequency
        // acceleration-error density that is itself sea-state dependent.
        const float p = rs_sigma_exponent_;
        if (p != 0.0f && sigma > 0.0f) {
            const float sigma_ref = std::sqrt(2.0f * r_a) /
                (k3 * R_S_coeff_ * rs_accel_noise_sigma_() *
                 std::sqrt(pseudo_update_fixed_period_s_));
            if (std::isfinite(sigma_ref) && sigma_ref > 0.0f) {
                rs *= std::pow(sigma / sigma_ref, p);
            }
        }
        return rs;
    }

    void apply_RS_tune_(float rs_scale = 1.0f) {
        if (!mekf_) return;
        const float s = (std::isfinite(rs_scale) && rs_scale > 0.0f)
                        ? std::min(rs_scale, 1.0f)
                        : 1.0f;
        const float RSbase = std::min(std::max(tune_.RS_applied, min_R_S_), max_R_S_);
        const float RSb = RSbase * pseudo_update_information_rate_scale_();
        const float rs_z = RSb * s;
        mekf_->set_RS_noise(Eigen::Vector3f(
            rs_z * R_S_x_factor_,
            rs_z * R_S_y_factor_,
            rs_z
        ));
    }

    void update_tuner(float dt, float a_vertical_measurement, float freq_hz_for_tuner) {
        stepVarianceChannel_(dt, a_vertical_measurement, freq_hz_for_tuner);

        if (!advanceStartupStage_(wavePeriodUsable())) return;

        // OU-III attenuates the wave variance with a 1 s time constant while
        // the stillness detector reports still water.
        constexpr float STILL_VAR_DECAY_SEC = 1.0f;
        const float sea_time_sec = measureOperatingPoint_(
            frontEndStill_(), frontEndStillTimeSec_(), STILL_VAR_DECAY_SEC,
            tau_coeff_, sigma_coeff_);

        // r_S is derived from the *live* tau and sigma_a estimates, before any
        // channel freeze is applied below.  Deriving it from frozen values
        // instead would make "freeze the OU channel" silently freeze r_S too,
        // which is precisely the confound the channel ablation exists to
        // remove.
        float RS_raw = rs_target_from_law_(tau_target_, sigma_target_);

        if (enable_clamp_) {
            RS_target_ = std::min(std::max(RS_raw, min_R_S_), max_R_S_);
        } else {
            RS_target_ = RS_raw;
        }

        if (freeze_ou_channel_) {
            tau_target_   = frozen_tau_s_;
            sigma_target_ = frozen_sigma_a_;
        }
        if (freeze_RS_channel_) {
            RS_target_ = frozen_RS_;
        }
        adapt_mekf(dt, tau_target_, sigma_target_, RS_target_, sea_time_sec);
    }

    void adapt_mekf(float dt, float tau_t, float sigma_t, float RS_t,
                    float sea_time_sec) {
        const float alpha = tauSigmaSmoothingAlpha_(dt, sea_time_sec);

        const float RS_sec = seastate::common::adaptiveSmoothingHorizonSec(
            adapt_RS_mult_, tau_t, RS_t, tune_.RS_applied, adapt_RS_slew_log_, dt);
        const float alpha_RS = 1.0f - std::exp(-dt / RS_sec);

        // Every valid physical sample contributes to every EMA.  The deployed
        // activation cadence is deliberately kept separate from that sampling:
        // it throttles parameter commits, not physical measurements or EMA input.
        tune_.tau_applied   += alpha    * (tau_t   - tune_.tau_applied);
        tune_.sigma_applied += alpha    * (sigma_t - tune_.sigma_applied);
        tune_.RS_applied    += alpha_RS * (RS_t    - tune_.RS_applied);

        // Commit this candidate at the beginning of updateTime(k+1).
        stageOnlineTuneCommit_();
    }

    void resetTrackingState_() {
        resetFrontEnd_();
        resetAdaptation_();
    }

    void enterCold_() {
        startup_stage_   = StartupStage::Cold;
        startup_stage_t_ = 0.0f;

        if (!mekf_) return;

        accel_bias_locked_   = with_mag_;
        mag_updates_applied_ = 0;
        first_mag_update_time_  = NAN;

        mekf_->set_acc_bias_updates_enabled(false);
    }

    // The accelerometer sigma the stage logic wants, before any vibration
    // inflation.  Returns a zero vector when it is not known, which is the
    // signal to leave the commanded covariance alone.
    Eigen::Vector3f racc_base_std_() const {
        if (Racc_nominal_.allFinite() && Racc_nominal_.minCoeff() > 0.0f) {
            return Racc_nominal_;
        }
        return Eigen::Vector3f::Zero();
    }

    // Vibration inflation, scaled by the optional low-wave vessel-response
    // noise weighting.  The tuner is updated later in this sample, so the
    // weighting reads its previous measurement-only schedule, independently
    // of direction RAO on/off.
    void apply_racc_vibration_inflation_() {
        if (!mekf_ || (!(racc_vibration_gain_ > 0.0f) &&
                       !(low_wave_noise_.max_std_scale > 1.0f) && !racc_inflated_)) return;

        const Eigen::Vector3f base = racc_base_std_();
        if (!(base.minCoeff() > 0.0f)) return;

        const Eigen::Vector3f scales = isAdaptiveLive()
            ? low_wave_noise_.scales(tune_.sigma_applied / sigma_coeff_,
                                     tuner_frequency_hz_(), base)
            : Eigen::Vector3f::Ones();
        commandRaccStd_(*mekf_, base, scales);
    }

    void enterLive_() {
        startup_stage_   = StartupStage::Live;
        startup_stage_t_ = 0.0f;

        if (!mekf_) return;
        apply_ou_tune_(true);

        // The linear block has been carried through the bootstrap without ever
        // being propagated -- the MEKF is not driven until now -- so its a_w
        // marginal is still the construction seed and its cross-covariances are
        // stale.  Seat the marginal on the operating point the tuner just
        // committed and drop the cross terms before the first prediction.
        mekf_->reset_aw_covariance_to_stationary();

        const bool allow_bias = !accel_bias_locked_ && !acc_bias_hold_;
        mekf_->set_acc_bias_updates_enabled(allow_bias);

        apply_RS_tune_();
    }

    Eigen::Vector3f Racc_nominal_     = Eigen::Vector3f::Constant(0.0f);

    bool accel_bias_locked_ = true;
    int  mag_updates_applied_ = 0;
    static constexpr int MAG_UPDATES_TO_UNLOCK = 250;
    int  mag_updates_to_unlock_ = MAG_UPDATES_TO_UNLOCK;
    bool acc_bias_hold_ = false;

    bool  with_mag_;
    float mag_delay_sec_ = MAG_DELAY_SEC;
    float first_mag_update_time_ = NAN;

    seastate::common::TiltResetWatchdog tilt_watchdog_{};

    // Per-channel freezes for the partial-adaptation ablation; see
    // setChannelFreeze().  Both false is the deployed fully adaptive filter.
    bool freeze_ou_channel_ = false;
    bool freeze_RS_channel_ = false;
    float frozen_tau_s_ = NAN;
    float frozen_sigma_a_ = NAN;
    float frozen_RS_ = NAN;

    bool congruent_aw_cov_sync_ = false;

    // SpectralMSE is the deployed law: it is the only one of the four that
    // answers which regularization corner minimizes displacement error rather
    // than merely how to hold a corner once chosen, and it is the only
    // configuration measured on this branch that improves the primary endpoint
    // against main (-0.0263 +/- 0.0137 %Hs, n = 90) while passing all eight
    // deterministic quality gates.  Its coefficient is analytical, not fitted.
    RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;
    float rs_pole_kappa_ = R_S_POLE_KAPPA_DEFAULT;
    float rs_accel_noise_density_ = R_S_ACCEL_NOISE_DENSITY_DEFAULT;
    float rs_mse_coeff_           = R_S_MSE_COEFF_DEFAULT;
    float rs_qeff_pow_            =
        std::pow(2.0f * R_S_ACCEL_NOISE_DENSITY_DEFAULT, 1.0f / 14.0f);
    float rs_sigma_exponent_ = 0.0f;

    float min_R_S_                = MIN_R_S;
    float max_R_S_                = MAX_R_S;
    float adapt_RS_mult_          = ADAPT_RS_MULT;
    float adapt_RS_slew_log_      = ADAPT_RS_SLEW_LOG;

    // Per-axis horizontal integral-regularization scale, against the vertical
    // one.  The factors are independent because surge and sway have different
    // vessel responses; the calibrated defaults are 0.72 for X and 0.50 for Y.
    //
    // A fixed +/-30 degree projection puts the same ratio in every record:
    // cos30/sin30 = 1.732 in RMS, independent of sea state.  The records do not
    // do that.  Horizontal displacement RMS x/y runs
    //
    //   H0.27  2.246 / 2.049      H4.0  1.880 / 1.858      (JONSWAP / PM-Stokes)
    //   H1.5   2.332 / 2.214      H8.5  1.613 / 1.684
    //
    // falling monotonically with wavelength and reaching the geometric value
    // only in the longest waves.  Geometry alone cannot vary with period; a
    // hull response can, and this is the shape of one.  The dataset's vessel
    // carries a keel: sway is resisted by it and surge is not, so the hull
    // suppresses lateral motion at short periods and follows the orbit at long
    // ones.  The direction-RAO profile for these records says the same thing
    // from the other side, carrying separate horizontal time constants (1.0 s
    // and 0.7 s) rather than one.
    //
    // Every record also holds yaw at exactly 0 -- mean 0.000, sd 0.000, all
    // eight -- so world x and y ARE the vessel's surge and sway here, and this
    // knob pair is the surge/sway split rather than a world-frame accident.
    // rho_y = 0.50 gives the keel-damped axis the
    // tighter integral anchor.  Pooled over four fresh IMU draws and the eight
    // records: pitch -17 percent, y accelerometer bias -21 percent, 3D
    // accelerometer bias -2 percent, yaw unchanged at 1.003, against roll +4
    // percent and x bias +5 percent.  The vertical channels do not move at all
    // (1.000).  An isotropic 0.5 costs yaw 5 percent on the same draws,
    // which the split does not.
    //
    // The frame matters: the keel's anisotropy
    // is a body-frame property and R_S is applied in world NED, so these two
    // numbers only coincide with surge/sway while the vessel heads north, as it
    // does in every record here.  Rotating the anisotropy into the body frame
    // by the estimated heading is the deployment-general form of this and is
    // not attempted here: with yaw pinned at 0 these records cannot tell the
    // rotation from the constant, so it would ship unmeasured.
    float R_S_x_factor_ = 0.72f;
    float R_S_y_factor_ = 0.50f;
    // Horizontal stationary acceleration scale relative to the vertical one.
    // S_factor = 1 gives an isotropic stationary acceleration prior.  For
    // isotropic r_S this also equalizes r_S/(sigma_aw tau^3) across axes.
    float S_factor_      = 1.0f;

    TuneState tune_;
    float     RS_target_ = NAN;

    // The three adaptation coefficients of Eq. (adapt-three-layer-summary):
    //     tau      = tau_coeff * T_z / 2
    //     sigma_aw = sigma_coeff * sigma_a,B
    //     r_S,base = R_S_coeff * sqrt(R_a) * tau^3        (Cubic law)
    // tau_coeff = 1 is the documented intent, tau equal to half the
    // zero-crossing period.  R_S_coeff is the coefficient of the Cubic law
    // (and the anchor of the Riccati tilt); the deployed SpectralMSE law uses
    // C_J = rs_mse_coeff_ instead.
    //
    // R_S_coeff is C_R, and it is not a bare gain: by Eq. (adapt-cR-kappa) it
    // places the normalized regularizer corner, kappa = omega_R tau =
    // (sqrt(2h/T_S,0)/C_R)^(1/3).  The applied value is the analytical
    // pole-placement value itself, and that is a measured outcome rather than a
    // deference to the theory: C_R = 17.11 is the argmin of every c_sigma row
    // of the coarse calibration sweep, and on the paired multi-seed harness it
    // beats the deterministic sweep argmin C_R = 14.4 on the primary endpoint
    // and on pitch.  The scalar reduction that predicts it omits
    // attitude/gravity leakage, residual bias, three-axis covariance coupling
    // and the cadence clamps, so this agreement was checked, not assumed.
    //
    // sigma_coeff sets the OU prior scale.  SpectralMSE divides that coefficient
    // back out when recovering physical wave RMS for its distortion cost.
    float R_S_coeff_    = R_S_COEFF_ANALYTICAL_REFERENCE;
    float tau_coeff_    = 1.0f;
    float sigma_coeff_  = 0.9f;

    std::unique_ptr<Kalman3D_Wave_OU_III<float>>  mekf_;
};

// Deployed OU-III front end: proxy bootstrap, two-stage magnetic acquisition,
// continuous hard iron and displacement output (shared with OU-II; see
// ProxyStartupFusion.h), around the OU-III filter above.
struct SeaStateFusionConfig_OU_III : seastate::common::ProxyStartupFusionConfig {
    // Initial integral pseudo-measurement variance, (m*s)^2.  The tuner
    // overwrites it at Live; it sets the pre-Live value and seeds the
    // commanded-parameter filter.
    float R_S_noise = 1.5f;
};

template<TrackerType trackerT>
class SeaStateFusion_OU_III
    : public seastate::common::ProxyStartupFusion<SeaStateFusionFilter_OU_III<trackerT>,
                                                  SeaStateFusionConfig_OU_III>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Config = SeaStateFusionConfig_OU_III;

    void begin(const Config& cfg) {
        this->beginStartup_(cfg, [](SeaStateFusionFilter_OU_III<trackerT>& impl, const Config& c) {
            impl.initialize_ext(c.sigma_a, c.sigma_g, c.sigma_m,
                                c.Pq0, c.Pb0, c.b0, c.R_S_noise,
                                c.gravity_magnitude);
        });
    }

    Eigen::Vector3f eulerNauticalDeg() const {
        return this->impl_.getEulerNautical();
    }
};
