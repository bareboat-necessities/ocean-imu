#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy

  Released under the MIT License

  SeaStateFusionFilter_OU_II

  Marine Inertial Navigational System (INS) Filter for IMU

  Combines multiple real-time estimators into a cohesive ocean-state tracker:

    • Quaternion-based attitude and linear motion estimation via
      Kalman3D_Wave_OU_II: the OU chain [v, p, a_w] regularized by the two
      zero pseudo-measurements p = 0 + n_p and v = 0 + n_v.

    • The shared measurement-only marine front end
      (seastate::common::MarineWaveFrontEnd): private Mahony observer,
      acceleration-band frequency tracker with dual-stage smoothing (wave
      direction carrier and reporting only), WavePeriodEstimator and the
      two-stage wave-direction stage.

    • Online auto-tuning of Kalman filter parameters (τ, σ_aw, r_p0, r_v0).
      The operating point (τ, σ_aw) is measured by the shared adaptation
      mechanics (seastate::common::SeaStateAdaptationCommon): a period-scaled
      measurement-only wave band ahead of SeaStateAutoTuner, whose variance is
      averaged over K wave periods, with the time scale taken from
      WavePeriodEstimator at every instant of the run -- a fixed wave-band
      prior before it has a value, never the acceleration-band tracker.  The
      pseudo-measurement pair then follows from the deployed PhysicalMSE law
      of the OU-II dual-regularization note (PseudoAdaptationLaw below),
          r_p = C_P q_eff^(1/10) σ_a,B^(4/5) τ^(12/5) / sqrt(T_S),
          r_v = r_p / (ratio · τ),
      with T_S = c_T τ the self-similar pseudo-measurement cadence.  The
      Empirical law r_p0 ~ σ_aw τ², r_v0 ~ σ_aw τ is also
      selectable.

  Where
  – τ (tau):  OU process time constant = c_τ · T_z/2, scaled from the zero-crossing period
              of the surface elevation as estimated by WavePeriodEstimator.
              It is deliberately *not* half the dominant period of acceleration:
              an ocean acceleration spectrum is the elevation spectrum weighted
              by (2πf)⁴, so its apparent frequency sits well above the spectral
              peak and barely moves as the sea grows, which pins τ to roughly
              one value across every sea state.
  – σ_aw:     Stationary acceleration scale from period-scaled wave-band RMS
  – r_p0:     Pseudo-measurement noise controlling p drift suppression
  – r_v0:     Pseudo-measurement noise controlling v drift suppression
  – r_p0,xy:  Anisotropic X/Y weight on the p pseudo-measurement

  Adaptive update: τ/σ_aw smoothed over 0.40·T_sea (T_sea = T_z/2; fixed
  seconds retained for ablation), r_p0/r_v0 over ADAPT_R_*_MULT·τ; the smoothed
  candidate is committed one IMU sample later.

  ------------------------------------------------------------------------
  WHAT IS SHARED WITH SeaStateFusionFilter_OU_III, AND WHAT IS NOT
  ------------------------------------------------------------------------

  The tuner and the attitude front end sit outside the estimator and do not
  know which of the two filters they are driving, so they are literally the
  same code (src/kalman_common/): the front end, the vibration conditioning,
  the adaptation mechanics -- exogenous one-sample staging, the JONSWAP-similar
  sigma band, the τ-scaled pseudo cadence T_S = (0.015/1.1)·τ clamped to
  [5 ms, T_S,max] -- and the proxy-startup wrapper, in which the proxy owns
  tilt and magnetic learning.  The two filters differ in the translational
  state structure -- OU-III carries the extra integral displacement state and
  regularizes it with a single r_S, this one regularizes p and v separately --
  and in everything that follows from it, which is what this header contains:
  the dual-channel PhysicalMSE law and its covariance commits, the p/v
  smoothing horizons, and the estimator-specific bounds (T_S,max = 250 ms here
  against OU-III's 150 ms, tuning-frequency ceiling 1.5 Hz against
  1.2 Hz, σ_aw ceiling 6 against 4, a configurable 5 s still-water decay
  against OU-III's fixed 1 s).

  Only the Empirical law is renormalized by sqrt(T_0/T_S) to hold the
  nominal 15 ms information rate; PhysicalMSE is derived on the continuous
  densities ρ = r² T_S and already contains the realized cadence.

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

#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
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
// not a tuning surface.  1.5 Hz is a 0.67 s zero-crossing period, shorter than
// the periods in the reference records.  OU-III uses a 1.2 Hz ceiling.
constexpr float MAX_TUNE_FREQ_HZ = 1.5f;

constexpr float MIN_TAU_S     = 0.02f;  // sec
// The tau ceiling admits developed seas and long swell when tuning from
// the zero-crossing wave period.
constexpr float MAX_TAU_S     = 12.0f;  // sec
constexpr float MAX_SIGMA_A   = 6.0f;
// Saturation safeguards for the position and velocity pseudo-measurement
// standard deviations, in m and m/s respectively.  regularizer_floor-test
// checks that the deployed law stays clear of the floors on its stress cases.
constexpr float MIN_R_p0_std  = 0.05f;
constexpr float MAX_R_p0_std  = 150.0f;
constexpr float MIN_R_v0_std  = 0.01f;
constexpr float MAX_R_v0_std  = 40.0f;

// Smoothing horizons of the two drift-correction channels, in units of
// tau_target.  Measured on the versioned records against synthesized sea-state
// transitions: the error during a transition falls monotonically as these
// shorten, the stationary worst-record vertical error rises monotonically, and
// 3.0 is where the two cross at an acceptable cost for both channels.
constexpr float ADAPT_R_p0_MULT            = 3.0f;   // dimensionless
constexpr float ADAPT_R_v0_MULT            = 3.0f;   // dimensionless
// Discrepancy, in natural-log units of the target-to-applied ratio, above
// which the smoothing horizon of either drift-correction channel shortens.
// Zero keeps the plain proportional horizon.
constexpr float ADAPT_R_SLEW_LOG           = 0.0f;   // ln units

// Self-similar drift-regularizer pseudo-measurement cadence.
//
// The S=0 pseudo updates of OU-III and the p=0/v=0 pseudo updates here share
// one mechanism and one 15 ms / 1.1 s nominal operating point, a single
// periodic tick inside the MEKF fires them, and the ratio and lower clamp are
// shared (SeaStateOUFamilyDefaults.h).  One pseudo update has covariance r^2,
// so at one update per T_S seconds the continuous-equivalent information rate
// goes as 1/(r^2 T_S); scaling T_S with tau while holding r fixed would
// therefore change the regularization strength with sea state as a side
// effect, which is why each law states how it refers its targets to T_S.
// The upper guard is inactive over the current tau <= 12 s operating envelope.
constexpr float PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT = 0.25f;

// Dual pseudo-measurement adaptation laws.  Both schedule the standard
// deviations of the two zero pseudo-measurements -- p = 0 and v = 0 -- that
// regularize the low-frequency end of the integration chain.  They differ in
// what they claim those numbers are.
//
//   Empirical    r_p0,base = c_p sigma_aw tau^2,  r_v0,base = c_v sigma_aw tau,
//                both renormalized by sqrt(T_0/T_S).  Dimensionally consistent
//                and swept to a sharp optimum on the reference records, but the
//                powers of sigma_aw and tau were never derived: dimensional
//                analysis alone does not identify them.  Effectively
//                sigma_aw tau^(3/2) and sigma_aw tau^(1/2) at the filter input.
//   PhysicalMSE  the joint displacement-MSE law of the OU-II dual-regularization
//                note, doc/kalman_ou_ii/ou2-dual-regularization-mse.tex.
//                DEPLOYED DEFAULT.
//
// The PhysicalMSE derivation is the OU-II analogue of the OU-III SpectralMSE
// law and rests on the same distinction: the Kalman covariance of a filter that
// believes a zero pseudo-measurement is not that filter's physical error.  The
// pseudo-measurements suppress integration drift and simultaneously distort the
// real wave, and only the second cost scales with physical wave energy.
//
// Treating the two channels jointly rather than as independent penalties, the
// reduced two-state Kalman-Bucy CARE for the double integrator
//     p_dot = v,  v_dot = n_a,   S_na = q_eff
// under continuous pseudo densities rho_p = r_p^2 T_S and rho_v = r_v^2 T_S has
// the closed-form stabilizing solution of Theorem (dual-channel CARE), which is
// parameterized by
//     omega_p = (q/rho_p)^(1/4),   chi = (omega_v/omega_p)^2,
//     omega_v = (q/rho_v)^(1/2).
// chi says how much of one second-order loop is supplied through the velocity
// channel; the two pseudo channels do not generate two independent corners.
// The closed loop reconstructs displacement through
//     G(s) = s^2 H_pa(s) = (1/(1+chi)) s^2 / (s^2 + omega_p sqrt(2+chi) s + omega_p^2),
// so the physical objective is
//     J = J_n + J_w,
//     J_n = q / (2 omega_p^3 (1+chi)^2 sqrt(2+chi)),          exact H_2 norm
//     J_w = (1/2pi) int |G(jw)-1|^2 S_p^(2) dw.
// Expanding both for chi << 1 and a regularization corner below the wave band,
//     J ~ q(1 - 9chi/4)/(2 sqrt2 omega_p^3) + 2 M_-2 omega_p^2 + M_0 chi^2,
// whose stationary point is
//     omega_p*^5 = 3q/(8 sqrt2 M_-2),   chi* = 9q/(16 sqrt2 M_0 omega_p*^3).
// For a self-similar sea M_0 ~ sigma_a^2 tau^4 and M_-2 ~ sigma_a^2 tau^6, so
//     r_p = C_P q_eff^(1/10) sigma_a,B^(4/5) tau^(12/5) / sqrt(T_S),
//     r_v = C_V q_eff^(1/10) sigma_a,B^(4/5) tau^(7/5)  / sqrt(T_S),
// i.e. tau^(19/10) and tau^(9/10) at the filter input away from cadence clamps,
// against the Empirical tau^(3/2) and tau^(1/2).  Note that both channels move
// by the same tau^(2/5): the theory does not ask for a different *relative*
// period law, only a common shape correction and a different amplitude power
// (4/5 rather than 1).
//
// sigma_a,B is the physical band-limited acceleration RMS, which is what
// carries the wave energy the pseudo-measurements distort.  The tuner passes
// the OU prior sigma_aw = c_sigma sigma_a,B, so the law divides c_sigma back
// out: the distortion penalty is a property of the sea and not of our choice of
// prior, and C_P and c_sigma must stay separately identifiable.
//
// Empirical stays selectable with setPseudoLaw(PseudoAdaptationLaw::Empirical),
// which also needs R_p0_coeff/R_v0_coeff (c_p, c_v) rather than C_P and the
// channel ratio.  It costs no transcendental per tuner update where PhysicalMSE
// costs one powf, so it remains the supported low-cost configuration for
// embedded targets without hardware transcendental support.
enum class PseudoAdaptationLaw : uint8_t {
    Empirical   = 0,
    PhysicalMSE = 1,
};

// Coefficient C_P of the PhysicalMSE position channel,
//     r_p = C_P q_eff^(1/10) sigma_a,B^(4/5) tau^(12/5) / sqrt(T_S).
// C_P absorbs the dimensionless spectral moment of the self-similar
// displacement spectrum, C_P = (8 sqrt2 / 3)^(2/5) mu_-2^(2/5) with
// mu_-2 = M_-2 / (sigma_a^2 tau^6). The fixed coefficient is retained;
// tools/ou2_dual_mse_coefficients.py reports the current vessel-CG spectral
// diagnostic separately. It is not a universal hull-independent identity.
constexpr float R_PSEUDO_MSE_COEFF_DEFAULT = 0.1116f;

// Channel ratio C_P/C_V of the PhysicalMSE law.  This is the corollary the
// note identifies as the more tightly constrained half of the prediction,
// because q_eff and the cadence normalization both cancel out of it:
//     (r_p / r_v)^2 = (3/2) M_-2 / M_0,
// hence r_p/r_v = C_P/C_V * tau for a fixed normalized sea shape.  The eight
// vessel-CG spectra motivate a separate ratio diagnostic. The deployed 0.5
// matches the vessel-RAO replay profile and its paired sensor-draw validation.
//
// Applying the ratio rather than a second independent power is exact, not an
// approximation: sigma^(4/5) tau^(7/5) = sigma^(4/5) tau^(12/5) / tau.  It is
// also what keeps the whole schedule down to a single transcendental.
constexpr float R_PSEUDO_MSE_RATIO_DEFAULT = 0.5f;

// q_eff = 2 r_a with r_a = R_a * h, the density of the residual acceleration
// error the integration chain actually sees.
//
// This is not the accelerometer's bench noise spec.  The note defines q as the
// residual error "presented to the integration chain after acceleration
// estimation", which carries attitude and gravity leakage, residual
// accelerometer bias and estimation error on top of the sensor floor -- and
// the reduced model that drops all three from its *dynamics* still needs their
// intensity here.  OU-III's SpectralMSE law uses the bench figure legitimately,
// because its strong-observation branch is a statement about the sensor; the
// OU-II objective instead includes the additional residual-error sources.
//
// The filter already carries a measured estimate of the right quantity:
// ACC_NOISE_FLOOR_SIGMA_DEFAULT, the pre-band vertical acceleration noise
// floor the tuner subtracts as non-wave energy, at 0.12 m/s^2 -- eight times
// the bench 0.0148 m/s^2.  Referring the law to it moves the schedule by
// (0.12/0.0148)^(1/5) = 1.52, and the eight-record scale sweep puts the
// complete-MEKF vertical optimum within 3 % of exactly that.  So the
// noise density and C_P serve distinct roles in the deployed law.
//
// It is a separate knob rather than a live read of acc_noise_floor_sigma_, so
// that the sweep axis stays clean and so that q stays sea-state independent,
// which is the assumption the 4/5 amplitude exponent rests on.  A platform
// re-characterization should move both: they are the same physical number.
//
// The law depends on it only as q_eff^(1/10), so ten times the noise *density*
// moves the schedule by 26 % and ten times the noise standard deviation by
// 58 %.  It cancels out of the channel ratio entirely.
constexpr float R_PSEUDO_ACCEL_NOISE_DENSITY_DEFAULT =
    ACC_NOISE_FLOOR_SIGMA_DEFAULT * ACC_NOISE_FLOOR_SIGMA_DEFAULT * FREQ_SMOOTHER_DT;

struct TuneState {
    float tau_applied      = 1.1f;   // s
    float sigma_applied    = 1e-2f;  // m/s²
    float R_p0_std_applied = 0.1f;   // m
    float R_v0_std_applied = 0.1f;   // m/s
};

// OU-II sea-state fusion filter: shared front end and adaptation mechanics
// around the OU-II MEKF and its dual p/v regularizer.
template<TrackerType trackerT>
class SeaStateFusionFilter_OU_II
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

    explicit SeaStateFusionFilter_OU_II(bool with_mag = true)
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
        mekf_ = std::make_unique<Kalman3D_Wave_OU_II<float>>(sigma_a, sigma_g, sigma_m);
        seastate::common::finalizeInitialization(
            mekf_,
            [this]() { enterCold_(); },
            [this]() { apply_ou_tune_(true); });
    }

    void initialize_ext(const Eigen::Vector3f& sigma_a,
                        const Eigen::Vector3f& sigma_g,
                        const Eigen::Vector3f& sigma_m,
                        float Pq0, float Pb0,
                        float b0, float R_p0_var_init, float R_v0_var_init,
                        float gravity_magnitude)
    {
        gravity_mps2_ = gravity_magnitude;
        mekf_ = std::make_unique<Kalman3D_Wave_OU_II<float>>(
            sigma_a, sigma_g, sigma_m, Pq0, Pb0, b0, R_p0_var_init, R_v0_var_init, gravity_magnitude);
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
    // goLive(), then call this.  SeaStateFusion_OU_II does exactly that.
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

        // R_p0/R_v0 are committed with the rest of the staged online
        // schedule at the beginning of the next IMU sample.

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
    // Magnetometer correction
    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!with_mag_ || !mekf_) return;
        if (time_ < mag_delay_sec_) return;

        mekf_->measurement_update_mag_only(mag_body_ned);
        // Unlike OU-III, this counter is not saturated.
        mag_updates_applied_++;

        if (!std::isfinite(first_mag_update_time_)) {
            first_mag_update_time_ = static_cast<float>(time_);
        }

        // We can "unlock" once mag has had a few updates, but we DO NOT
        // enable accel-bias learning unless we're already Live.
        if (accel_bias_locked_ &&
            startup_stage_ == StartupStage::Live &&
            mag_updates_applied_ >= mag_updates_to_unlock_ &&
            std::isfinite(first_mag_update_time_) &&
            (static_cast<float>(time_) - first_mag_update_time_) > 1.0f)
        {
            accel_bias_locked_ = false;

            // Only allow accel bias to start learning once no external hold is
            // in force.
            if (!acc_bias_hold_) {
                mekf_->set_acc_bias_updates_enabled(true);
            }
        }
    }

    void setWithMag(bool with_mag) { with_mag_ = with_mag; }

    // Anisotropy configuration (runtime)
    // P-factor scales horizontal vs vertical stationary std of a_w.
    // The R_p0 x/y factors scale the position pseudo-measurement noise per
    // horizontal axis against Z.
    void setPFactor(float p) {
        if (std::isfinite(p) && p > 0.0f) P_factor_ = p;
    }

    // The ceiling is 4, matching OU-III's.  A ceiling of 1 would encode the
    // assumption that a horizontal position anchor can only ever be tighter
    // than the vertical one, which is the assumption these knobs exist to
    // measure rather than one they should impose.
    void setR_p0_XFactor(float k) {
        if (std::isfinite(k)) {
            R_p0_x_factor_ = std::min(std::max(k, 0.0f), 4.0f);
        }
    }
    void setR_p0_YFactor(float k) {
        if (std::isfinite(k)) {
            R_p0_y_factor_ = std::min(std::max(k, 0.0f), 4.0f);
        }
    }

    // tau = tau_coeff * T_z/2 and sigma_aw = sigma_coeff * sigma_a,B.
    void setTauCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) tau_coeff_ = c;
    }
    void setSigmaCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) sigma_coeff_ = c;
    }

    // c_p of the Empirical law.  The in-place rescale of the staged and applied
    // values only makes sense while that law is the one generating them; under
    // PhysicalMSE the coefficient is stored for a later switch and nothing is
    // rescaled.
    void setR_p0_Coeff(float c) {
        if (std::isfinite(c) && c > 0.0f) {
            const float prev = R_p0_coeff_;
            R_p0_coeff_ = c;

            if (pseudo_law_ == PseudoAdaptationLaw::Empirical &&
                std::isfinite(prev) && prev > 0.0f) {
                const float scale = c / prev;

                if (std::isfinite(tune_.R_p0_std_applied) && tune_.R_p0_std_applied > 0.0f) {
                    tune_.R_p0_std_applied *= scale;
                }
                if (std::isfinite(R_p0_std_target_) && R_p0_std_target_ > 0.0f) {
                    R_p0_std_target_ *= scale;
                }
                apply_R_p0_tune_();
            }
        }
    }

    // c_v of the Empirical law; see setR_p0_Coeff().
    void setR_v0_Coeff(float c) {
        if (std::isfinite(c) && c > 0.0f) {
            const float prev = R_v0_coeff_;
            R_v0_coeff_ = c;

            if (pseudo_law_ == PseudoAdaptationLaw::Empirical &&
                std::isfinite(prev) && prev > 0.0f) {
                const float scale = c / prev;

                if (std::isfinite(tune_.R_v0_std_applied) && tune_.R_v0_std_applied > 0.0f) {
                    tune_.R_v0_std_applied *= scale;
                }
                if (std::isfinite(R_v0_std_target_) && R_v0_std_target_ > 0.0f) {
                    R_v0_std_target_ *= scale;
                }
                apply_R_v0_tune_();
            }
        }
    }

    // Select which law generates the two pseudo-measurement targets.  Changing
    // it while Live restages the schedule on the next adapt tick; the active
    // covariances are refreshed here so ablations switch cleanly.
    void setPseudoLaw(PseudoAdaptationLaw law) {
        pseudo_law_ = law;
        if (mekf_ && startup_stage_ == StartupStage::Live) {
            apply_R_p0_tune_();
            apply_R_v0_tune_();
        }
    }
    PseudoAdaptationLaw getPseudoLaw() const noexcept { return pseudo_law_; }

    // C_P of the PhysicalMSE position channel.
    void setPseudoMseCoeff(float c) {
        if (std::isfinite(c) && c > 0.0f) pseudo_mse_coeff_ = c;
    }
    float getPseudoMseCoeff() const noexcept { return pseudo_mse_coeff_; }

    // C_P/C_V, the PhysicalMSE channel ratio: r_p/r_v = ratio * tau.
    void setPseudoMseRatio(float r) {
        if (std::isfinite(r) && r > 0.0f) pseudo_mse_ratio_ = r;
    }
    float getPseudoMseRatio() const noexcept { return pseudo_mse_ratio_; }

    // r_a = R_a * dt of the residual acceleration error entering the
    // integration chain; the law uses q_eff = 2 r_a.  This is the same
    // physical number as the acceleration noise floor, not the accelerometer
    // bench spec; see R_PSEUDO_ACCEL_NOISE_DENSITY_DEFAULT.
    void setPseudoAccelNoiseDensity(float r_a) {
        if (!(std::isfinite(r_a) && r_a > 0.0f)) return;
        pseudo_accel_noise_density_ = r_a;
        refresh_pseudo_qeff_pow_();
    }
    float getPseudoAccelNoiseDensity() const noexcept {
        return pseudo_accel_noise_density_;
    }

    // Self-similar drift-regularizer pseudo-measurement cadence
    // T_S = c_T * tau_applied.  Enabled by default; disabling restores the
    // fixed 15 ms cadence for ablation.
    // Whenever the cadence changes while Live, reapply r_p0 and r_v0 so their
    // per-update covariances stay information-rate matched.
    void setTauScaledPseudoUpdateCadence(bool flag) {
        tau_scaled_pseudo_cadence_ = flag;
        if (!mekf_) return;
        if (flag) apply_pseudo_update_cadence_(mekf_.get(), tune_.tau_applied);
        else mekf_->set_pseudo_update_period_s(pseudo_update_fixed_period_s_);
        if (startup_stage_ == StartupStage::Live) {
            apply_R_p0_tune_();
            apply_R_v0_tune_();
        }
    }
    void setPseudoUpdateTauRatio(float ratio) {
        if (!(std::isfinite(ratio) && ratio > 0.0f)) return;
        pseudo_update_tau_ratio_ = ratio;
        if (tau_scaled_pseudo_cadence_) {
            apply_pseudo_update_cadence_(mekf_.get(), tune_.tau_applied);
            if (startup_stage_ == StartupStage::Live) {
                apply_R_p0_tune_();
                apply_R_v0_tune_();
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
                apply_R_p0_tune_();
                apply_R_v0_tune_();
            }
        }
    }
    float getPseudoUpdatePeriodSec() const noexcept {
        return mekf_ ? mekf_->get_pseudo_update_period_s() : NAN;
    }

    // Freeze the online tuner at an externally supplied operating point. This
    // is primarily useful for controlled ablations (fixed-nominal and
    // fixed-oracle) after the normal startup sequence has reached Live.
    bool setFixedTuning(float tau_s,
                        float sigma_a,
                        float R_p0_std,
                        float R_v0_std)
    {
        if (!(std::isfinite(tau_s) && tau_s > 0.0f &&
              std::isfinite(sigma_a) && sigma_a > 0.0f &&
              std::isfinite(R_p0_std) && R_p0_std > 0.0f &&
              std::isfinite(R_v0_std) && R_v0_std > 0.0f))
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
        R_p0_std_target_ = enable_clamp_
            ? std::min(std::max(R_p0_std, MIN_R_p0_std_), MAX_R_p0_std_)
            : R_p0_std;
        R_v0_std_target_ = enable_clamp_
            ? std::min(std::max(R_v0_std, MIN_R_v0_std_), MAX_R_v0_std_)
            : R_v0_std;

        tune_.tau_applied = tau_target_;
        tune_.sigma_applied = sigma_target_;
        tune_.R_p0_std_applied = R_p0_std_target_;
        tune_.R_v0_std_applied = R_v0_std_target_;
        apply_ou_tune_(true);
        if (startup_stage_ == StartupStage::Live) {
            apply_R_p0_tune_();
            apply_R_v0_tune_();
        }
        return true;
    }

    void setR_p0_Bounds(float min_R_p0_std, float max_R_p0_std) {
        if (!std::isfinite(min_R_p0_std) || !std::isfinite(max_R_p0_std)) return;
        if (min_R_p0_std <= 0.0f || max_R_p0_std <= min_R_p0_std) return;
        MIN_R_p0_std_ = min_R_p0_std;
        MAX_R_p0_std_ = max_R_p0_std;
    }

    void setR_v0_Bounds(float min_R_v0_std, float max_R_v0_std) {
        if (!std::isfinite(min_R_v0_std) || !std::isfinite(max_R_v0_std)) return;
        if (min_R_v0_std <= 0.0f || max_R_v0_std <= min_R_v0_std) return;
        MIN_R_v0_std_ = min_R_v0_std;
        MAX_R_v0_std_ = max_R_v0_std;
    }

    // Smoothing-horizon multipliers for the two drift-correction channels.
    // The EMA time constant of each channel is mult * tau_target, so the
    // horizon follows the sea state instead of being pinned to one second
    // count.  r_p0 ~ sigma_aw * tau^2 and r_v0 ~ sigma_aw * tau amplify a tau
    // error by different powers, so the two channels need not share a
    // multiplier.
    void setR_p0_AdaptMult(float m) {
        if (std::isfinite(m) && m > 0.0f) adapt_R_p0_mult_ = m;
    }

    void setR_v0_AdaptMult(float m) {
        if (std::isfinite(m) && m > 0.0f) adapt_R_v0_mult_ = m;
    }

    // Size, in natural-log units of the target-to-applied ratio, of a
    // discrepancy the smoothers should treat as a real sea-state move rather
    // than tuner jitter.  Both channels are driven by the same tau and sigma
    // estimates, so their discrepancies move together and share one threshold.
    // Zero or negative leaves the plain proportional horizon.  See
    // seastate::common::adaptiveSmoothingHorizonSec.
    void setR_AdaptSlewLog(float d) {
        if (std::isfinite(d)) adapt_R_slew_log_ = d;
    }

    float getR_p0_AdaptMult() const noexcept { return adapt_R_p0_mult_; }
    float getR_v0_AdaptMult() const noexcept { return adapt_R_v0_mult_; }
    float getR_AdaptSlewLog() const noexcept { return adapt_R_slew_log_; }

    void setMagDelaySec(float delay_sec) {
        if (std::isfinite(delay_sec) && delay_sec >= 0.0f) mag_delay_sec_ = delay_sec;
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
    // the bias state has a long correlation time, so a wrong value learned
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

    void setNominalRaccStd(const Eigen::Vector3f& r) { Racc_nominal_std_ = r; }

    inline float getTauApplied()       const noexcept { return mekf_ ? mekf_->get_aw_time_constant() : NAN; }
    inline float getSigmaApplied()     const noexcept { return mekf_ ? mekf_->get_aw_stationary_std().z() : NAN; }
    inline float getR_p0_std_applied() const noexcept { return mekf_ ? mekf_->get_Rp0_noise_std().z() : NAN; }
    inline float getR_v0_std_applied() const noexcept { return mekf_ ? mekf_->get_Rv0_noise_std().z() : NAN; }
    inline float getR_p0_std_target()  const noexcept { return R_p0_std_target_; }
    inline float getR_v0_std_target()  const noexcept { return R_v0_std_target_; }

    inline float getHeaveAbs() const noexcept {
        if (!mekf_) return NAN;
        return std::fabs(mekf_->get_position().z());
    }

    inline float getDisplacementScale(bool smoothed = true) const noexcept {
        const float tau   = smoothed ? tune_.tau_applied   : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;
        if (!std::isfinite(sigma) || !std::isfinite(tau)) return NAN;
        constexpr float C_HS = 2.0f * std::numbers::sqrt2_v<float> / (std::numbers::pi_v<float> * std::numbers::pi_v<float>);
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

    // Time constant of the still-water attenuation of the wave variance.
    // OU-II keeps it as a knob at 5 s; OU-III fixes 1 s.
    void setSigmaStillnessDecaySec(float sec) {
        if (std::isfinite(sec) && sec > 0.0f) sigma_stillness_decay_sec_ = sec;
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
            apply_R_p0_tune_();
            apply_R_v0_tune_();
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
        // Commit the pseudo-update cadence with the same applied tau so
        // T_S/tau stays constant apart from explicit safety clamps.
        apply_pseudo_update_cadence_(mekf_.get(), tune_.tau_applied);

        const float sigma_floor = std::max(0.05f, band_noise_floor_sigma_());
        const float sZ = std::max(sigma_floor, tune_.sigma_applied);
        const float sH = sZ * P_factor_;
        const Eigen::Vector3f aw_std(sH, sH, sZ);
        mekf_->set_aw_stationary_std(aw_std);
        if (sync_covariance) {
            mekf_->synchronize_aw_covariance_to_stationary();
            last_aw_cov_sync_sec_ = time_;
        }
    }

    // Re-align the posterior a_w marginal with the stationary prior at the
    // adaptation cadence. Runs independently of the tuner so that
    // fixed-tuning modes apply the same policy and remain matched controls.
    void periodic_aw_cov_sync_tick_() {
        if (!mekf_ || !periodicAwCovSyncDue_()) return;
        mekf_->synchronize_aw_covariance_to_stationary();
        last_aw_cov_sync_sec_ = time_;
    }

    // set_Rp0_noise_std()/set_Rv0_noise_std() accept standard deviations, so one
    // pseudo update has covariance r^2.  With updates every T_S seconds, the
    // continuous-equivalent information rate is proportional to 1/(r^2 T_S).
    // The Empirical law preserves the nominal 15 ms information rate by
    // normalizing the filter-input standard deviations:
    //     r_filter = r_base * sqrt(T_0/T_S).
    // The base tuner values remain clamped to their configured bounds.  Do not
    // clamp again after this normalization: the smallest-sea operating point
    // may already sit on a base floor and must be allowed below it when
    // T_S > 15 ms, or the clipping the floor exists to avoid comes straight
    // back.  With unclamped T_S proportional to tau this turns the base
    // sigma_aw*tau^2 and sigma_aw*tau schedules into effective
    // sigma_aw*tau^(3/2) and sigma_aw*tau^(1/2) schedules at the filter input.
    float pseudo_update_information_rate_scale_() const noexcept {
        // The PhysicalMSE law already contains the realized T_S, so
        // renormalizing it again would double-count the cadence.
        if (pseudo_law_ != PseudoAdaptationLaw::Empirical) return 1.0f;
        if (!tau_scaled_pseudo_cadence_ || !mekf_) return 1.0f;
        return cadenceRenormalization_(mekf_->get_pseudo_update_period_s());
    }

    // q_eff^(1/10) for the current acceleration-error density, cached because
    // it is constant between setPseudoAccelNoiseDensity() calls.
    void refresh_pseudo_qeff_pow_() noexcept {
        float r_a = pseudo_accel_noise_density_;
        if (!(std::isfinite(r_a) && r_a > 0.0f))
            r_a = R_PSEUDO_ACCEL_NOISE_DENSITY_DEFAULT;
        pseudo_qeff_pow_ = std::pow(2.0f * r_a, 0.1f);
    }

    // The two pseudo targets for the selected law.
    //
    // Empirical returns the *base* pair of Eq. (implemented-base); the caller's
    // cadence renormalization sqrt(T_0/T_S) then carries it to the filter
    // input.  PhysicalMSE returns the filter input directly, because the law is
    // derived on the continuous densities rho = r^2 T_S and therefore already
    // contains the realized cadence; pseudo_update_information_rate_scale_()
    // must not renormalize it a second time.
    void pseudo_targets_from_law_(float tau, float sigma,
                                  float& r_p, float& r_v) const noexcept {
        if (pseudo_law_ == PseudoAdaptationLaw::Empirical) {
            r_p = R_p0_coeff_ * sigma * tau * tau;
            r_v = R_v0_coeff_ * sigma * tau;
            return;
        }

        const float TS = pseudo_update_period_for_(tau);
        const float c_sigma = (std::isfinite(sigma_coeff_) && sigma_coeff_ > 0.0f)
                              ? sigma_coeff_ : 1.0f;
        const float sigma_aB = std::max(sigma / c_sigma, 1e-6f);
        const float ratio = pseudo_mse_ratio_;
        if (!(TS > 0.0f) || !(tau > 0.0f) ||
            !(std::isfinite(ratio) && ratio > 0.0f)) {
            // degenerate config: fall back to the Empirical base pair
            r_p = R_p0_coeff_ * sigma * tau * tau;
            r_v = R_v0_coeff_ * sigma * tau;
            return;
        }
        // sigma^(4/5) tau^(12/5) == (sigma tau^3)^(4/5) exactly, so the
        // position channel needs one transcendental rather than two, and
        // q_eff^(1/10) is constant for a given sensor and is cached.  The
        // velocity channel is then r_p / (ratio tau), which is exact and free:
        // Corollary (optimal pseudo-channel ratio) fixes r_p/r_v = ratio * tau
        // and the two channels share every other factor.
        const float u = sigma_aB * tau * tau * tau;
        r_p = pseudo_mse_coeff_ * pseudo_qeff_pow_
            * std::pow(u, 0.8f) / std::sqrt(TS);
        r_v = r_p / (ratio * tau);
    }

    void apply_R_p0_tune_(float rp_scale = 1.0f) {
        if (!mekf_) return;
        const float p = (std::isfinite(rp_scale) && rp_scale > 0.0f) ? std::min(rp_scale, 1.0f) : 1.0f;
        const float R_p0_base = std::min(std::max(tune_.R_p0_std_applied, MIN_R_p0_std_), MAX_R_p0_std_);
        const float R_p0_b = R_p0_base * pseudo_update_information_rate_scale_();
        const float rp_z = R_p0_b * p;
        mekf_->set_Rp0_noise_std(Eigen::Vector3f(
            rp_z * R_p0_x_factor_,
            rp_z * R_p0_y_factor_,
            rp_z));
    }

    void apply_R_v0_tune_(float rv_scale = 1.0f) {
        if (!mekf_) return;
        const float p = (std::isfinite(rv_scale) && rv_scale > 0.0f) ? std::min(rv_scale, 1.0f) : 1.0f;
        const float R_v0_base = std::min(std::max(tune_.R_v0_std_applied, MIN_R_v0_std_), MAX_R_v0_std_);
        const float R_v0_b = R_v0_base * pseudo_update_information_rate_scale_();
        mekf_->set_Rv0_noise_std(Eigen::Vector3f::Constant(R_v0_b * p));
    }

    void update_tuner(float dt, float a_vertical_measurement, float freq_hz_for_tuner) {
        stepVarianceChannel_(dt, a_vertical_measurement, freq_hz_for_tuner);

        if (!advanceStartupStage_(wavePeriodUsable())) return;

        // This attenuation follows the statistical variance EMA.  Its horizon
        // must preserve wave energy across brief quiet intervals.
        const float sea_time_sec = measureOperatingPoint_(
            frontEndStill_(), frontEndStillTimeSec_(), sigma_stillness_decay_sec_,
            tau_coeff_, sigma_coeff_);

        float R_p0_raw = NAN;
        float R_v0_raw = NAN;
        pseudo_targets_from_law_(tau_target_, sigma_target_, R_p0_raw, R_v0_raw);

        if (enable_clamp_) {
            R_p0_std_target_ = std::min(std::max(R_p0_raw, MIN_R_p0_std_), MAX_R_p0_std_);
            R_v0_std_target_ = std::min(std::max(R_v0_raw, MIN_R_v0_std_), MAX_R_v0_std_);
        } else {
            R_p0_std_target_ = R_p0_raw;
            R_v0_std_target_ = R_v0_raw;
        }

        adapt_mekf(dt, tau_target_, sigma_target_, R_p0_std_target_, R_v0_std_target_,
                   sea_time_sec);
    }

    void adapt_mekf(float dt, float tau_t, float sigma_t, float R_p0_t, float R_v0_t,
                    float sea_time_sec) {
        const float alpha = tauSigmaSmoothingAlpha_(dt, sea_time_sec);

        const float R_p0_sec = seastate::common::adaptiveSmoothingHorizonSec(
            adapt_R_p0_mult_, tau_t, R_p0_t, tune_.R_p0_std_applied,
            adapt_R_slew_log_, dt);
        const float R_v0_sec = seastate::common::adaptiveSmoothingHorizonSec(
            adapt_R_v0_mult_, tau_t, R_v0_t, tune_.R_v0_std_applied,
            adapt_R_slew_log_, dt);
        const float alpha_R_p0 = 1.0f - std::exp(-dt / R_p0_sec);
        const float alpha_R_v0 = 1.0f - std::exp(-dt / R_v0_sec);

        tune_.tau_applied      += alpha      * (tau_t   - tune_.tau_applied);
        tune_.sigma_applied    += alpha      * (sigma_t - tune_.sigma_applied);
        tune_.R_p0_std_applied += alpha_R_p0 * (R_p0_t  - tune_.R_p0_std_applied);
        tune_.R_v0_std_applied += alpha_R_v0 * (R_v0_t  - tune_.R_v0_std_applied);

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
        first_mag_update_time_ = NAN;

        mekf_->set_acc_bias_updates_enabled(false);
    }

    // The accelerometer sigma the stage logic wants, before any vibration
    // inflation.  Returns a zero vector when it is not known, which is the
    // signal to leave the commanded covariance alone.
    Eigen::Vector3f racc_base_std_() const {
        if (Racc_nominal_std_.allFinite() && Racc_nominal_std_.minCoeff() > 0.0f) {
            return Racc_nominal_std_;
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

        // Accel-bias learning stays gated on the magnetometer-update count.
        const bool allow_bias = !accel_bias_locked_ && !acc_bias_hold_;
        mekf_->set_acc_bias_updates_enabled(allow_bias);

        apply_R_p0_tune_();
        apply_R_v0_tune_();
    }

    Eigen::Vector3f Racc_nominal_std_ = Eigen::Vector3f::Constant(0.0f);

    bool accel_bias_locked_ = true;
    int  mag_updates_applied_ = 0;
    static constexpr int MAG_UPDATES_TO_UNLOCK = 250;
    int  mag_updates_to_unlock_ = MAG_UPDATES_TO_UNLOCK;
    bool acc_bias_hold_ = false;

    bool  with_mag_;
    float mag_delay_sec_ = MAG_DELAY_SEC;
    float first_mag_update_time_ = NAN;

    seastate::common::TiltResetWatchdog tilt_watchdog_{};

    float MIN_R_p0_std_           = MIN_R_p0_std;
    float MAX_R_p0_std_           = MAX_R_p0_std;
    float MIN_R_v0_std_           = MIN_R_v0_std;
    float MAX_R_v0_std_           = MAX_R_v0_std;
    float sigma_stillness_decay_sec_ = 5.0f;
    float adapt_R_p0_mult_        = ADAPT_R_p0_MULT;
    float adapt_R_v0_mult_        = ADAPT_R_v0_MULT;
    float adapt_R_slew_log_       = ADAPT_R_SLEW_LOG;


    // Per-axis horizontal p-regularization scale, against the vertical one.
    // In a sweep over
    // the eight scored records and three IMU seed triplets, 3D displacement RMS
    // has an interior minimum at 0.65 -- -3.56 percent against 1, with the same
    // sign in all 24 record x seed cells -- and 0.8 and 0.55 bracket it at
    // -2.95 and -2.80.  Vertical is unchanged at -0.01 percent.  That is the
    // same shape and nearly the same optimum OU-III's own regularizer has, on a
    // different pseudo-measurement, which is the strongest evidence either
    // family offers that the effect is structural rather than a fit to one
    // wrapper.
    //
    // Both axes use 0.72.  This is not the OU-II
    // family's own argmin -- 0.65 is -- but the basin is flat enough that the
    // difference does not matter: measured on the same 24 cells, 0.72 scores
    // -3.50 percent of 3D RMS against 1 where 0.65 scores -3.56, a difference
    // of six hundredths of a point.
    float R_p0_x_factor_ = 0.72f;
    float R_p0_y_factor_ = 0.72f;
    float P_factor_       = 1.5f;

    TuneState tune_;

    float R_p0_std_target_ = NAN;
    float R_v0_std_target_ = NAN;

    // Empirical law: r_p0 = R_p0_coeff * sigma_aw * tau^2 and
    // r_v0 = R_v0_coeff * sigma_aw * tau, before cadence normalization.
    // The time scale is tau = tau_coeff * T_z / 2.  These coefficients apply
    // only when Empirical is selected, not to the deployed PhysicalMSE law.
    float R_p0_coeff_  = 0.65f;
    float R_v0_coeff_  = 1.3f;

    // PhysicalMSE is the deployed law: it is the one of the two whose powers
    // of sigma_a and tau are derived from a physical displacement-MSE
    // criterion rather than selected empirically.  Empirical stays selectable
    // as the low-cost embedded configuration.
    PseudoAdaptationLaw pseudo_law_ = PseudoAdaptationLaw::PhysicalMSE;
    float pseudo_mse_coeff_ = R_PSEUDO_MSE_COEFF_DEFAULT;
    float pseudo_mse_ratio_ = R_PSEUDO_MSE_RATIO_DEFAULT;
    float pseudo_accel_noise_density_ = R_PSEUDO_ACCEL_NOISE_DENSITY_DEFAULT;
    float pseudo_qeff_pow_ =
        std::pow(2.0f * R_PSEUDO_ACCEL_NOISE_DENSITY_DEFAULT, 0.1f);
    float tau_coeff_   = 0.95f;
    float sigma_coeff_ = 0.85f;

    std::unique_ptr<Kalman3D_Wave_OU_II<float>> mekf_;
};

// Deployed OU-II front end: proxy bootstrap, two-stage magnetic acquisition,
// continuous hard iron and displacement output (shared with OU-III; see
// ProxyStartupFusion.h), around the OU-II filter above.
struct SeaStateFusionConfig_OU_II : seastate::common::ProxyStartupFusionConfig {
    // Empirical covariance weight, applied equally to supplied sensor
    // specifications in simulations and firmware. This does not rescale
    // measurements or claim a lower physical sensor-noise bound.  OU-III
    // applies no such weight.
    float gyro_noise_scale = 0.2f;

    // Initial position and velocity pseudo-measurement variances.  The tuner
    // overwrites both at Live; they set the pre-Live values.  The values
    // reproduce the Kalman3D_Wave_OU_II header defaults exactly.
    float R_p0_noise = 1.5f;
    float R_v0_noise = 0.3f;
};

template<TrackerType trackerT>
class SeaStateFusion_OU_II
    : public seastate::common::ProxyStartupFusion<SeaStateFusionFilter_OU_II<trackerT>,
                                                  SeaStateFusionConfig_OU_II>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Config = SeaStateFusionConfig_OU_II;

    void begin(const Config& cfg) {
        this->beginStartup_(cfg, [](SeaStateFusionFilter_OU_II<trackerT>& impl, const Config& c) {
            impl.initialize_ext(c.sigma_a, c.sigma_g * c.gyro_noise_scale, c.sigma_m,
                                c.Pq0, c.Pb0, c.b0,
                                c.R_p0_noise, c.R_v0_noise,
                                c.gravity_magnitude);
        });
    }
};
