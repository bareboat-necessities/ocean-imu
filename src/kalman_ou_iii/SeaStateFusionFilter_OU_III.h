#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy
  Released under the MIT License

  SeaStateFusionFilter_OU_III

  Marine Inertial Navigational System (INS) Filter for IMU

  Combines multiple real-time estimators into a cohesive ocean-state tracker:

    • Quaternion-based attitude and linear motion estimation via Kalman3D_Wave_OU_III

    • Acceleration-band dominant frequency tracking using one of:
          – AranovskiyFreqTracker     (frequency estimator)
          – KalmANFFreqTracker        (adaptive notch / Kalman frequency tracker)
          - PLLFreqTracker            (PLL frequency tracker)
          – SchmittTrigger            (zero-cross event detector)

      This drives the wave-direction demodulator and nothing else.  The OU
      operating point is a property of the *wave* band, which the acceleration
      band does not measure: an ocean acceleration spectrum is the elevation
      spectrum weighted by (2πf)⁴, so its apparent frequency sits well above
      the spectral peak and barely moves as the sea grows.  Tuning therefore
      reads WavePeriodEstimator, and reads it at every instant of the run --
      before that estimator has a value, a fixed wave-band prior stands in
      rather than the tracker.  The tracker has no path into the adaptation.

    • Dual-stage frequency smoothing of the tracker:
          – Fast 1st-order IIR (≈ few s, ~90% step) for demodulation / direction
          – Slow 1st-order IIR (≈ longer s, ~90% step) reported as getPeriodSec()

    • Online auto-tuning of Kalman filter parameters (τ, σₐ, Rₛ) through
      SeaStateAutoTuner.  The OU-III variance channel first applies a
      period-scaled measurement-only wave band, then estimates acceleration
      variance over a horizon of K wave periods, and applies the σₐ·τ³
      regularization law.

  Where
  – τ (tau):  OU process time constant ≈ ½ · T  (half the wave zero-crossing period)
  – σₐ:       Stationary acceleration scale from period-scaled wave-band RMS
  – Rₛ:       Pseudo-measurement noise controlling integral drift suppression
  – Rₛ_xy:    Reduced in X/Y (anisotropic weighting for vertical-dominant seas)

  Adaptive update:  exponential smoothing toward targets over ADAPT_TAU_SEC

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

#include "freq/FirstOrderIIRSmoother.h"
#include "freq/FrequencyTrackerPolicy.h"
#include "tuner/AdaptiveWaveBandPass.h"
#include "tuner/SeaStateAutoTuner.h"
#include "tuner/WavePeriodEstimator.h"
#include "tuner/VerticalAccelComplementary.h"
#include "tuner/MagAutoTuner.h"
#include "tuner/ContinuousMagHardIronEstimator.h"
#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"
#include "wave_dir/KalmanWaveDirection.h"
#include "wave_dir/WaveDirectionDetector.h"
#include "wave_dir/WaveDirectionFrame.h"
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

// Estimated pre-band vertical accel noise floor (1σ), m/s².
// The adaptive sigma band propagates this white-noise variance through its
// actual time-varying coefficients before subtraction.
constexpr float ACC_NOISE_FLOOR_SIGMA_DEFAULT = 0.12f;

constexpr float MIN_FREQ_HZ = 0.2f;
constexpr float MAX_FREQ_HZ = 6.0f;

// Floor for the wave-band tuning frequency.  MIN_FREQ_HZ bounds the
// acceleration-band tracker and is far too high for a zero-crossing period:
// 0.03 Hz admits a 33 s swell.
constexpr float MIN_TUNE_FREQ_HZ = 0.03f;

// Ceiling for the wave-band tuning frequency.  MAX_FREQ_HZ is the tracker's
// bound and does not belong in the adaptation path: with the tuning frequency
// read from the wave band, setFreqBounds() is a wave-direction knob and must
// not move the OU operating point.  Neither bound binds on any reference
// record -- the wave band spans 0.12 to 0.40 Hz -- so this is a safety limit,
// not a tuning surface.  1.5 Hz is a 0.67 s zero-crossing period, shorter than
// any sea a hull responds to.
constexpr float MAX_TUNE_FREQ_HZ = 1.5f;

// Wave-band tuning frequency used before WavePeriodEstimator has a value.
//
// It has to be a constant: the whole point of the wave-band source is that no
// part of the adaptation path reads the acceleration-band tracker, and a
// constant is trivially exogenous.  0.2 Hz is a 5 s zero-crossing period, and
// the estimator reports 2.3-8.4 s across the reference family, so the prior is
// never worse than a factor of two off; the estimator replaces it after about
// 50 s in any case.  It is the same constant SeaStateFusionFilter_TFG has always
// used for this, which is where the value comes from rather than a fit.
// Sensitivity to it is measured in docs/ou-sigma-horizon.md: sweeping it over
// 0.1-0.4 Hz leaves every scored 900 s metric unchanged to four decimal places
// in both families, because the estimator has replaced it 250 s before the
// window opens.
constexpr float TUNE_FREQ_PRIOR_HZ = 0.2f;

// JONSWAP-similar acceleration-variance band.  Away from the absolute safety
// clamps, its transfer shape is fixed in f/f_tune.
constexpr float SIGMA_BAND_LOW_RATIO_DEFAULT  = 0.5f;
constexpr float SIGMA_BAND_HIGH_RATIO_DEFAULT = 4.0f;
constexpr float SIGMA_BAND_MIN_HZ_DEFAULT     = 0.01f;
constexpr float SIGMA_BAND_MAX_HZ_DEFAULT     = 6.0f;

constexpr float MIN_TAU_S   = 0.02f;
// tau now scales with the zero-crossing wave period rather than with an
// acceleration-band frequency, so the ceiling has to admit a developed sea:
// T_z reaches 8.6 s at H_s = 8.5 m and a long swell goes further.  The old
// 3.0 s ceiling was reached at H_s = 8.5 m, which clipped the operating point
// exactly where the filter was losing.
constexpr float MAX_TAU_S   = 12.0f;
constexpr float MAX_SIGMA_A = 6.0f;
// The old 0.4 floor was not a safety limit in practice, it was the binding
// constraint on every low-motion sea.  The schedule asks for 0.24 m*s at the
// calibrated H_s = 0.27 m point, so the floor clipped it and pinned both
// low-motion scenarios at an error the tuner multipliers could not move: a
// full sweep of c_tau, c_sigma and c_R left J0.27 and P0.27 constant to three
// decimals.  Dropping the floor below the schedule's own demand recovers
// -8.3 % on the worst sea, -9.5 % on PM-Stokes 0.27 m and -2.5 % on the
// eight-sea mean, and cuts the near-still H_s = 0.05 m stress case from
// 27.0 to 17.6 percent of H_s (13.5 mm to 8.8 mm absolute).  The gain saturates below ~0.25; 0.15 keeps a
// real guard against r_S collapsing toward zero while staying clear of the
// calibrated envelope.
constexpr float MIN_R_S     = 0.15f;
// r_S ~ sqrt(R_a) tau^3 inherits that range.  The old 35 m*s ceiling was the
// binding constraint at H_s = 8.5 m: the calibrated fixed-oracle point sat at
// 34.66 and the error was still falling monotonically against it.
constexpr float MAX_R_S     = 400.0f;

// Legacy fixed-second horizon is retained for ablation/backward compatibility.
// A positive sea-period multiplier makes the common tau/sigma EMA self-similar:
//     tau_ema = ADAPT_TAU_SEA_PERIODS * T_sea,  T_sea = T_z/2.
// The measured default 0.40 is the 0.20*T_z winner from the paired sweep.
// Zero selects the legacy ADAPT_TAU_SEC path for controlled ablations.
constexpr float ADAPT_TAU_SEC              = 1.8f;
constexpr float ADAPT_TAU_SEA_PERIODS          = 0.40f;
constexpr float TUNER_FREQ_EMA_SEA_PERIODS     = 0.10f;
// This cadence belongs only to posterior a_w covariance maintenance.  Online
// tuner candidates are staged after every valid physical measurement.
constexpr float ADAPT_EVERY_SECS           = 0.1f;
// Smoothing horizon of the r_S channel, in units of tau_target.  Measured on
// the versioned records against synthesized sea-state transitions: the error
// during a transition falls monotonically as this shortens while the worst
// single realization degrades monotonically the other way, and the value is
// where the mean gain and the worst-record loss cross.  That crossing sat at
// 3.0 while the transition instrument was a 360 s crossfade; re-measured on
// the 120 s crossfade the study now uses -- fast enough for the horizon to
// actually lag it -- it sits at 1.5, and 1.5 also wins on the stationary
// ensemble, on attitude and on the accelerometer bias.  r_S ~ tau^3 amplifies
// tau noise by the third power, which is why this stays above the horizon a
// tau^1 channel would want.  See docs/ou-ema-adaptation-tuning.md.
constexpr float ADAPT_RS_MULT              = 1.5f;   // dimensionless
// Discrepancy, in natural-log units of the r_S target-to-applied ratio, above
// which the smoothing horizon shortens.  Zero keeps the plain proportional
// horizon; see docs/ou-ema-adaptation-tuning.md.
constexpr float ADAPT_RS_SLEW_LOG          = 0.0f;   // ln units
constexpr float ONLINE_TUNE_WARMUP_SEC     = 5.0f;
constexpr float MAG_DELAY_SEC              = 7.0f;

// Gains of the private Mahony observer, which levels the vertical channel and
// solves the startup attitude.
//
// two_kp is the accelerometer-to-gyro correction corner and carries the same
// constraint it does in the vertical channel: it must stay an order of
// magnitude below the wave band, or the observer levels itself against the
// orbital specific force instead of gravity.
//
// two_ki estimates the gyro bias.  The vertical channel ran at zero for years
// because everything downstream of it is high-passed, but an attitude seed
// keeps whatever static tilt the bias leaves, so the observer that serves both
// has to estimate it; see the member declaration for the measurement.
constexpr float STARTUP_PROXY_TWO_KP_DEFAULT = 0.2f;
constexpr float STARTUP_PROXY_TWO_KI_DEFAULT = 0.02f;

// Frequency smoother dt (SeaStateFusionFilter_OU_III is designed for 200 Hz)
constexpr float FREQ_SMOOTHER_DT = 1.0f / 200.0f;

// Self-similar S=0 pseudo-measurement cadence.  The deployed wrapper historically
// used T_S=15 ms while its initial applied OU time constant is tau=1.1 s, so
// this ratio preserves that exact operating point and scales cadence thereafter.
constexpr float PSEUDO_UPDATE_PERIOD_NOMINAL_S = 0.015f;
constexpr float PSEUDO_UPDATE_TAU_NOMINAL_S = 1.1f;
constexpr float PSEUDO_UPDATE_TAU_RATIO_DEFAULT =
    PSEUDO_UPDATE_PERIOD_NOMINAL_S / PSEUDO_UPDATE_TAU_NOMINAL_S;
// A pseudo update cannot occur more often than the nominal 200 Hz IMU schedule.
// The upper guard is inactive over the current tau <= 12 s operating envelope.
constexpr float PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT = FREQ_SMOOTHER_DT;
constexpr float PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT = 0.25f;

// Integral-regularizer adaptation laws.  All three place the drift-band
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
//   PosteriorRiccati the full transition law, reducing to a sigma_aw tau^3
//                    schedule as zeta -> 0 and to StrongRiccati as
//                    zeta -> infinity.
//
// The deployed envelope has zeta ~ 1e5..1e7 against the bench sensor floor, so
// the strong-observation branch q_eff ~ 2 r_a is the applicable one and the
// drift-driving error is set by the *sensor*, not by the sea.  The base
// schedule takes its acceleration scale accordingly: R_a is the accelerometer
// measurement-noise variance, and sqrt(R_a) tau^3 carries the required units of
// m*s.  The wave amplitude still enters the filter through the OU prior,
// sigma_aw = c_sigma sigma_a,B; it no longer sets pseudo-measurement strength.
//
// Writing the base this way makes C_R a pole placement rather than a gain.
// With r_a = R_a h, the cadence-normalized base is
//     C_R sqrt(R_a) tau^3 sqrt(T_0/T_S) = sqrt(2 r_a) tau^3 / (kappa^3 sqrt(T_S))
// exactly when C_R = sqrt(2h/T_0)/kappa^3, so C_R and the normalized corner
// kappa = omega_R tau are two spellings of one number and the Cubic and
// StrongRiccati laws coincide there.  R_S_COEFF_ANALYTICAL_REFERENCE is that
// value for kappa = 0.3627.  The applied C_R comes from a complete-MEKF sweep
// against it, because the scalar reduction omits attitude/gravity leakage,
// residual bias, three-axis covariance coupling and the cadence clamps; see
// tools/ou_rs_amplitude_retune_sweep.py and
// docs/ou-iii-rs-amplitude-retune.md.
//
// The amplitude tilt of the Riccati laws spans the removed multiplier as a
// one-parameter family: p = 0 is the deployed schedule and p = 1 restores the
// sigma_aw factor, so tools/ou_rs_law_ablation.py still measures exactly that
// one degree of freedom.  The non-Cubic laws exist for that ablation, not as
// candidate defaults.
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
// The ratio between the tilted (p = 1) and deployed (p = 0) members of the
// Riccati family is independent of tau,
//     r_S,p=1 / r_S,p=0 = sigma_aw / sigma_ref,
//     sigma_ref = sqrt(2 r_a) / (kappa^3 c_R sqrt(T_0)),
// because the tau^3/sqrt(T_S) factor is common to both.  sigma_ref is the
// amplitude at which the two members cross, so the ablation measures only the
// sigma_aw dependence of the schedule rather than an overall change of
// regularizer gain.  kappa was calibrated when the deployed base still carried
// sigma_aw, so that sigma_ref matched the sigma_aw of the Hs = 1.5 m nominal
// sea; it is left as it was, and the anchor now sits wherever the re-fitted
// c_R puts it.  Nothing in the ablation depends on where that anchor is.
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
// against, not as the applied value; see docs/ou-iii-rs-amplitude-retune.md.
constexpr float R_S_COEFF_ANALYTICAL_REFERENCE = 17.112f;

// Coefficient C_J of the SpectralMSE law,
//     r_S = C_J q_eff^(1/14) sigma_a,B^(6/7) tau^(24/7) / sqrt(T_S).
// It absorbs the dimensionless spectral moment int x^-8 Phi_a(x) dx of the
// self-similar acceleration spectrum, raised to 3/7.  Evaluating that moment on
// the eight reference spectra and solving the exact balance record by record
// gives C_J ~ 0.054, so the default is an analytical prediction rather than a
// fitted number; the sweep is what decides whether it is the complete-MEKF
// optimum.  See docs/ou-iii-rs-amplitude-retune.md.
constexpr float R_S_MSE_COEFF_DEFAULT = 0.0538f;

struct TuneState {
    float tau_applied   = 1.1f;    // s
    float sigma_applied = 1e-2f;   // m/s²
    float RS_applied    = 0.5f;    // m*s
};

//  Unified SeaState fusion filter
template<TrackerType trackerT>
class SeaStateFusionFilter_OU_III {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using TrackingPolicy = TrackerPolicy<trackerT>;

    enum class StartupStage {
        Cold,        // just booted or just had a big tilt reset
        TunerWarm,   // MEKF + freq running, tuner collecting stats
        TunerReady,  // tuner trusted, MEKF still held by an external bootstrap
        Live         // tuner is trusted; full adaptation & extras allowed
    };

    // Who solves the attitude the filter starts from.
    //
    // StagedMekf is the original behaviour: the MEKF is fed from the first
    // sample and learns its own tilt while running degraded -- linear block
    // off, accelerometer bias frozen, Racc inflated -- and the caller reads
    // *that* attitude back to gate the magnetometer and to build the frame the
    // magnetic world reference is averaged in.  Those reads are the problem.
    // The reference and the yaw gauge are locked once, so whatever tilt error
    // the warming MEKF has at that moment becomes a standing attitude and
    // heading bias, and it is exactly the moment the MEKF is least able to
    // level: the linear block that absorbs orbital acceleration is switched
    // off, so its accelerometer update is fighting the waves with no state to
    // put them in.
    //
    // MahonyProxy takes those jobs away from the MEKF and gives them to the
    // private Mahony observer this filter already runs (see
    // VerticalAccelComplementary).  That observer is a pure function of the raw
    // gyro and accelerometer, and its correction corner sits an order of
    // magnitude below the wave band, so it rejects the orbital specific force
    // by construction rather than chasing it.  The measurement-only front end
    // -- proxy, frequency tracker, wave-period estimator, sigma band, tuner --
    // runs from the first sample without the MEKF, so by the time the tilt has
    // settled and the magnetometer has gauged north, the operating point is
    // converged too.  The MEKF is then seeded with that attitude and goes Live
    // directly, never occupying the degraded warmup configuration at all.
    //
    // See docs/ou-iii-startup-init.md for the measured comparison, including
    // the one metric that moves the wrong way.
    enum class StartupInitPolicy {
        StagedMekf,   // legacy: MEKF learns its own tilt while warming
        MahonyProxy,  // default: proxy owns tilt + mag learning, MEKF starts Live
    };

    explicit SeaStateFusionFilter_OU_III(bool with_mag = true)
        : with_mag_(with_mag),
          time_(0.0),
          last_adapt_time_sec_(0.0),
          freq_hz_(FREQ_GUESS),
          freq_hz_slow_(FREQ_GUESS)
    {
        // Default cutoff ~max_freq_hz_ Hz: passes waves, kills 8–37 Hz engine band
        freq_input_lpf_.setCutoff(max_freq_hz_);
        freq_stillness_.setTargetFreqHz(min_freq_hz_);
        // OU-III opts into the measured sea-scaled tuner-frequency EMA.
        // SeaStateAutoTuner itself remains fixed-seconds by default so OU-II/TFG
        // are unchanged by this OU-III-only calibration.
        tuner_.setFrequencySmoothingSeaPeriods(TUNER_FREQ_EMA_SEA_PERIODS);
        startup_stage_   = StartupStage::Cold;
        startup_stage_t_ = 0.0f;
    }

    StartupStage getStartupStage() const noexcept { return startup_stage_; }
    bool isAdaptiveLive() const noexcept { return startup_stage_ == StartupStage::Live; }

    // The operating point is trustworthy.  Under StagedMekf this coincides
    // with going Live; under MahonyProxy it is reached first, while the MEKF
    // is still held, and it is one of the conditions the caller waits on
    // before handing the attitude over.
    bool isTunerReady() const noexcept {
        return startup_stage_ == StartupStage::TunerReady ||
               startup_stage_ == StartupStage::Live;
    }

    void setStartupInitPolicy(StartupInitPolicy policy) noexcept {
        startup_init_policy_ = policy;
    }
    StartupInitPolicy startupInitPolicy() const noexcept {
        return startup_init_policy_;
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

    // Time update (IMU integration + frequency tracking)
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

    // Attitude of the startup Mahony observer, BODY -> NED.  Only its tilt is
    // meaningful; yaw is unobservable to it and drifts.
    Eigen::Quaternionf startupProxyQuat() const noexcept {
        return vertical_accel_comp_.quaternion();
    }

    // Tilt-only form of the same attitude, safe to use as a magnetometer
    // accumulation frame because no heading can leak through it.
    Eigen::Quaternionf startupProxyTiltQuat() const noexcept {
        return vertical_accel_comp_.tiltQuaternion();
    }

    bool startupProxyInitialized() const noexcept {
        return vertical_accel_comp_.isInitialized();
    }

    // Gains of the private Mahony observer, which serves both the vertical
    // channel and the startup attitude.  Same knob as
    // setWavePeriodComplementaryGains(); kept under this name because the
    // startup path is the one with an opinion about two_ki.
    void setStartupProxyGains(float two_kp, float two_ki) {
        vertical_accel_comp_.setGains(two_kp, two_ki);
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

        // Private Mahony observer for the wave-period estimator and default
        // sigma channel. It is fed raw gyro and accelerometer before the MEKF
        // sees them, so its levelled vertical acceleration remains a pure
        // function of the measurements.
        vertical_accel_comp_.update(dt, gyro, acc, g_std);

        // MEKF updates first (attitude + latent a_w)
        if (drive_mekf) {
            mekf_->time_update(gyro, dt);
            mekf_->measurement_update_acc_only(acc, tempC);
        }

        // The tilt watchdog reads and rewrites MEKF attitude, so it only has
        // meaning while the MEKF is the one propagating it.  A bootstrap that
        // has not handed over yet has no attitude here to run away.
        if (drive_mekf) {
            Eigen::Quaternionf q_bw = mekf_->quaternion_boat();
            q_bw.normalize();

            const Eigen::Vector3f z_body_down_world = q_bw * Eigen::Vector3f(0.0f, 0.0f, 1.0f);
            const Eigen::Vector3f z_world_down(0.0f, 0.0f, 1.0f);

            float cos_tilt = z_body_down_world.normalized().dot(z_world_down);
            cos_tilt = std::max(-1.0f, std::min(1.0f, cos_tilt));
            const float tilt_deg = std::acos(cos_tilt) * 57.295779513f;

            constexpr float TILT_RESET_DEG = 70.0f;
            constexpr float TILT_RESET_HOLD_SEC = 0.35f;
            constexpr float TILT_RESET_COOLDOWN_SEC = 3.0f;

            if (tilt_reset_cooldown_sec_ > 0.0f) {
                tilt_reset_cooldown_sec_ = std::max(0.0f, tilt_reset_cooldown_sec_ - dt);
            }

            if (tilt_deg > TILT_RESET_DEG) {
                tilt_over_limit_sec_ += dt;
            } else {
                // decay quickly on recovery so brief transients do not trigger resets
                tilt_over_limit_sec_ = std::max(0.0f, tilt_over_limit_sec_ - 2.0f * dt);
            }

            if (tilt_over_limit_sec_ >= TILT_RESET_HOLD_SEC && tilt_reset_cooldown_sec_ <= 0.0f) {
                if (startup_stage_ == StartupStage::Live) {
                    // In Live, re-lock only tilt while preserving yaw/north frame.
                    mekf_->initialize_from_acc_preserve_yaw(acc);
                } else {
                    // During startup stages, accel-only re-lock is acceptable.
                    mekf_->initialize_from_acc(acc);
                    enterCold_();
                    resetTrackingState_();
                }

                tilt_over_limit_sec_ = 0.0f;
                tilt_reset_cooldown_sec_ = TILT_RESET_COOLDOWN_SEC;
            }
        }

        // The one levelled vertical measurement every consumer in this filter
        // reads: the frequency tracker and its stillness detector, the
        // wave-period estimator, and the sigma channel.  It comes from the
        // private Mahony observer and therefore reads no MEKF state.  The
        // observer is seeded from the first accelerometer sample, so the value
        // is usable immediately; its isReady() gate is about settled period
        // *statistics*, which is a stricter requirement than a usable tilt.
        const float a_vert_measurement = vertical_accel_comp_.verticalAccelUpMs2();

        // LPF on the tracker input
        const float a_vert_lp = freq_input_lpf_.step(a_vert_measurement, dt);

        // Raw freq from tracker
        const float f_tracker = static_cast<float>(tracker_policy_.run(a_vert_lp, dt));
        f_raw = f_tracker;

        // Stillness detector shares the tracker's input, as it always has.
        const float f_after_still = freq_stillness_.step(a_vert_lp, dt, f_tracker);

        // Fast & slow smoothed frequencies
        float f_fast = freq_fast_smoother_.update(f_after_still);
        float f_slow = freq_slow_smoother_.update(f_fast);

        f_fast = std::min(std::max(f_fast, min_freq_hz_), max_freq_hz_);
        f_slow = std::min(std::max(f_slow, min_freq_hz_), max_freq_hz_);

        freq_hz_      = f_fast;   // demod / direction
        freq_hz_slow_ = f_slow;   // tuner / moments

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

        const float omega = 2.0f * static_cast<float>(M_PI) * freq_hz_;

        // Resolve direction in a leveled frame aligned with boat heading.
        // This removes roll/pitch mixing while preserving 0 deg = bow and
        // positive angles toward starboard.  Direction may use the MEKF frame;
        // the default tuner channels do not.
        //
        // Before handoff the MEKF has no attitude to offer, so the proxy's is
        // used.  heading_frame_acceleration() resolves into the projected bow
        // axis and is therefore invariant under q -> Rz(psi) q, so only the
        // tilt of whichever quaternion is supplied can reach the result.
        const auto direction_accel = wave_direction::heading_frame_acceleration<float>(
            drive_mekf ? mekf_->quaternion_boat() : startupProxyQuat(), acc, g_std);

        // Stage 1 estimates the apparent propagation plane as an unsigned axis
        // relative to boat heading.  Stage 2 resolves propagation sense along
        // that same axis from horizontal/vertical orbital phase.
        // Zero-crossing wave period.  This runs beside the frequency tracker
        // rather than replacing it: the tracker supplies the acceleration-band
        // carrier the direction demodulator needs, while the OU operating point
        // needs the wave band.
        //
        // The input must be levelled and must not read estimator state, and
        // those two requirements pulled against each other for a while.
        //
        // Levelled, because double integration weights a spectrum by
        // 1/omega^4, so the sub-band gravity leakage a tilting platform puts
        // into a raw body-Z residual dominates the elevation proxy.  Fed that
        // residual the estimator reports 6.8-10.0 s whatever the sea does,
        // against a truth of 2.4-8.7 s, and vertical RMS degrades 2.5x.  That
        // is why no body-Z path is left in this filter.
        //
        // Exogenous, because levelling with the filter's own attitude closes a
        // loop: a 0.25 rad attitude displacement moved the reported period
        // 8.05 -> 10.1 s and tau by 1.25x, and it reached the linear block too,
        // since displacing v, p, S or a_w perturbs attitude through the
        // filter's cross-covariances.  The stability appendix carried that as
        // its open interconnection.
        //
        // VerticalAccelComplementary satisfies both: it levels, but with a
        // private Mahony observer reading only the raw gyro and accelerometer,
        // so it is a pure function of the measurements.  It costs nothing --
        // over the eight reference records it matches the old attitude-levelled
        // input to within 0.2% of vertical RMS -- and it is the default.
        // setWavePeriodInput() still reaches the attitude-levelled one for
        // ablation;
        // tests/kalman_ou_iii/tuner_coupling-test.cpp asserts the default is
        // exogenous bit-for-bit and bounds the Leveled path's gain.
        wave_period_.update(dt, wave_period_input_ms2_(direction_accel));

        dir_filter_.update(direction_accel.forward_ms2,
                           direction_accel.starboard_ms2,
                           omega, dt);
        const Eigen::Vector2f propagation_axis_boat = dir_filter_.getAxis();
        dir_sign_state_ = dir_sign_.update(
            direction_accel.forward_ms2,
            direction_accel.starboard_ms2,
            direction_accel.up_ms2,
            propagation_axis_boat.x(), propagation_axis_boat.y(),
            dt, dir_filter_.getLastStableConfidence());
    }

public:
    //  Magnetometer correction
    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!with_mag_ || !mekf_) return;
        if (time_ < mag_delay_sec_) return;

        mekf_->measurement_update_mag_only(mag_body_ned);
        mag_updates_applied_++;

        if (!std::isfinite(first_mag_update_time_)) {
            first_mag_update_time_ = static_cast<float>(time_);
        }

        // We can "unlock" once mag has had a few updates, but we DO NOT
        // enable accel-bias learning or restore Racc unless we're already Live.
        if (accel_bias_locked_ &&
            startup_stage_ == StartupStage::Live &&
            mag_updates_applied_ >= mag_updates_to_unlock_ &&
            std::isfinite(first_mag_update_time_) &&
            (static_cast<float>(time_) - first_mag_update_time_) > 1.0f) // 1s guard
        {
            accel_bias_locked_ = false;

            // Only allow accel bias to start learning once the system is Live
            // and no external hold is in force.
            if (freeze_acc_bias_until_live_ && startup_stage_ == StartupStage::Live &&
                !acc_bias_hold_) {
                mekf_->set_acc_bias_updates_enabled(true);

                // Restore nominal Racc only when bias learning is allowed.
                if (warmup_Racc_active_) {
                    if (Racc_nominal_.allFinite() && Racc_nominal_.maxCoeff() > 0.0f) {
                        mekf_->set_Racc_std(Racc_nominal_);
                        warmup_Racc_active_ = false;
                    }
                }
            }
        }
    }

    void setWithMag(bool with_mag) {
        with_mag_ = with_mag;
    }

    // Anisotropy configuration (runtime)
    // S-factor scales horizontal vs vertical stationary std of a_w.
    // RS XY factor scales pseudo-measurement noise in X/Y vs Z.
    void setSFactor(float s) {
        if (std::isfinite(s) && s > 0.0f) {
            S_factor_ = s;
        }
    }
    // The upper bound is 4, not 1.  A ceiling of 1 encoded the assumption that
    // the horizontal integral anchor can only ever be tighter than the vertical
    // one, which silently excluded the value the similarity law asks for: with
    // sigma_aw,H = S_factor * sigma_aw,Z the natural scale of the horizontal S
    // state is S_factor times the vertical one (Theorem "Nondimensional
    // similarity of the OU-III chain", sigma_S ~ sigma_aw tau^3), so the
    // dimensionless regularizer r_S/(sigma_aw tau^3) is equal on all three axes
    // only at R_S_xy_factor = S_factor > 1.  Under the old clamp that
    // configuration was accepted and then silently reduced to 1, so an ablation
    // of it measured nothing and looked like a null result.
    void setRSXYFactor(float k) {
        if (std::isfinite(k)) {
            R_S_xy_factor_ = std::min(std::max(k, 0.0f), 4.0f);
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

                if (enable_linear_block_) {
                    apply_RS_tune_();
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

    // Dimensionless sigma-band shape.  Defaults correspond to the theorem's
    // measurement-only soft band [0.5,4] in units of f_tune.  Changing these
    // ratios preserves the similarity structure as long as they remain fixed
    // across sea states; absolute limits below are safety clamps only.
    void setSigmaWaveBandRatios(float low_ratio, float high_ratio) {
        sigma_wave_band_.setRatios(low_ratio, high_ratio);
    }
    void setSigmaWaveBandLimitsHz(float min_hz, float max_hz) {
        sigma_wave_band_.setLimitsHz(min_hz, max_hz);
    }
    float getSigmaWaveBandLowHz() const noexcept { return sigma_wave_band_.lowHz(); }
    float getSigmaWaveBandHighHz() const noexcept { return sigma_wave_band_.highHz(); }
    float getSigmaWaveBandLowRatio() const noexcept { return sigma_wave_band_.lowRatio(); }
    float getSigmaWaveBandHighRatio() const noexcept { return sigma_wave_band_.highRatio(); }
    float getSigmaBandNoiseStd() const noexcept { return band_noise_floor_sigma_(); }

    // Configure LPF on BODY-Z proxy for tracker input
    void setFreqInputCutoffHz(float fc) {
        freq_input_lpf_.setCutoff(fc);
    }

    void enableClamp(bool flag = true) {
        enable_clamp_ = flag;
    }
    void enableTuner(bool flag = true) {
        enable_tuner_ = flag;
        if (!flag) online_tune_apply_pending_ = false;
    }
    bool tunerEnabled() const noexcept { return enable_tuner_; }

    // Policy for the latent-acceleration marginal P_{a_w a_w}.
    //
    // Default (true): once per adaptation period the marginal is re-aligned
    // with the stationary OU covariance, keeping the cross-covariances the
    // filter has learned.  This is a deliberate bounded covariance inflation.
    // It stops the a_w marginal from settling far below the level the process
    // model considers stationary, which keeps the accelerometer gain
    // responsive when the sea state changes.  It is not free -- it discards
    // posterior information at the adaptation cadence -- so the alternative
    // is available and measured rather than assumed.
    //
    // With false, the marginal is aligned only at discrete reconfiguration
    // events (construction, the transition to Live, setFixedTuning()) and a
    // changed stationary scale reaches the filter solely through the discrete
    // OU process covariance.
    //
    // The policy is applied independently of the tuner so that fixed-tuning
    // modes run it too.  Otherwise an adaptive-versus-fixed comparison would
    // confound whether the parameters adapt with whether part of the
    // covariance is periodically re-aligned.
    void setPeriodicAwCovarianceSync(bool flag) {
        periodic_aw_cov_sync_ = flag;
        last_aw_cov_sync_sec_ = time_;
    }
    bool periodicAwCovarianceSync() const noexcept { return periodic_aw_cov_sync_; }

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
    // Enabled by default; disabling restores the historical fixed 15 ms cadence
    // for direct old-versus-new ablation. Whenever cadence changes while Live,
    // reapply R_S so its per-update covariance stays information-rate matched.
    void setTauScaledPseudoUpdateCadence(bool flag) {
        tau_scaled_pseudo_cadence_ = flag;
        if (!mekf_) return;
        if (flag) apply_pseudo_update_cadence_();
        else mekf_->set_pseudo_update_period_s(pseudo_update_fixed_period_s_);
        if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
            apply_RS_tune_();
        }
    }
    bool tauScaledPseudoUpdateCadence() const noexcept { return tau_scaled_pseudo_cadence_; }
    void setPseudoUpdateTauRatio(float ratio) {
        if (!(std::isfinite(ratio) && ratio > 0.0f)) return;
        pseudo_update_tau_ratio_ = ratio;
        if (tau_scaled_pseudo_cadence_) {
            apply_pseudo_update_cadence_();
            if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
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
            apply_pseudo_update_cadence_();
            if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
                apply_RS_tune_();
            }
        }
    }
    // Select which drift-band regularization law generates the r_S target.
    // Changing it while Live restages the schedule on the next adapt tick; the
    // active covariance is refreshed here so ablations switch cleanly.
    void setRSLaw(RSAdaptationLaw law) {
        rs_law_ = law;
        if (mekf_ && startup_stage_ == StartupStage::Live && enable_linear_block_) {
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

    float getPseudoUpdateTauRatio() const noexcept { return pseudo_update_tau_ratio_; }
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
        if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
            apply_RS_tune_();
        }
        return true;
    }

    // Freeze one adaptation channel while the other keeps tracking the sea.
    //
    // The deployed law derives r_S from tau: r_S = clip(C_R sqrt(R_a) tau^3), so
    // simply freezing tau and sigma_a with setFixedTuning() freezes r_S as
    // well, and an ablation built that way cannot say which channel carries
    // the benefit.  Here the tuner keeps running and keeps deriving r_S from
    // its *live* tau and sigma_a estimates; only the channels named below are
    // held at the supplied operating point on the way to the filter.
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
            if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
                apply_RS_tune_();
            }
        }
        return true;
    }

    bool frozenOUChannel() const noexcept { return freeze_ou_channel_; }
    bool frozenRSChannel() const noexcept { return freeze_RS_channel_; }

    // Enable/disable use of the extended linear block [v,p,S,a_w] in Kalman3D_Wave_OU_III.
    void enableLinearBlock(bool flag = true) {
        enable_linear_block_ = flag;
        if (mekf_) {
            const bool on_now = flag && (startup_stage_ == StartupStage::Live);
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

    void setRSBounds(float min_RS, float max_RS) {
        if (!std::isfinite(min_RS) || !std::isfinite(max_RS)) return;
        if (min_RS <= 0.0f || max_RS <= min_RS) return;
        min_R_S_ = min_RS;
        max_R_S_ = max_RS;
    }

    // Fixed-seconds compatibility mode.  Explicitly selecting it disables
    // the period-scaled common tau/sigma EMA.
    void setAdaptationTimeConstants(float tau_sec) {
        if (std::isfinite(tau_sec) && tau_sec > 0.0f) {
            adapt_tau_sec_ = tau_sec;
            adapt_tau_sea_periods_ = 0.0f;
        }
    }

    // Common tau/sigma EMA horizon as a fraction/multiple of measured T_sea=T_z/2.
    // This changes only averaging memory; tau_coeff_, sigma_coeff_, the sigma
    // variance K-period horizon, and the already swept r_S horizon are untouched.
    void setAdaptationSeaPeriods(float periods) {
        if (std::isfinite(periods) && periods > 0.0f) {
            adapt_tau_sea_periods_ = periods;
        }
    }

    // Frequency EMA inside SeaStateAutoTuner, likewise expressed in T_sea=T_z/2.
    void setTunerFreqSmoothingSeaPeriods(float periods) {
        if (std::isfinite(periods) && periods > 0.0f) {
            tuner_freq_ema_sea_periods_ = periods;
            tuner_.setFrequencySmoothingSeaPeriods(periods);
        }
    }

    // Fixed-seconds frequency-EMA compatibility mode for controlled ablations.
    void setTunerFreqSmoothingTimeConstant(float tau_sec) {
        if (std::isfinite(tau_sec) && tau_sec > 0.0f) {
            tuner_freq_ema_sea_periods_ = 0.0f;
            tuner_.setTauFreq(tau_sec);
        }
    }

    // Smoothing-horizon multiplier for the r_S channel.  The EMA time
    // constant is mult * tau_target, so the horizon follows the sea state
    // instead of being pinned to one second count.  r_S ~ sigma_aw * tau^3
    // amplifies a tau error by the third power, which is what sets the
    // multiplier apart from the tau/sigma one; see
    // docs/ou-ema-adaptation-tuning.md.
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

    void setFreezeAccBiasUntilLive(bool en) { freeze_acc_bias_until_live_ = en; }

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

    void setWarmupRaccStd(float r) { if (std::isfinite(r) && r > 0.0f) Racc_warmup_std_ = r; }

    // For SeaStateFusionFilter_OU_III to restore Racc automatically
    void setNominalRaccStd(const Eigen::Vector3f& r) { Racc_nominal_ = r; }

    //  Exposed getters
    inline float getFreqHz()        const noexcept { return freq_hz_; }        // fast branch
    inline float getFreqSlowHz()    const noexcept { return freq_hz_slow_; }   // slow branch
    inline float getFreqRawHz()     const noexcept { return f_raw; }
    inline float getTauApplied()    const noexcept { return tune_.tau_applied; }
    inline float getSigmaApplied()  const noexcept { return tune_.sigma_applied; }
    inline float getRSApplied()     const noexcept { return tune_.RS_applied; }
    inline float getTauTarget()     const noexcept { return tau_target_;   }
    inline float getSigmaTarget()   const noexcept { return sigma_target_; }
    inline float getRSTarget()      const noexcept { return RS_target_;    }

    // Apparent period of the *acceleration* band, from the slow tracker
    // branch.  This is a reporting channel: the OU operating point uses
    // getWavePeriodSec(), the zero-crossing period of the elevation, which is
    // a different and much longer quantity.
    inline float getPeriodSec() const noexcept {
        return (freq_hz_slow_ > 1e-6f) ? 1.0f / freq_hz_slow_ : NAN;
    }

    // Variance measured after the period-scaled sigma band.
    inline float getAccelVariance() const noexcept { return tuner_.getAccelVariance(); }

    // Up-positive vertical acceleration from the private Mahony observer: the
    // signal the tracker, the wave-period estimator and the sigma channel all
    // run on.  Measurement-only -- it reads no filter state.
    inline float getAccelVertical() const noexcept {
        return vertical_accel_comp_.verticalAccelUpMs2();
    }

    inline float getHeaveAbs() const noexcept { if (!mekf_) return NAN; return std::fabs(mekf_->get_position().z()); }

    inline float getDisplacementScale(bool smoothed = true) const noexcept {
        const float tau = smoothed ? tune_.tau_applied : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;
        if (!std::isfinite(sigma) || !std::isfinite(tau)) return NAN;
        constexpr float C_HS  = 2.0f * std::sqrt(2.0f) / (std::numbers::pi_v<float> * std::numbers::pi_v<float>);
        return C_HS * sigma * tau * tau / 2.0f;
    }

    float getVerticalSpeedEnvelopeMps(bool smoothed = true) const noexcept {
        const float tau   = smoothed ? tune_.tau_applied   : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;
        if (!(tau > 1e-6f) || !std::isfinite(tau) || !std::isfinite(sigma)) return NAN;
        constexpr float K = std::sqrt(2.0f) / std::numbers::pi_v<float>;
        const float v_env = K * sigma * tau;
        return std::isfinite(v_env) ? v_env : NAN;
    }

    // Zero-crossing wave period [s] from the independent accelerometer-only
    // estimator; NaN until it settles.
    inline float getWavePeriodSec() const noexcept { return wave_period_.getPeriodSec(); }
    inline bool wavePeriodReady() const noexcept { return wave_period_.isReady(); }

    // Select which vertical acceleration drives the wave-period estimator.
    // Complementary (default) levels with the private Mahony observer and is
    // measurement-only, so the tuner is outside the estimator's loop.  Leveled
    // is the older behaviour, which levels with the main filter's attitude and
    // closes that loop.  See the call site in updateTime for what it costs.
    void setWavePeriodInput(WavePeriodInputSource source) {
        wave_period_input_ = source;
    }
    WavePeriodInputSource wavePeriodInput() const noexcept {
        return wave_period_input_;
    }

    // Gains of the private Mahony observer that levels the default input.
    // two_kp sets the accelerometer-to-gyro correction corner, which must stay
    // below the wave band; see VerticalAccelComplementary.h.
    void setWavePeriodComplementaryGains(float two_kp, float two_ki) {
        vertical_accel_comp_.setGains(two_kp, two_ki);
    }

    // Wave-band frequency used before WavePeriodEstimator has a value.
    void setTuneFreqPriorHz(float hz) {
        if (std::isfinite(hz) && hz > 0.0f) tune_freq_prior_hz_ = hz;
    }
    float tuneFreqPriorHz() const noexcept { return tune_freq_prior_hz_; }

    // Bounds on the wave-band tuning frequency.  Distinct from setFreqBounds(),
    // which bounds the acceleration-band tracker that drives wave direction.
    void setTuneFreqBounds(float min_hz, float max_hz) {
        if (!(std::isfinite(min_hz) && std::isfinite(max_hz))) return;
        if (!(min_hz > 0.0f && max_hz > min_hz)) return;
        min_tune_freq_hz_ = min_hz;
        max_tune_freq_hz_ = max_hz;
    }

    // sigma_a averaging horizon, in periods of the tuning frequency, and its
    // absolute clamps in seconds.  The horizon moved with the operating point
    // when tuning moved to the wave band -- the same K is now 5-17 s instead of
    // about 4 -- so it is a tuning surface; docs/ou-sigma-horizon.md measures it.
    void setSigmaVarianceKPeriods(float k) { tuner_.setKPeriods(k); }
    float getSigmaVarianceKPeriods() const noexcept { return tuner_.getKPeriods(); }
    void setSigmaVarianceHorizonBounds(float min_s, float max_s) {
        tuner_.setVarianceHorizonBounds(min_s, max_s);
    }
    // Horizon currently in force [s]; the variance is a two-stage EWMA at it.
    float getSigmaVarianceHorizonSec() const noexcept {
        return tuner_.getVarianceHorizonSec();
    }

    inline WaveDirection getDirSignState() const noexcept { return dir_sign_state_; }

    // Propagation-plane angle relative to boat +X, modulo 180 degrees.
    inline float getWaveAxisDeg() const noexcept { return dir_filter_.getAxisDegrees(); }
    inline float getWaveDirectionDeg() const noexcept { return getWaveAxisDeg(); }

    // Fully directed apparent propagation angles observed by the moving boat.
    // These are encounter/apparent directions unless vessel-motion correction
    // is applied externally (see wave_dir/WaveEncounter.h).
    inline float getApparentWaveDirectionToDeg() const noexcept {
        return dir_sign_.getDirectedAngleDegrees();
    }
    inline float getApparentWaveDirectionFromDeg() const noexcept {
        return dir_sign_.getWaveFromAngleDegrees();
    }
    inline float getDirSenseCoherence() const noexcept {
        return dir_sign_.getCoherence();
    }

    inline auto& mekf() noexcept { return *mekf_; }
    inline const auto& mekf() const noexcept { return *mekf_; }

    inline KalmanWaveDirection& dir() noexcept { return dir_filter_; }
    inline const KalmanWaveDirection& dir() const noexcept { return dir_filter_; }

    inline WaveDirectionDetector<float>& dir_sign() noexcept { return dir_sign_; }
    inline const WaveDirectionDetector<float>& dir_sign() const noexcept { return dir_sign_; }

private:

    // Simple first-order low-pass filter for vertical accel → tracker input
    using FreqInputLPF = seastate::common::FreqInputLPF;
    using StillnessAdapter = seastate::common::StillnessAdapter;

    float band_noise_floor_sigma_() const noexcept {
        if (!sigma_wave_band_.isReady()) {
            return acc_noise_floor_sigma_;
        }
        const float gain = sigma_wave_band_.whiteNoiseVarianceGain();
        if (!(std::isfinite(gain) && gain >= 0.0f)) return acc_noise_floor_sigma_;
        return acc_noise_floor_sigma_ * std::sqrt(gain);
    }

    // Apply the online tuner output only at the next IMU-sample boundary.
    // adapt_mekf() may consume y_k and smooth its candidate during step k,
    // but this function runs before y_{k+1} reaches the MEKF. Therefore the
    // active schedule at step k+1 is measurable with respect to data through k.
    void apply_pending_online_tune_() {
        if (!online_tune_apply_pending_ || !mekf_) return;
        apply_ou_tune_(false);
        if (startup_stage_ == StartupStage::Live && enable_linear_block_) {
            apply_RS_tune_();
        }
        online_tune_apply_pending_ = false;
    }

    void apply_pseudo_update_cadence_() {
        if (!mekf_ || !tau_scaled_pseudo_cadence_) return;
        const float tau = tune_.tau_applied;
        if (!(std::isfinite(tau) && tau > 0.0f)) return;
        const float requested = pseudo_update_tau_ratio_ * tau;
        const float period = std::min(
            std::max(requested, pseudo_update_period_min_s_),
            pseudo_update_period_max_s_);
        mekf_->set_pseudo_update_period_s(period);
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
        apply_pseudo_update_cadence_();

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
        if (!periodic_aw_cov_sync_ || !mekf_) return;
        if (startup_stage_ != StartupStage::Live) return;
        if (time_ - last_aw_cov_sync_sec_ <= adapt_every_secs_) return;
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
    // information rate is proportional to 1/(r_S^2 T_S). Preserve the historical
    // 15 ms information rate by normalizing the filter-input standard deviation:
    //     r_S,filter = r_S,base * sqrt(T_0/T_S).
    // The base tuner value remains clamped to [min_R_S_, max_R_S_]. Do not clamp
    // again after this normalization: the smallest-sea operating point is already
    // on the 0.4 base floor and must move below 0.4 when T_S > 15 ms to preserve
    // r_S^2 T_S. With unclamped T_S proportional to tau this turns the existing
    // base C_R*sqrt(R_a)*tau^3 schedule into an effective tau^(5/2) schedule
    // at the filter input.
    float pseudo_update_information_rate_scale_() const noexcept {
        // The Riccati laws already contain the realized T_S, so renormalizing
        // them again would double-count the cadence.
        if (rs_law_ != RSAdaptationLaw::Cubic) return 1.0f;
        if (!tau_scaled_pseudo_cadence_ || !mekf_) return 1.0f;
        const float period = mekf_->get_pseudo_update_period_s();
        if (!(std::isfinite(period) && period > 0.0f)) return 1.0f;
        return std::sqrt(pseudo_update_fixed_period_s_ / period);
    }

    // Pseudo-update period the cadence scheduler would select for a given tau.
    // Used by the Riccati laws so that the target r_S and the target cadence
    // refer to the same operating point, including the safety clamps.
    float pseudo_update_period_for_(float tau) const noexcept {
        if (!tau_scaled_pseudo_cadence_) return pseudo_update_fixed_period_s_;
        if (!(std::isfinite(tau) && tau > 0.0f)) return pseudo_update_fixed_period_s_;
        return std::min(std::max(pseudo_update_tau_ratio_ * tau,
                                 pseudo_update_period_min_s_),
                        pseudo_update_period_max_s_);
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
        // p = 0 is the pure strong-observation law, which is the deployed
        // Cubic schedule up to the constant gain sigma_ref; p = 1 restores the
        // sigma_aw multiplier the deployed base no longer carries.  At the
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
        const float rs_xy = RSb * s * R_S_xy_factor_;
        mekf_->set_RS_noise(Eigen::Vector3f(
            rs_xy,
            rs_xy,
            RSb * s
        ));
    }

    void update_tuner(float dt, float a_vertical_measurement, float freq_hz_for_tuner) {
        // Use the previous smoothed tuner frequency for the current sample's
        // band corners.  This keeps the band motion smooth and one-sample
        // predictable while remaining measurement-only.  Until the frequency
        // EMA is ready, fall back to the current external frequency estimate.
        float f_for_sigma_band = tuner_.isFreqReady()
            ? tuner_.getFrequencyHz()
            : freq_hz_for_tuner;
        const float f_tune_floor = min_tune_freq_hz_;
        const float f_tune_ceil = max_tune_freq_hz_;
        if (!std::isfinite(f_for_sigma_band) || f_for_sigma_band < f_tune_floor) {
            f_for_sigma_band = f_tune_floor;
        }
        f_for_sigma_band = std::min(f_for_sigma_band, f_tune_ceil);

        const float a_for_variance =
            sigma_wave_band_.step(a_vertical_measurement, dt, f_for_sigma_band);

        tuner_.update(dt, a_for_variance, freq_hz_for_tuner);

        // Startup stage logic
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
                  if (startup_init_policy_ == StartupInitPolicy::MahonyProxy) {
                      // The operating point is trusted, but the attitude is
                      // not this filter's to decide.  Park here and let the
                      // bootstrap call goLive() once it has tilt and north.
                      startup_stage_   = StartupStage::TunerReady;
                      startup_stage_t_ = 0.0f;
                  } else {
                      enterLive_();
                  }
              }
              break;

           case StartupStage::TunerReady:
           case StartupStage::Live:
               break;
        }

        // The tuning frequency is a wave-band quantity and is bounded by the
        // wave band, not by the tracker's bounds: a developed sea has
        // T_z = 8.6 s, i.e. 0.12 Hz, well under the 0.2 Hz the tracker is
        // bounded to, and setFreqBounds() moves the demodulator carrier rather
        // than the OU operating point.
        float f_tune = tuner_.getFrequencyHz();
        if (!std::isfinite(f_tune) || f_tune < f_tune_floor) {
            f_tune = f_tune_floor;
        }
        if (f_tune > f_tune_ceil) {
            f_tune = f_tune_ceil;
        }

        const float band_noise_sigma = band_noise_floor_sigma_();
        const float var_noise = band_noise_sigma * band_noise_sigma;
        float var_total = var_noise;
        if (tuner_.isVarReady()) {
            var_total = std::max(0.0f, tuner_.getAccelVariance());
        }
        float var_wave = var_total - var_noise;
        if (var_wave < 0.0f) var_wave = 0.0f;

        if (freq_stillness_.isStill()) {
            const float still_t = std::max(0.0f, freq_stillness_.getStillTime());
            constexpr float STILL_VAR_DECAY_SEC = 1.0f;
            float atten = std::exp(-still_t / STILL_VAR_DECAY_SEC);
            atten = std::min(std::max(atten, 0.0f), 1.0f);
            var_wave *= atten;
        }

        var_wave = std::max(var_wave, 1e-6f);
        float sigma_wave = std::sqrt(var_wave);
        float tau_raw = tau_coeff_ * 0.5f / f_tune;

        if (enable_clamp_) {
            tau_target_   = std::min(std::max(tau_raw,  min_tau_s_), max_tau_s_);
            sigma_target_ = std::min(sigma_wave * sigma_coeff_,      max_sigma_a_);
        } else {
            tau_target_   = tau_raw;
            sigma_target_ = sigma_wave;
        }
        if (!tuner_.isVarReady()) {
            sigma_target_ = std::max(sigma_target_, std::max(0.05f, band_noise_sigma));
        }

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
        const float sea_time_sec = 0.5f / f_tune;  // T_sea = T_z/2
        adapt_mekf(dt, tau_target_, sigma_target_, RS_target_, sea_time_sec);
    }

    void adapt_mekf(float dt, float tau_t, float sigma_t, float RS_t,
                    float sea_time_sec) {
        float adapt_sec = adapt_tau_sec_;
        if (adapt_tau_sea_periods_ > 0.0f &&
            std::isfinite(sea_time_sec) && sea_time_sec > 0.0f) {
            const float safe_sea_time =
                seastate::tuner::limits::clampDynamicEmaTimeScaleSec(sea_time_sec);
            adapt_sec = seastate::tuner::limits::clampDynamicEmaHorizonSec(
                adapt_tau_sea_periods_ * safe_sea_time, dt);
        }
        const float alpha = 1.0f - std::exp(-dt / adapt_sec);

        const float RS_sec = seastate::common::adaptiveSmoothingHorizonSec(
            adapt_RS_mult_, tau_t, RS_t, tune_.RS_applied, adapt_RS_slew_log_, dt);
        const float alpha_RS = 1.0f - std::exp(-dt / RS_sec);

        // Every valid physical sample contributes to every EMA.  The deployed
        // activation cadence is deliberately kept separate from that sampling:
        // it throttles parameter commits, not physical measurements or EMA input.
        tune_.tau_applied   += alpha    * (tau_t   - tune_.tau_applied);
        tune_.sigma_applied += alpha    * (sigma_t - tune_.sigma_applied);
        tune_.RS_applied    += alpha_RS * (RS_t    - tune_.RS_applied);

        if (time_ - last_adapt_time_sec_ > adapt_every_secs_) {
            // y_k may change the smoothed candidate here, but the MEKF keeps
            // the schedule that was active before y_k arrived. Commit this
            // candidate at the beginning of updateTime(k+1).
            online_tune_apply_pending_ = true;
            last_adapt_time_sec_ = time_;
        }
    }

    // Vertical acceleration the wave-period estimator is driven by.  The
    // leveled ablation falls back to the complementary observer while heading
    // is not yet resolved, so the estimator is never fed a body-frame residual.
    float wave_period_input_ms2_(
        const wave_direction::HeadingFrameAcceleration<float>& leveled) const
    {
        const float a_comp = vertical_accel_comp_.verticalAccelUpMs2();
        switch (wave_period_input_) {
            case WavePeriodInputSource::Leveled:
                return leveled.heading_valid ? leveled.up_ms2 : a_comp;
            case WavePeriodInputSource::Complementary:
            default:
                return a_comp;
        }
    }

    // The frequency the whole adaptation path runs on: the sigma band's
    // corners, the sigma_a averaging horizon, and tau.  It is a wave-band
    // quantity and nothing else; the acceleration-band tracker never reaches
    // it, at any instant of the run.
    //
    // The estimator's own value is taken as soon as it exists rather than at
    // isReady(): the readiness gate wants a settled *statistic*, and it does
    // not clear until 60-85 s into a run, whereas a value that has survived the
    // integrator settling transient is already a far better wave-band estimate
    // than a constant.  Before that the fixed wave-band prior stands in.  Both
    // are exogenous, so the schedule is a pure function of the measurements at
    // every instant of the run rather than only after the gate clears.
    float tuner_frequency_hz_() const {
        const float wave_hz = wave_period_.getFrequencyHz();
        if (std::isfinite(wave_hz) && wave_hz > 0.0f) return wave_hz;
        return tune_freq_prior_hz_;
    }

    void resetTrackingState_() {
        tracker_policy_       = TrackingPolicy{};
        wave_period_          = WavePeriodEstimator{};
        vertical_accel_comp_.reset();
        sigma_wave_band_.reset();
        freq_input_lpf_       = FreqInputLPF{};
        freq_stillness_       = StillnessAdapter(g_std, min_freq_hz_, FREQ_GUESS);
        freq_input_lpf_.setCutoff(max_freq_hz_);
        freq_stillness_.setTargetFreqHz(min_freq_hz_);

        tuner_.reset();

        freq_fast_smoother_   = FirstOrderIIRSmoother<float>(FREQ_SMOOTHER_DT, 3.5f);
        freq_slow_smoother_   = FirstOrderIIRSmoother<float>(FREQ_SMOOTHER_DT, 10.0f);

        freq_hz_      = FREQ_GUESS;
        freq_hz_slow_ = FREQ_GUESS;
        f_raw         = FREQ_GUESS;

        dir_filter_ = KalmanWaveDirection(2.0f * static_cast<float>(M_PI) * FREQ_GUESS);
        dir_sign_.reset();
        dir_sign_state_ = UNCERTAIN;

        last_adapt_time_sec_ = time_;
        last_aw_cov_sync_sec_ = time_;
        online_tune_apply_pending_ = false;
    }

    void enterCold_() {
        startup_stage_   = StartupStage::Cold;
        startup_stage_t_ = 0.0f;

        if (!mekf_) return;
        mekf_->set_linear_block_enabled(false);

        accel_bias_locked_   = with_mag_;
        mag_updates_applied_ = 0;
        first_mag_update_time_  = NAN;

        if (freeze_acc_bias_until_live_) {
            mekf_->set_acc_bias_updates_enabled(false);
            mekf_->set_Racc_std(Eigen::Vector3f::Constant(Racc_warmup_std_));
            warmup_Racc_active_ = true;
        }
    }

    void enterLive_() {
        startup_stage_   = StartupStage::Live;
        startup_stage_t_ = 0.0f;

        if (!mekf_) return;
        apply_ou_tune_(true);
        mekf_->set_linear_block_enabled(enable_linear_block_);

        if (freeze_acc_bias_until_live_) {
            const bool allow_bias = !accel_bias_locked_ && !acc_bias_hold_;
            mekf_->set_acc_bias_updates_enabled(allow_bias);

            if (warmup_Racc_active_ &&
                Racc_nominal_.allFinite() &&
                Racc_nominal_.maxCoeff() > 0.0f)
            {
                mekf_->set_Racc_std(Racc_nominal_);
            }
            warmup_Racc_active_ = false;
        }

        if (enable_linear_block_) apply_RS_tune_();
    }

    StartupStage startup_stage_    = StartupStage::Cold;
    float        startup_stage_t_  = 0.0f;

    // Warmup behavior
    bool  freeze_acc_bias_until_live_ = true;
    float Racc_warmup_std_            = 0.6f;
    bool  warmup_Racc_active_         = false;
    Eigen::Vector3f Racc_nominal_     = Eigen::Vector3f::Constant(0.0f);

    bool accel_bias_locked_ = true;
    int  mag_updates_applied_ = 0;
    static constexpr int MAG_UPDATES_TO_UNLOCK = 250;
    int  mag_updates_to_unlock_ = MAG_UPDATES_TO_UNLOCK;
    bool acc_bias_hold_ = false;

    // StagedMekf here, MahonyProxy in SeaStateFusion_OU_III::Config, and the
    // asymmetry is deliberate.  Only something outside this class can perform
    // the handoff, so a filter driven directly through updateTime() with no
    // bootstrap above it would park at TunerReady forever if it defaulted to
    // the proxy policy.  This class's standalone contract -- feed it and it
    // becomes live on its own -- is therefore left exactly as it was, and the
    // policy is chosen by the wrapper that is actually able to honour it.
    StartupInitPolicy startup_init_policy_ = StartupInitPolicy::StagedMekf;

    bool   with_mag_;
    double time_;
    double last_adapt_time_sec_;

    float first_mag_update_time_ = NAN;

    float tilt_over_limit_sec_ = 0.0f;
    float tilt_reset_cooldown_sec_ = 0.0f;

    float freq_hz_       = FREQ_GUESS;
    float freq_hz_slow_  = FREQ_GUESS;
    float f_raw          = FREQ_GUESS;

    bool enable_clamp_ = true;
    bool enable_tuner_ = true;
    bool online_tune_apply_pending_ = false;

    // Per-channel freezes for the partial-adaptation ablation; see
    // setChannelFreeze().  Both false is the deployed fully adaptive filter.
    bool freeze_ou_channel_ = false;
    bool freeze_RS_channel_ = false;
    float frozen_tau_s_ = NAN;
    float frozen_sigma_a_ = NAN;
    float frozen_RS_ = NAN;

    bool enable_linear_block_ = true;

    // Covariance-inflation policy; see setPeriodicAwCovarianceSync.
    bool   periodic_aw_cov_sync_ = true;
    bool   congruent_aw_cov_sync_ = false;
    double last_aw_cov_sync_sec_ = 0.0;

    bool  tau_scaled_pseudo_cadence_ = true;
    float pseudo_update_tau_ratio_ = PSEUDO_UPDATE_TAU_RATIO_DEFAULT;
    float pseudo_update_period_min_s_ = PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT;
    float pseudo_update_period_max_s_ = PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT;
    float pseudo_update_fixed_period_s_ = PSEUDO_UPDATE_PERIOD_NOMINAL_S;

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

    WavePeriodInputSource wave_period_input_ = WavePeriodInputSource::Complementary;
    float tune_freq_prior_hz_     = TUNE_FREQ_PRIOR_HZ;
    float min_tune_freq_hz_       = MIN_TUNE_FREQ_HZ;
    float max_tune_freq_hz_       = MAX_TUNE_FREQ_HZ;
    float min_freq_hz_            = MIN_FREQ_HZ;
    float max_freq_hz_            = MAX_FREQ_HZ;
    float min_tau_s_              = MIN_TAU_S;
    float max_tau_s_              = MAX_TAU_S;
    float max_sigma_a_            = MAX_SIGMA_A;
    float min_R_S_                = MIN_R_S;
    float max_R_S_                = MAX_R_S;
    float adapt_tau_sec_          = ADAPT_TAU_SEC;
    float adapt_tau_sea_periods_      = ADAPT_TAU_SEA_PERIODS;
    float tuner_freq_ema_sea_periods_ = TUNER_FREQ_EMA_SEA_PERIODS;
    float adapt_RS_mult_          = ADAPT_RS_MULT;
    float adapt_RS_slew_log_      = ADAPT_RS_SLEW_LOG;
    float adapt_every_secs_       = ADAPT_EVERY_SECS;
    float online_tune_warmup_sec_ = ONLINE_TUNE_WARMUP_SEC;
    float mag_delay_sec_          = MAG_DELAY_SEC;

    // Horizontal integral-regularization scale relative to the vertical one.
    // 0.36 made the horizontal high-pass 2.8x stronger than the vertical one,
    // which was a small-sea optimum applied to every sea state: with the
    // operating point now tied to the wave band, every stationary record scores
    // better with the two equal, by 7 to 27 percent of 3D RMS in the two
    // largest seas.  Retained as a setter because the bound is a real one.
    float R_S_xy_factor_ = 1.0f;
    // Horizontal stationary acceleration scale relative to the vertical one.
    // 1.87 was carried from the acceleration-band operating point and never
    // re-measured against the records, which put it at 0.81 (x) and 0.55 (y)
    // of vertical per axis -- 0.99 combined, as deep-water theory requires.
    // Swept at the deployed isotropic r_S over the eight scored records and
    // five seeds (tools/ou_anisotropy_ablation.py), 3D displacement RMS has a
    // flat minimum exactly at 1: -1.28 percent pooled against 1.87, -7.7
    // percent on JONSWAP Hs = 8.5 m and -2.5 percent on Hs = 4.0 m with every
    // seed agreeing, PM-Stokes flat to 0.15 percent, vertical unchanged at
    // +0.16 percent pooled, and no record worse than +0.14 percent.  0.9 and
    // 1.1 both score -1.19 percent, so the minimum is genuinely here and not a
    // fitted edge.
    //
    // Note that 1 is also the axis-consistent value: with sigma_aw isotropic,
    // an isotropic r_S is what the similarity law sigma_S ~ sigma_aw tau^3
    // asks for.  Closing that gap the other way -- keeping 1.87 and setting
    // R_S_xy_factor to match -- was measured and is 13.7 percent worse; see
    // docs/ou-iii-anisotropy-consistency.md.
    float S_factor_      = 1.0f;

    TrackingPolicy                  tracker_policy_{};
    FirstOrderIIRSmoother<float>    freq_fast_smoother_{FREQ_SMOOTHER_DT, 3.5f};
    FirstOrderIIRSmoother<float>    freq_slow_smoother_{FREQ_SMOOTHER_DT, 10.0f};
    SeaStateAutoTuner               tuner_;
    WavePeriodEstimator             wave_period_;
    // One private Mahony observer, serving both the vertical channel and the
    // startup attitude.
    //
    // These were briefly two instances, because the two jobs disagree about
    // the integral term.  The vertical channel had always run at two_ki = 0
    // and accepted the ~2b/two_kp static tilt a gyro bias leaves, on the
    // grounds that its two high-pass stages reject anything static.  Nothing
    // high-passes an *attitude seed*: the same error is a standing roll and
    // pitch bias for the whole run, and seeding from a zero-integral observer
    // measurably costs 0.17 deg of roll RMS on jonswap H1.5.
    //
    // Turning the integral term on for both settles that disagreement in the
    // only way that leaves one observer.  It does change the vertical channel,
    // which is why it was checked rather than assumed: over the eight
    // reference records the wave-period, sigma and r_S channels are unmoved
    // and every scored metric is at least as good.  Estimating the bias is
    // also the better answer for the vertical channel on its own terms -- the
    // static tilt it used to tolerate was leaking gravity into the levelled
    // acceleration, it was simply being high-passed away afterwards.
    //
    // two_kp stays at 0.2 for the reason it always did: the correction corner
    // must sit an order of magnitude below the wave band, or the observer
    // levels itself against the orbital specific force instead of gravity.
    VerticalAccelComplementary      vertical_accel_comp_{
        STARTUP_PROXY_TWO_KP_DEFAULT,
        STARTUP_PROXY_TWO_KI_DEFAULT};

    AdaptiveWaveBandPass            sigma_wave_band_{
        SIGMA_BAND_LOW_RATIO_DEFAULT,
        SIGMA_BAND_HIGH_RATIO_DEFAULT,
        SIGMA_BAND_MIN_HZ_DEFAULT,
        SIGMA_BAND_MAX_HZ_DEFAULT};
    TuneState                       tune_;

    float tau_target_   = NAN;
    float sigma_target_ = NAN;
    float RS_target_    = NAN;

    float acc_noise_floor_sigma_ = ACC_NOISE_FLOOR_SIGMA_DEFAULT;

    // The three adaptation coefficients of Eq. (adapt-three-layer-summary):
    //     tau      = tau_coeff * T_z / 2
    //     sigma_aw = sigma_coeff * sigma_a,B
    //     r_S,base = R_S_coeff * sqrt(R_a) * tau^3
    // tau_coeff = 1 is the documented intent, tau equal to half the
    // zero-crossing period.
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
    // sigma_coeff no longer enters r_S at all, and the filter is correspondingly
    // insensitive to it: the whole axis from 0.27 to 1.8 spans 3 % of the
    // vertical endpoint at the optimal C_R.  It stays at the 0.9 that shipped
    // with the old law.  The deterministic sweep prefers ~0.32, but that
    // preference does not survive the multi-seed protocol and costs 0.027 deg
    // of pitch; the band-matching prediction c_sigma = F_OU^(-1/2) ~ 1.80 is
    // significantly worse on the vertical endpoint, if only by 0.037 %Hs.  See
    // docs/ou-iii-rs-amplitude-retune.md and
    // tools/ou_rs_amplitude_retune_sweep.py.
    float R_S_coeff_    = R_S_COEFF_ANALYTICAL_REFERENCE;
    float tau_coeff_    = 1.0f;
    float sigma_coeff_  = 0.9f;

    std::unique_ptr<Kalman3D_Wave_OU_III<float>>  mekf_;
    KalmanWaveDirection                    dir_filter_{2.0f * static_cast<float>(M_PI) * FREQ_GUESS};

    FreqInputLPF        freq_input_lpf_;
    StillnessAdapter    freq_stillness_;

    WaveDirectionDetector<float> dir_sign_{0.002f, 0.005f};
    WaveDirection                dir_sign_state_ = UNCERTAIN;
};

template<TrackerType trackerT>
class SeaStateFusion_OU_III {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using StartupInitPolicy =
        typename SeaStateFusionFilter_OU_III<trackerT>::StartupInitPolicy;

    struct Config {
        bool with_mag = true;

        bool enable_linear_block = true;
        bool require_mag_lock_for_linear_block = false;

        // Who solves the startup attitude; see StartupInitPolicy.
        //
        // MahonyProxy is the default.  The measurement-only front end runs
        // from the first sample, the private Mahony observer supplies the tilt
        // that gates the magnetometer and frames the world-reference average,
        // and the MEKF is seeded with the finished solution and starts live.
        // StagedMekf restores the previous staged behaviour, in which the MEKF
        // is fed from the first sample and those same reads come back out of
        // it while it is still warming.
        StartupInitPolicy startup_init_policy = StartupInitPolicy::MahonyProxy;

        // Earliest and latest the proxy bootstrap may hand over.
        //
        // The normal exit is by quality: proxy tilt holding gravity agreement,
        // magnetic north gauged, and the tuner ready.  The floor keeps a
        // record whose first seconds happen to look calm from handing over on
        // a tilt the observer has barely integrated, and the ceiling
        // guarantees the filter always starts -- a platform that never satisfies
        // the gate still gets a live filter, on the best attitude available,
        // rather than sitting in bootstrap forever.
        float proxy_startup_min_sec     = 8.0f;
        float proxy_startup_timeout_sec = 150.0f;

        // Magnetic acquisition runs in two stages, because the two things it
        // has to deliver want opposite schedules.
        //
        // A usable heading is wanted within seconds of power-on.  A *good*
        // reference wants the startup observer to have settled first: its
        // correction corner sits below the wave band by design, which is what
        // stops it chasing orbital acceleration, and the same low corner means
        // it needs tens of seconds to converge from its accelerometer seed.
        // Measured against truth on the eight reference records, mean tilt
        // error over a 23 s averaging window starting at 7 s runs 0.33 to 2.76
        // deg; starting at 40 s it is 0.13 to 0.85 deg.
        //
        // Waiting for the good one before reporting anything would put first
        // heading around 105 s, which is not a usable device.  So the first
        // stage locks a provisional reference as soon as the gravity gate
        // allows -- heading and a live filter in roughly 20 s, as before --
        // and the second stage re-learns the reference once the estimate has
        // actually converged, in the *MEKF's* own tilt frame rather than the
        // observer's, since by then the MEKF has the linear block, the
        // magnetometer and the bias states and its tilt is the better of the
        // two.  The correction lands long before the scored window opens.
        //
        // proxy_mag_settle_sec holds the provisional stage off; 0 means "as
        // soon as the gravity gate is happy" and is the default, because the
        // refinement is what carries the accuracy now.
        float proxy_mag_settle_sec = 0.0f;

        bool  mag_refine_enabled    = true;
        float mag_refine_start_sec  = 90.0f;
        float mag_refine_window_sec = 30.0f;

        // Covariance seeded at handoff.  Tilt has been integrated through the
        // wave band by an observer whose correction corner is below it, so it
        // is worth about the accel-only default; yaw is either gauged by the
        // magnetometer or entirely unknown, and those two cases are an order
        // of magnitude apart, which is the whole reason the seed is split.
        float proxy_handoff_tilt_sigma_rad      = 0.035f;  // ~2 deg
        float proxy_handoff_yaw_sigma_rad       = 0.087f;  // ~5 deg, north gauged
        float proxy_handoff_yaw_sigma_free_rad  = 1.5708f; // ~90 deg, no lock

        float mag_delay_sec          = MAG_DELAY_SEC;
        float online_tune_warmup_sec = 10.0f;

        bool  freeze_acc_bias_until_live = true;
        float Racc_warmup_std = 0.5f;

        // Magnetometer updates that must land after the filter goes live
        // before the accelerometer-bias gate opens.
        //
        // Accelerometer bias and a tilt error are only weakly separable in
        // waves -- a roll error tips gravity into body Y and reads as a Y bias
        // -- so opening this gate while the attitude is still settling lets the
        // bias absorb the error and hold it.  Under the proxy policy the filter
        // reaches live far earlier than it used to, which moves this gate
        // earlier in absolute terms unless it is set to account for that.
        int acc_bias_unlock_mag_updates = 250;

        Eigen::Vector3f sigma_a = Eigen::Vector3f(0.2f, 0.2f, 0.2f);
        Eigen::Vector3f sigma_g = Eigen::Vector3f(0.01f, 0.01f, 0.01f);
        Eigen::Vector3f sigma_m = Eigen::Vector3f(0.3f, 0.3f, 0.3f);

        // The remaining MEKF variances the Kalman3D_Wave_OU_III constructor
        // takes.  They were reachable only through initialize_ext(), which
        // this wrapper never called, so every deployment ran on the header
        // defaults; docs/ou-iii-qmekf-variances.md is the sweep that gauged
        // them.  The values here reproduce those defaults exactly.
        //
        //   Pq0      initial attitude-error variance, rad^2.  Under the
        //            default MahonyProxy startup the handoff overwrites the
        //            attitude block, so this reaches the filter only under
        //            StagedMekf.
        //   Pb0      initial gyro-bias variance, (rad/s)^2.
        //   b0       gyro-bias random-walk variance density, (rad/s)^2/s.
        //   R_S_noise  initial integral pseudo-measurement variance, (m*s)^2.
        //            The tuner overwrites it at Live; it sets the pre-Live
        //            value and seeds the commanded-parameter filter.
        float Pq0       = 5e-4f;
        float Pb0       = 1e-6f;
        float b0        = 1e-11f;
        float R_S_noise = 1.5f;
        float gravity_magnitude = g_std;

        // Period-scaled sigma-band shape.  These are dimensionless wave-band
        // ratios except for the absolute safety clamps.  Keeping the ratios
        // fixed is what gives the JONSWAP sigma channel its similarity law.
        float sigma_band_low_ratio  = SIGMA_BAND_LOW_RATIO_DEFAULT;
        float sigma_band_high_ratio = SIGMA_BAND_HIGH_RATIO_DEFAULT;
        float sigma_band_min_hz     = SIGMA_BAND_MIN_HZ_DEFAULT;
        float sigma_band_max_hz     = SIGMA_BAND_MAX_HZ_DEFAULT;

        // Mag-start gate: gravity-direction agreement using current tilt.
        float mag_gravity_align_max_sin   = 0.075f; // sin(deg)
        float mag_gravity_align_hold_sec  = 2.0f;
        float mag_gravity_align_lpf_tau   = 1.0f;
        float mag_tilt_fallback_sec       = 30.0f;
        float mag_extreme_gyro_dps        = 30.0f; // veto only truly violent motion
        float mag_init_min_mag_norm       = 1e-3f;

        // Real-device mag acquisition.
        //
        // The important rule:
        //
        //   MagAutoTuner must never receive the MEKF's yaw, which is arbitrary
        //   and unobservable before mag lock.  It receives the MEKF quaternion
        //   with yaw divided out instead: that tilt is invariant under
        //   q_bw -> Rz(psi) q_bw, so no heading can leak through it, and it is
        //   a far better level frame in waves than one rebuilt from accel.
        //
        // The reference is an average of the field in that tilt frame, so
        // whatever tilt error survives the window survives in the reference.
        // In waves the error is periodic, so what the window has to buy is
        // whole wave periods, not samples: 128 samples is 5.1 s at a 25 Hz mag
        // ODR, short enough to lock in the phase it started on rather than
        // cancel it.  15 s covers a couple of periods across the band these
        // filters work in and captures most of what a much longer window would,
        // at a startup cost of 15 s rather than 40 s.  Held in seconds so it
        // does not silently shorten at a higher ODR.
        int   mag_min_samples              = 128;
        float mag_min_window_sec           = 15.0f;
        float mag_max_window_sec           = 0.0f;
        float mag_sample_dt_sec            = 1.0f / 200.0f;

        // Keep off in waves. Accel/gyro weighting can phase-select wave motion.
        // Body-frame hard-iron offset, learned during startup alongside the
        // reference and then subtracted from every magnetometer sample.
        //
        // Off by default.  The MEKF has no mag-bias state, so an offset left in
        // the stream is heading error one-for-one against the horizontal field;
        // but the offset is only weakly separable from the reference at a fixed
        // heading, and a wrong one subtracted everywhere is worse than none.
        // Turn it on where the platform changes heading during startup.
        bool  mag_estimate_hard_iron = false;

        // Continuous hard-iron estimation, and the reference that goes with it.
        //
        // The startup estimate above is a single window, and a single window is
        // where the offset is least identifiable: the excitation is whatever
        // tilt the hull happened to take in fifteen seconds.  This one never
        // closes its accumulation, so the question stops being "can the offset
        // be read out of these fifteen seconds" and becomes "keep watching, and
        // correct when the data finally say something".
        //
        // Two things make that safe to leave running.  The estimator is
        // exogenous -- gravity-referenced tilt from the private Mahony observer
        // and the raw magnetometer, never a filter state, so no loop is closed
        // through the MEKF and the ISS argument is untouched.  And the applied
        // offset and the magnetic reference move together, out of the same
        // statistics, so the filter is never subtracting one offset while
        // steering to a reference that belongs to another.
        //
        // Nothing here starts until the two-stage startup acquisition has
        // finished, so first heading, handoff and refinement are exactly as
        // they were.
        bool  mag_continuous_hard_iron        = true;
        float mag_hi_memory_sec               = 600.0f;
        // Absolute floor only; the relative term below carries the
        // calibration.  See ContinuousMagHardIronEstimator::Config and
        // docs/continuous-mag-hard-iron.md for why it came down from 4e-3.
        float mag_hi_model_ridge              = 5.0e-4f;
        float mag_hi_model_ridge_relative     = 0.5f;
        float mag_hi_min_information          = 2.0f;
        float mag_hi_min_effective_weight     = 500.0f;
        float mag_hi_max_residual_rms_uT      = 3.0f;
        float mag_hi_max_bias_fraction        = 0.35f;

        // Fraction of the fitted offset the filter is willing to apply, and the
        // time constant it moves over.  The ridge already shrinks the fit for
        // what the model cannot see; this is the separate, blunter statement
        // that a calibration nobody has checked should not arrive as a step.
        float mag_hi_apply_fraction           = 1.0f;
        float mag_hi_slew_tau_sec             = 45.0f;

        bool  mag_enable_quality_weighting = false;
        float mag_min_effective_weight     = 0.0f;
        float mag_acc_norm_rel_soft        = 0.22f;
        float mag_gyro_soft_dps            = 45.0f;

        // Bootstrap tilt observer for dynamic motion in waves.
        float bootstrap_tilt_obs_acc_tau_sec  = 2.15f; // accel correction time constant
        float bootstrap_gravity_slow_tau_sec  = 6.0f; // slow gravity reference LPF
        float bootstrap_gravity_align_max_sin = 0.070f; // sin(deg)
        float bootstrap_gravity_hold_sec      = 2.0f;
        float bootstrap_gravity_min_sec       = 6.87f;
        float bootstrap_gravity_timeout_sec   = 15.0f;
        float bootstrap_gravity_norm_frac     = 0.22f; // downweight accel when |a| departs from g

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
        mag_gravity_aligned_branch_ = false;
        mag_init_eligible_t0_ = NAN;
        last_mag_sample_t_ = NAN;

        mag_ref_set_ = false;

        last_mag_tilt_frame_yaw_rad_ = NAN;
        last_mag_startup_yaw_correction_rad_ = NAN;

        MagAutoTuner::Config mag_cfg;
        mag_cfg.mag_norm_min = cfg_.mag_init_min_mag_norm;
        mag_cfg.min_samples  = cfg_.mag_min_samples;
        mag_cfg.min_window_sec = cfg_.mag_min_window_sec;
        mag_cfg.max_window_sec = cfg_.mag_max_window_sec;
        mag_cfg.sample_dt_sec  = cfg_.mag_sample_dt_sec;

        mag_cfg.gravity_ref = g_std;
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

        mag_hi_startup_body_uT_.setZero();
        mag_hi_applied_body_uT_.setZero();
        mag_hi_anchor_bias_body_uT_.setZero();
        mag_hi_anchor_world_ref_uT_.setZero();
        mag_hi_anchored_ = false;
        last_hi_sample_t_ = NAN;
        last_hi_apply_t_  = NAN;
        mag_hard_iron_body_uT_.setZero();

        mag_world_ref_uT_.setZero();
        mag_world_ref_valid_ = false;

        resetTiltInit_();

        last_acc_body_ned_.setZero();
        last_gyro_body_ned_.setZero();
        have_last_imu_ = false;

        pending_yaw_abs_rad_ = NAN;

        mag_refine_started_  = false;
        mag_refine_done_     = false;
        mag_refine_time_sec_ = NAN;
        mag_north_lock_time_sec_ = NAN;
        live_time_sec_           = NAN;

        impl_.setWithMag(cfg_.with_mag);
        impl_.setStartupInitPolicy(cfg_.startup_init_policy);
        impl_.setFreezeAccBiasUntilLive(cfg_.freeze_acc_bias_until_live);
        impl_.setMagUpdatesToUnlockAccBias(cfg_.acc_bias_unlock_mag_updates);

        // The provisional reference is deliberately cheap and early, so the
        // accelerometer bias must not be allowed to fit itself to it; see
        // SeaStateFusionFilter_OU_III::setAccBiasHold().
        impl_.setAccBiasHold(usingProxyInit_() && cfg_.with_mag &&
                             cfg_.mag_refine_enabled);
        impl_.setWarmupRaccStd(cfg_.Racc_warmup_std);
        impl_.setMagDelaySec(0.0f); // outer wrapper owns mag delay
        impl_.setOnlineTuneWarmupSec(cfg_.online_tune_warmup_sec);
        impl_.setSigmaWaveBandRatios(cfg_.sigma_band_low_ratio,
                                     cfg_.sigma_band_high_ratio);
        impl_.setSigmaWaveBandLimitsHz(cfg_.sigma_band_min_hz,
                                       cfg_.sigma_band_max_hz);

        impl_.initialize_ext(cfg_.sigma_a, cfg_.sigma_g, cfg_.sigma_m,
                             cfg_.Pq0, cfg_.Pb0, cfg_.b0, cfg_.R_S_noise,
                             cfg_.gravity_magnitude);
        last_impl_startup_stage_ = impl_.getStartupStage();

        syncLinearBlockGate_();

        impl_.setNominalRaccStd(cfg_.sigma_a);

        displacement_up_m_.setZero();
        displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};

        if (cfg_.enable_displacement_detrend) {
            if (cfg_.use_custom_displacement_detrend_cfg) {
                displacement_detrender_.setConfig(cfg_.displacement_detrend_cfg);
            } else {
                displacement_detrender_.setConfig(
                    seastate::common::defaultDisplacementDetrenderConfig<AdaptiveWaveDetrender3D::Config>(FREQ_GUESS));
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

        const bool proxy_init = usingProxyInit_();

        // Under the proxy policy the gravity-lock bootstrap is the Mahony
        // observer inside the front end, which runs from the first sample, so
        // there is no phase in which the wrapper withholds IMU data.
        if (!proxy_init && stage_ == Stage::Uninitialized) {
            const bool tilt_ready = seastate::common::runStartupGravityInit(
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
                });

            if (tilt_ready) {
                stage_ = Stage::Warming;
            }
        }

        last_acc_body_ned_  = acc_body_ned;
        last_gyro_body_ned_ = gyro_body_ned;
        have_last_imu_      = true;

        if (proxy_init || stage_ != Stage::Uninitialized) {
            if (proxy_init && stage_ != Stage::Live) {
                // Bootstrap: front end only, MEKF held.
                impl_.updateFrontEnd(dt, gyro_body_ned, acc_body_ned);
            } else {
                impl_.updateTime(dt, gyro_body_ned, acc_body_ned, tempC);
            }

            const Eigen::Vector3f acc_gate_lp =
                gravity_gate_acc_lpf_.step(
                    acc_body_ned,
                    dt,
                    cfg_.mag_gravity_align_lpf_tau);

            // Whose tilt the magnetometer gate is judged against.  Under the
            // proxy policy this is the observer's before handoff, so the gate
            // measures the attitude that will actually frame the magnetic
            // reference rather than one the MEKF is still converging toward.
            const float align_sin =
                seastate::common::gravityAlignResidualSin(
                    attitudeReferenceQuat_(),
                    acc_gate_lp);

            // The sine residual is the same at an angle and at its supplement,
            // so it accepts an attitude flipped through 180 deg just as
            // readily as the right one.  The branch is the sign of the dot
            // product, and the gate has to carry it: this gate is what
            // certifies the tilt that frames the magnetic reference and that
            // is handed to the MEKF.
            const bool aligned_branch =
                seastate::common::gravityAlignedBranch(
                    attitudeReferenceQuat_(),
                    acc_gate_lp);

            mag_gravity_aligned_branch_ = aligned_branch;

            const float gyro_dps =
                gyro_body_ned.norm() * 57.295779513f;

            const bool extreme_motion =
                !std::isfinite(gyro_dps) ||
                (gyro_dps > cfg_.mag_extreme_gyro_dps);

            const bool gravity_good_now =
                std::isfinite(align_sin) &&
                (align_sin <= cfg_.mag_gravity_align_max_sin) &&
                aligned_branch &&
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

            const Eigen::Vector3f pos_ned_m = impl_.mekf().get_position();

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
                    (wave_hz >= displacement_detrender_.config().min_wave_freq_hz) &&
                    (wave_hz <= displacement_detrender_.config().max_wave_freq_hz);

                displacement_det_out_ =
                    displacement_detrender_.update(
                        displacement_up_m_,
                        dt,
                        wave_hz,
                        ext_freq_valid);
            } else {
                displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};
                displacement_det_out_.input = displacement_up_m_;
                displacement_det_out_.baseline_slow = Eigen::Vector3f::Zero();
                displacement_det_out_.wave_raw = displacement_up_m_;
                displacement_det_out_.wave_clean = displacement_up_m_;
            }
        }

        const auto cur_stage = impl_.getStartupStage();

        if (cur_stage != last_impl_startup_stage_) {
            if (cur_stage == SeaStateFusionFilter_OU_III<trackerT>::StartupStage::Cold) {
                mag_ref_set_ = false;
                mag_auto_tuner_.reset();

                gravity_gate_acc_lpf_.reset();
                mag_gravity_good_sec_ = 0.0f;
                mag_gravity_aligned_branch_ = false;
                mag_init_eligible_t0_ = NAN;
                last_mag_sample_t_ = NAN;

                last_mag_tilt_frame_yaw_rad_ = NAN;
                last_mag_startup_yaw_correction_rad_ = NAN;
                
                syncLinearBlockGate_();

                if (stage_ != Stage::Live) {
                    stage_ = Stage::Warming;

                    displacement_up_m_.setZero();
                    displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};

                    if (cfg_.enable_displacement_detrend) {
                        displacement_detrender_.reset(0.0f, 0.0f, 0.0f);
                    }
                }
            }

            last_impl_startup_stage_ = cur_stage;
        }

        if (proxy_init) {
            if (stage_ != Stage::Live) {
                maybeHandOffToMekf_();
            }
        } else if (stage_ == Stage::Warming && impl_.isAdaptiveLive()) {
            stage_ = Stage::Live;
            live_time_sec_ = t_;
        }

        // Re-apply gate every update.
        // Inner impl only enables the actual MEKF linear block when its own stage is Live.
        syncLinearBlockGate_();
    }

    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!begun_ || !cfg_.with_mag) return;
        // Under the proxy policy there is no withheld stage to wait out: the
        // observer has been levelling since the first sample, and learning
        // north before the MEKF starts is the entire point.  The quality gate
        // below still decides when accumulation may begin.
        if (!usingProxyInit_() && stage_ == Stage::Uninitialized) return;
        if (t_ < cfg_.mag_delay_sec) return;

        // Ahead of every gate below, and deliberately.  The continuous
        // estimator wants the whole magnetometer record, not the part the
        // startup machinery was willing to average, and it is reading a frame
        // the startup machinery does not own.
        accumulateContinuousHardIron_(mag_body_ned);

        // Hold the whole magnetometer path off until the startup observer has
        // settled.  This sits ahead of the eligibility clock deliberately, so
        // the tilt-fallback timer cannot start running and then wave the
        // accumulation through on an attitude that is still converging.
        //
        // What this window does *not* fix is the standing yaw error, which the
        // proxy policy leaves about 2 percent worse than the staged one
        // (1.795 -> 1.839 deg mean over the reference records).  That residual
        // is a property of the learned reference vector rather than of the
        // one-time gauge: seeding the handoff with a yaw variance from 5 deg
        // to 90 deg moves the scored yaw by under 1e-4 deg, so the filter is
        // not converging to the gauge, it is converging to the reference.
        // Improving it further means re-learning that reference in the live
        // MEKF's own tilt frame once it has converged, which is a runtime
        // re-acquisition rather than an initialization change and is
        // deliberately not attempted here.
        if (usingProxyInit_() && !mag_ref_set_ &&
            t_ < cfg_.proxy_mag_settle_sec) {
            return;
        }

        if (!std::isfinite(mag_init_eligible_t0_)) {
            mag_init_eligible_t0_ = t_;
        }

        const bool gravity_trusted =
            (mag_gravity_good_sec_ >= cfg_.mag_gravity_align_hold_sec);

        const bool fallback_ok =
            ((t_ - mag_init_eligible_t0_) >= cfg_.mag_tilt_fallback_sec);

        if (!mag_ref_set_) {
            if (!gravity_trusted && !fallback_ok) {
                return;
            }

            if (have_last_imu_) {
                const float dt_mag =
                    (std::isfinite(last_mag_sample_t_) && t_ > last_mag_sample_t_)
                        ? (t_ - last_mag_sample_t_)
                        : cfg_.mag_sample_dt_sec;

                last_mag_sample_t_ = t_;

                // Accumulate in a tilt frame with yaw removed.
                //
                // Stripping yaw makes the frame invariant to the estimator's
                // arbitrary startup heading, so this leaks no yaw into the
                // learned reference -- q_bw and Rz(psi) q_bw give the same
                // tilt.  A gravity-only frame rebuilt from low-passed accel
                // would be yaw-free too, but in waves that accel is gravity
                // plus a phase-lagged remnant of the orbital specific force,
                // so its tilt is wrong by a wave-correlated angle that the
                // averaging window is too short to cancel.
                //
                // Which estimator supplies that tilt is the substance of the
                // startup policy.  StagedMekf reads it back out of the warming
                // MEKF, whose accelerometer update is at that moment fighting
                // the orbital acceleration with its linear block switched off.
                // MahonyProxy reads the private observer instead, which is
                // gyro-propagated through the wave band and corrected below it.
                // The reference and the yaw gauge are locked exactly once, so
                // whichever tilt error survives this window is a standing bias
                // for the rest of the run.
                const Eigen::Quaternionf q_tilt_bw =
                    tiltOnlyQuatFromBoatQuat_(attitudeReferenceQuat_());

                if (mag_auto_tuner_.addSampleWithTiltQuatDt(
                        dt_mag,
                        q_tilt_bw,
                        last_acc_body_ned_,
                        last_gyro_body_ned_,
                        mag_body_ned))
                {
                    Eigen::Vector3f mag_world_ref_uT;

                    if (mag_auto_tuner_.getMagWorldRef(mag_world_ref_uT) &&
                        mag_world_ref_uT.allFinite() &&
                        mag_world_ref_uT.norm() > cfg_.mag_init_min_mag_norm)
                    {
                        // This reference was learned in a yaw-stripped tilt
                        // frame, so it carries no estimator heading.  It is a
                        // model parameter, so writing it before the MEKF has
                        // been handed the attitude is safe and leaves it ready
                        // to use the magnetometer from its first live sample.
                        setMagWorldRef_(mag_world_ref_uT);

                        const float mag_tilt_yaw_rad =
                            mag_auto_tuner_.getYawGaugeCorrectionRad();

                        if (std::isfinite(mag_tilt_yaw_rad)) {
                            // One-time yaw-gauge lock.
                            //
                            // mag_tilt_yaw_rad is the heading of averaged
                            // magnetic north in the accumulation frame, so
                            // driving the boat's absolute yaw to its negative
                            // puts the learned north on +X.
                            //
                            // Only yaw is written.  The MEKF's tilt is the
                            // gyro-propagated, accel-corrected estimate and is
                            // strictly better than any instantaneous
                            // accel-derived level frame, so overwriting the
                            // whole quaternion here would inject the wave tilt
                            // error the filter has already rejected.
                            const float yaw_abs_rad =
                                wrapPi_(-mag_tilt_yaw_rad);

                            // Before handoff there is no MEKF attitude to
                            // rewrite; the gauge is carried to the handoff
                            // instead and composed with the proxy tilt there,
                            // so the filter's very first attitude already has
                            // north in it.
                            if (usingProxyInit_() && stage_ != Stage::Live) {
                                pending_yaw_abs_rad_ = yaw_abs_rad;

                                last_mag_tilt_frame_yaw_rad_ =
                                    wrapPi_(mag_tilt_yaw_rad);

                                last_mag_startup_yaw_correction_rad_ =
                                    yaw_abs_rad;
                            } else {
                                Eigen::Quaternionf q_bw =
                                    impl_.mekf().quaternion_boat();
                                q_bw.normalize();

                                const Eigen::Quaternionf q_new =
                                    boatQuatWithAbsoluteYaw_(q_bw, yaw_abs_rad);

                                if (q_new.coeffs().allFinite()) {
                                    impl_.mekf().set_quaternion_boat(q_new);

                                    last_mag_tilt_frame_yaw_rad_ =
                                        wrapPi_(mag_tilt_yaw_rad);

                                    last_mag_startup_yaw_correction_rad_ =
                                        yaw_abs_rad;
                                }
                            }
                        }

                        Eigen::Vector3f hard_iron_uT;
                        mag_hi_startup_body_uT_ =
                            mag_auto_tuner_.getHardIronBodyUT(hard_iron_uT)
                                ? hard_iron_uT
                                : Eigen::Vector3f::Zero();
                        mag_hard_iron_body_uT_ = mag_hi_startup_body_uT_;

                        mag_ref_set_ = true;
                        mag_north_lock_time_sec_ = t_;
                        syncLinearBlockGate_();
                    }
                }
            }
        }

        maybeRefineMagReference_(mag_body_ned);
        maybeApplyContinuousHardIron_();

        // Magnetometer corrections go to the MEKF only once it owns the
        // attitude.  Before handoff its state is not the one being solved.
        if (mag_ref_set_ && (!usingProxyInit_() || stage_ == Stage::Live)) {
            impl_.updateMag(mag_body_ned - mag_hard_iron_body_uT_);
        }
    }

    bool hasMagNorthLock() const noexcept {
        return mag_ref_set_;
    }

    // True once the second-stage acquisition has replaced the provisional
    // reference; see maybeRefineMagReference_().
    bool hasRefinedMagReference() const noexcept { return mag_refine_done_; }
    float magRefineTimeSec() const noexcept { return mag_refine_time_sec_; }

    // Wall-clock marks for the startup sequence, for anyone measuring
    // time-to-first-fix rather than steady-state accuracy.
    float magNorthLockTimeSec() const noexcept { return mag_north_lock_time_sec_; }
    float liveTimeSec() const noexcept { return live_time_sec_; }

    // Body-frame hard-iron offset removed from the magnetometer stream.  Zero
    // unless Config::mag_estimate_hard_iron asked for it and the startup window
    // constrained it well enough to use.
    const Eigen::Vector3f& magHardIronBodyUT() const noexcept {
        return mag_hard_iron_body_uT_;
    }

    // Continuous estimator, for diagnostics and for tests that need to see why
    // it did or did not act.  estimate().valid is the gate; information() is
    // how much attitude excitation the memory window has actually collected.
    const ContinuousMagHardIronEstimator& magContinuousHardIron() const noexcept {
        return mag_hi_estimator_;
    }

    // The part of magHardIronBodyUT() the continuous estimator is responsible
    // for, as opposed to the startup solve.
    const Eigen::Vector3f& magContinuousHardIronAppliedUT() const noexcept {
        return mag_hi_applied_body_uT_;
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

    // Best available boat attitude, BODY -> WORLD (NED).
    //
    // Under the proxy policy the MEKF holds its initial quaternion until the
    // handoff, so reading it directly during the bootstrap would report a
    // level identity attitude rather than the platform's.  The bootstrap
    // observer's solution is what is actually known at that point, and it is
    // what this returns; after handoff both policies return the MEKF's.
    //
    // Heading is only meaningful once hasMagNorthLock() is true (or the build
    // has no magnetometer); before that the bootstrap yaw is arbitrary.  The
    // linear outputs -- displacement, velocity -- stay gated on isLive().
    Eigen::Quaternionf attitudeQuat() const {
        return attitudeReferenceQuat_();
    }

    Eigen::Vector3f eulerNauticalDeg() const {
        return impl_.getEulerNautical();
    }

    const Eigen::Vector3f& displacementUpMeters() const {
        return displacement_up_m_;
    }

    const AdaptiveWaveDetrender3D::Output& displacementDetrend() const {
        return displacement_det_out_;
    }

    SeaStateFusionFilter_OU_III<trackerT>& raw() {
        return impl_;
    }

    const SeaStateFusionFilter_OU_III<trackerT>& raw() const {
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

    float magTiltFrameYawDeg() const noexcept {
        return std::isfinite(last_mag_tilt_frame_yaw_rad_)
            ? last_mag_tilt_frame_yaw_rad_ * 57.29577951308232f
            : NAN;
    }

    float magStartupYawCorrectionDeg() const noexcept {
        return std::isfinite(last_mag_startup_yaw_correction_rad_)
            ? last_mag_startup_yaw_correction_rad_ * 57.29577951308232f
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

    using StartupTiltObserver = seastate::common::StartupTiltObserver;

    bool usingProxyInit_() const noexcept {
        return cfg_.startup_init_policy == StartupInitPolicy::MahonyProxy;
    }

    // The attitude the startup machinery judges and frames things against.
    //
    // Once the MEKF is live it is the answer under either policy -- it has the
    // magnetometer, the linear block and the bias states, and the proxy has
    // none of them.  Before that, under the proxy policy, the MEKF has nothing
    // to say and the observer does.
    Eigen::Quaternionf attitudeReferenceQuat_() const {
        if (usingProxyInit_() && stage_ != Stage::Live) {
            return impl_.startupProxyQuat();
        }
        return impl_.mekf().quaternion_boat();
    }

    // Hand the bootstrap attitude to the MEKF and start it live.
    //
    // The quality exit needs three things at once: a tilt that has held
    // agreement with gravity long enough to be trusted, a magnetic north gauge
    // (when a magnetometer is fitted), and an operating point the tuner
    // stands behind.  Waiting for all three is what lets the MEKF skip the
    // staged warmup entirely -- there is nothing left for it to converge.
    // Every write of the magnetometer's world reference goes through here, so
    // the wrapper always knows the vector the MEKF is steering to.  The MEKF
    // does not offer it back, and the continuous correction below moves the
    // reference by a delta rather than replacing it.
    void setMagWorldRef_(const Eigen::Vector3f& mag_world_ref_uT) {
        impl_.mekf().set_mag_world_ref(mag_world_ref_uT);
        mag_world_ref_uT_ = mag_world_ref_uT;
        mag_world_ref_valid_ = true;
    }

    // Feed the exogenous accumulation.  Raw magnetometer -- not the corrected
    // stream -- because the estimator is fitting the offset itself and must
    // not be shown data with its own answer already subtracted.
    void accumulateContinuousHardIron_(const Eigen::Vector3f& mag_body_ned) {
        if (!cfg_.mag_continuous_hard_iron) return;
        if (!usingProxyInit_()) return;
        if (!impl_.startupProxyInitialized()) return;

        const float dt_mag =
            (std::isfinite(last_hi_sample_t_) && t_ > last_hi_sample_t_)
                ? (t_ - last_hi_sample_t_)
                : cfg_.mag_sample_dt_sec;

        last_hi_sample_t_ = t_;

        mag_hi_estimator_.update(dt_mag,
                                 impl_.startupProxyTiltQuat(),
                                 mag_body_ned);
    }

    // Move the applied offset toward the fit, and re-gauge the reference.
    //
    // The reference has to be rebuilt in MagAutoTuner's canonical form --
    // horizontal magnitude on +X, vertical below it -- and not merely shifted
    // by the same amount as the measurement.  A shift that tracks the offset
    // exactly is a no-op: subtracting b from every sample and subtracting the
    // matching mean(R) b from the reference leaves the innovation identical at
    // the attitude the filter already holds, so nothing moves and the standing
    // yaw error survives the correction that was meant to remove it.
    //
    // The standing error is a *gauge*: the startup acquisition put the world
    // frame's north along the average of the uncorrected field, which is
    // magnetic north rotated by whatever the offset contributes.  Leaving the
    // canonical reference in place while the offset comes out of the stream
    // asks the filter for the heading the corrected field implies, and the
    // magnetometer update walks the yaw there over its own time constant.  No
    // attitude state is written, so the correction remains a change of
    // measurement-model parameters.
    //
    // Only the horizontal magnitude and the vertical component move with the
    // offset, and only by the amount the offset changes them.  They are not
    // recomputed from the estimator's own window: that window is longer and
    // less selective than the one the startup acquisition gated, and simply
    // adopting its magnitude and dip costs a fifth of the roll accuracy on the
    // moderate seas while the offset correction itself costs none of it.
    void maybeApplyContinuousHardIron_() {
        if (!cfg_.mag_continuous_hard_iron) return;
        if (!usingProxyInit_()) return;
        if (!mag_ref_set_ || stage_ != Stage::Live) return;
        if (cfg_.mag_refine_enabled && !mag_refine_done_) return;
        if (!mag_world_ref_valid_) return;

        const auto& est = mag_hi_estimator_.estimate();
        if (!est.valid) return;

        if (!mag_hi_anchored_) {
            mag_hi_anchor_bias_body_uT_ = mag_hi_applied_body_uT_;
            mag_hi_anchor_world_ref_uT_ = mag_world_ref_uT_;
            mag_hi_anchored_ = true;
        }

        const Eigen::Vector3f target =
            cfg_.mag_hi_apply_fraction * est.bias_body_uT;
        if (!target.allFinite()) return;

        const float dt_apply =
            (std::isfinite(last_hi_apply_t_) && t_ > last_hi_apply_t_)
                ? (t_ - last_hi_apply_t_)
                : cfg_.mag_sample_dt_sec;

        last_hi_apply_t_ = t_;

        const float tau = cfg_.mag_hi_slew_tau_sec;
        const float alpha = (std::isfinite(tau) && tau > 1.0e-3f)
                                ? (1.0f - std::exp(-dt_apply / tau))
                                : 1.0f;

        const Eigen::Vector3f applied =
            mag_hi_applied_body_uT_ + alpha * (target - mag_hi_applied_body_uT_);
        if (!applied.allFinite()) return;

        // Both evaluated against the statistics as they stand now, so the
        // difference is the offset's doing and nothing else.
        Eigen::Vector3f level_new;
        Eigen::Vector3f level_anchor;
        if (!mag_hi_estimator_.levelReferenceForBias(applied, level_new) ||
            !mag_hi_estimator_.levelReferenceForBias(mag_hi_anchor_bias_body_uT_,
                                                     level_anchor)) {
            return;
        }

        const float h_new = level_new.head<2>().norm();
        const float h_anchor = level_anchor.head<2>().norm();
        if (!std::isfinite(h_new) || !std::isfinite(h_anchor)) return;

        const float h = mag_hi_anchor_world_ref_uT_.x() + (h_new - h_anchor);
        const float z = mag_hi_anchor_world_ref_uT_.z() +
                        (level_new.z() - level_anchor.z());
        if (!(h > cfg_.mag_init_min_mag_norm) || !std::isfinite(h) ||
            !std::isfinite(z)) {
            return;
        }

        const Eigen::Vector3f ref(h, 0.0f, z);
        if (!ref.allFinite()) return;

        mag_hi_applied_body_uT_ = applied;
        mag_hard_iron_body_uT_ = mag_hi_startup_body_uT_ + applied;
        setMagWorldRef_(ref);
    }

    // Second-stage magnetic acquisition.
    //
    // The provisional reference was averaged in a tilt frame the startup
    // observer had barely converged, and it is what the filter has been
    // steering to ever since.  This re-runs the same acquisition once the MEKF
    // is live and settled, in the MEKF's own yaw-stripped tilt frame -- which
    // by then carries the linear block, the magnetometer and the bias states,
    // and is the better frame of the two.
    //
    // Both the reference vector and the heading gauge are replaced.  The yaw
    // write is a step, and deliberately so: it is the coarse-to-fine alignment
    // correction, it happens once, and it lands well before the scored window
    // opens.
    void maybeRefineMagReference_(const Eigen::Vector3f& mag_body_ned) {
        // Second-stage acquisition belongs to the proxy policy.  StagedMekf is
        // the ablation against the previous behaviour and has to stay exactly
        // that, so it keeps its single locked reference.
        if (!usingProxyInit_()) return;
        if (!cfg_.mag_refine_enabled) return;
        if (mag_refine_done_) return;
        if (!mag_ref_set_) return;
        if (stage_ != Stage::Live) return;
        if (t_ < cfg_.mag_refine_start_sec) return;
        if (!have_last_imu_) return;

        if (!mag_refine_started_) {
            MagAutoTuner::Config refine_cfg = mag_auto_tuner_.config();
            refine_cfg.min_window_sec = cfg_.mag_refine_window_sec;
            refine_cfg.min_samples    = cfg_.mag_min_samples;
            mag_auto_tuner_.setConfig(refine_cfg);
            mag_auto_tuner_.reset();
            mag_refine_started_  = true;
            last_mag_sample_t_   = NAN;
        }

        const float dt_mag =
            (std::isfinite(last_mag_sample_t_) && t_ > last_mag_sample_t_)
                ? (t_ - last_mag_sample_t_)
                : cfg_.mag_sample_dt_sec;

        last_mag_sample_t_ = t_;

        // The observer's tilt, not the MEKF's, for the same reason the first
        // stage used it -- and here the reason has teeth.
        //
        // By now the MEKF has the linear block and the bias states, so its
        // tilt looks like the better frame.  It is not usable: the MEKF has
        // been steering to the provisional reference this pass exists to
        // replace, so its tilt carries that reference's error, and averaging
        // the field in it re-derives the error it was meant to remove.  Tried
        // that way the refinement is self-confirming -- reference and yaw come
        // back within 1e-3 deg of the provisional ones, and the standing roll
        // bias is untouched.
        //
        // The observer never saw the reference, so its tilt is independent of
        // it, and by refinement time it has long since converged.
        const Eigen::Quaternionf q_tilt_bw = impl_.startupProxyTiltQuat();

        // Feed the same corrected stream the MEKF sees, so a hard-iron offset
        // already removed is not re-learned into the new reference.
        const Eigen::Vector3f mag_corrected = mag_body_ned - mag_hard_iron_body_uT_;

        if (!mag_auto_tuner_.addSampleWithTiltQuatDt(
                dt_mag, q_tilt_bw, last_acc_body_ned_,
                last_gyro_body_ned_, mag_corrected)) {
            return;
        }

        Eigen::Vector3f mag_world_ref_uT;
        if (!mag_auto_tuner_.getMagWorldRef(mag_world_ref_uT) ||
            !mag_world_ref_uT.allFinite() ||
            !(mag_world_ref_uT.norm() > cfg_.mag_init_min_mag_norm)) {
            return;
        }

        setMagWorldRef_(mag_world_ref_uT);

        const float mag_tilt_yaw_rad = mag_auto_tuner_.getYawGaugeCorrectionRad();

        if (std::isfinite(mag_tilt_yaw_rad)) {
            const float yaw_abs_rad = wrapPi_(-mag_tilt_yaw_rad);

            Eigen::Quaternionf q_bw = impl_.mekf().quaternion_boat();
            q_bw.normalize();

            const Eigen::Quaternionf q_new =
                boatQuatWithAbsoluteYaw_(q_bw, yaw_abs_rad);

            if (q_new.coeffs().allFinite()) {
                impl_.mekf().set_quaternion_boat(q_new);
                last_mag_tilt_frame_yaw_rad_ = wrapPi_(mag_tilt_yaw_rad);
                last_mag_startup_yaw_correction_rad_ = yaw_abs_rad;
            }
        }

        mag_refine_done_    = true;
        mag_refine_time_sec_ = t_;

        // The reference the bias would have been fitting is now the good one.
        impl_.setAccBiasHold(false);
    }

    void maybeHandOffToMekf_() {
        if (!begun_) return;

        const bool proxy_ready = impl_.startupProxyInitialized();

        const bool tilt_trusted =
            mag_gravity_aligned_branch_ &&
            (mag_gravity_good_sec_ >= cfg_.mag_gravity_align_hold_sec);

        const bool north_ready = !cfg_.with_mag || mag_ref_set_;

        const bool ready_by_quality =
            proxy_ready &&
            (t_ >= cfg_.proxy_startup_min_sec) &&
            tilt_trusted &&
            north_ready &&
            impl_.isTunerReady();

        // The timeout still requires an attitude to hand over; without one
        // there is nothing to seed and waiting costs nothing.
        //
        // It is also held clear of the magnetometer acquisition it would
        // otherwise cut short.  A timeout that fires while the reference is
        // still averaging hands over with no yaw gauge at all, which is a far
        // worse start than simply waiting: the gauge is the one chance to put
        // the filter on north before it goes live.  So the floor is the settle
        // time plus room for the averaging window to close.
        const float mag_acquire_deadline =
            cfg_.with_mag
                ? cfg_.proxy_mag_settle_sec +
                      2.0f * std::max(cfg_.mag_min_window_sec, 1.0f) +
                      cfg_.mag_tilt_fallback_sec
                : 0.0f;

        const float timeout_sec =
            std::max(cfg_.proxy_startup_timeout_sec, mag_acquire_deadline);

        // The timeout bounds how long startup may take; it does not license a
        // handoff onto the antipodal branch.  Seeding the MEKF with an
        // attitude that disagrees with measured gravity by more than a right
        // angle is worse than waiting out the extra samples it takes the
        // observer to leave a set it is not attracted to in the first place.
        const bool ready_by_timeout =
            proxy_ready &&
            (t_ >= timeout_sec) &&
            mag_gravity_aligned_branch_;

        if (!ready_by_quality && !ready_by_timeout) return;

        handOffToMekf_();
    }

    void handOffToMekf_() {
        const bool have_yaw_gauge = std::isfinite(pending_yaw_abs_rad_);

        // boatQuatWithAbsoluteYaw_ strips the incoming heading before writing
        // the new one, so the observer's drifted yaw cannot survive this even
        // though the quaternion is passed in whole.
        const Eigen::Quaternionf q_proxy = impl_.startupProxyQuat();

        const Eigen::Quaternionf q_seed =
            have_yaw_gauge
                ? boatQuatWithAbsoluteYaw_(q_proxy, pending_yaw_abs_rad_)
                : q_proxy;

        if (!q_seed.coeffs().allFinite()) return;

        const float yaw_sigma = have_yaw_gauge
            ? cfg_.proxy_handoff_yaw_sigma_rad
            : cfg_.proxy_handoff_yaw_sigma_free_rad;

        // The accelerometer-bias gate is deliberately left closed here.  Going
        // live early does not make that bias any more observable in waves, so
        // it keeps waiting for its count of magnetometer updates exactly as it
        // did before; see setMagUpdatesToUnlockAccBias().
        impl_.goLive(q_seed,
                     cfg_.proxy_handoff_tilt_sigma_rad,
                     yaw_sigma,
                     /*allow_acc_bias=*/false);

        stage_ = Stage::Live;
        live_time_sec_ = t_;
        last_impl_startup_stage_ = impl_.getStartupStage();

        syncLinearBlockGate_();
    }

    void resetTiltInit_() {
        bootstrap_tilt_obs_.reset();
        bootstrap_gravity_slow_lpf_.reset();
        bootstrap_gravity_good_sec_ = 0.0f;
    }

    bool linearBlockAllowed_() const {
        if (!cfg_.enable_linear_block) {
            return false;
        }
    
        if (!cfg_.with_mag) {
            return true;
        }
    
        if (!cfg_.require_mag_lock_for_linear_block) {
            return true;
        }
    
        return mag_ref_set_;
    }
    
    void syncLinearBlockGate_() {
        impl_.enableLinearBlock(linearBlockAllowed_());
    }

    static float wrapPi_(float a)
    {
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

    static float yawFromBoatQuatRad_(const Eigen::Quaternionf& q_bw_in) {
        if (!q_bw_in.coeffs().allFinite()) return NAN;

        Eigen::Quaternionf q_bw = q_bw_in;
        const float qn = q_bw.norm();

        if (!(qn > 1.0e-6f) || !std::isfinite(qn)) {
            return NAN;
        }

        q_bw.normalize();

        const Eigen::Matrix3f R = q_bw.toRotationMatrix();

        const float c = R(0, 0);
        const float s = R(1, 0);

        if (!std::isfinite(c) || !std::isfinite(s)) {
            return NAN;
        }

        return std::atan2(s, c);
    }

    // Tilt part of a BODY->WORLD quaternion, with the heading divided out.
    // Invariant under q_bw -> Rz(psi) q_bw, which is what makes it usable as a
    // mag accumulation frame that cannot leak the MEKF's arbitrary yaw.
    static Eigen::Quaternionf yawRemovedBoatQuat_(
        const Eigen::Quaternionf& q_bw_in)
    {
        if (!q_bw_in.coeffs().allFinite()) {
            return Eigen::Quaternionf::Identity();
        }

        Eigen::Quaternionf q_bw = q_bw_in;
        const float qn = q_bw.norm();

        if (!(qn > 1.0e-6f) || !std::isfinite(qn)) {
            return Eigen::Quaternionf::Identity();
        }

        q_bw.normalize();

        const float yaw = yawFromBoatQuatRad_(q_bw);

        if (!std::isfinite(yaw)) {
            return Eigen::Quaternionf::Identity();
        }

        const Eigen::Quaternionf q_yaw_inv(
            Eigen::AngleAxisf(-yaw, Eigen::Vector3f::UnitZ()));

        Eigen::Quaternionf q_tilt = q_yaw_inv * q_bw;
        q_tilt.normalize();

        if (!q_tilt.coeffs().allFinite()) {
            return Eigen::Quaternionf::Identity();
        }

        return q_tilt;
    }

    // Rewrites heading only, keeping the estimated tilt untouched.
    static Eigen::Quaternionf boatQuatWithAbsoluteYaw_(
        const Eigen::Quaternionf& q_bw_in,
        float yaw_abs_rad)
    {
        if (!std::isfinite(yaw_abs_rad)) {
            return q_bw_in;
        }

        const Eigen::Quaternionf q_tilt =
            yawRemovedBoatQuat_(q_bw_in);

        const Eigen::Quaternionf q_yaw(
            Eigen::AngleAxisf(
                yaw_abs_rad,
                Eigen::Vector3f::UnitZ()));

        Eigen::Quaternionf q_out = q_yaw * q_tilt;
        q_out.normalize();

        if (!q_out.coeffs().allFinite()) {
            return q_bw_in;
        }

        return q_out;
    }

    static Eigen::Quaternionf tiltOnlyQuatFromBoatQuat_(
        const Eigen::Quaternionf& q_bw_in)
    {
        return yawRemovedBoatQuat_(q_bw_in);
    }

private:
    Config cfg_{};
    SeaStateFusionFilter_OU_III<trackerT> impl_{false};

    bool begun_ = false;

    Stage stage_ = Stage::Uninitialized;
    float t_ = 0.0f;

    typename SeaStateFusionFilter_OU_III<trackerT>::StartupStage last_impl_startup_stage_ =
        SeaStateFusionFilter_OU_III<trackerT>::StartupStage::Cold;

    Eigen::Vector3f last_acc_body_ned_  = Eigen::Vector3f::Zero();
    Eigen::Vector3f last_gyro_body_ned_ = Eigen::Vector3f::Zero();
    bool have_last_imu_ = false;

    bool mag_ref_set_ = false;
    Eigen::Vector3f mag_hard_iron_body_uT_ = Eigen::Vector3f::Zero();
    MagAutoTuner mag_auto_tuner_{};

    ContinuousMagHardIronEstimator mag_hi_estimator_{};
    Eigen::Vector3f mag_hi_startup_body_uT_     = Eigen::Vector3f::Zero();
    Eigen::Vector3f mag_hi_anchor_bias_body_uT_ = Eigen::Vector3f::Zero();
    Eigen::Vector3f mag_hi_anchor_world_ref_uT_ = Eigen::Vector3f::Zero();
    bool  mag_hi_anchored_ = false;
    Eigen::Vector3f mag_hi_applied_body_uT_     = Eigen::Vector3f::Zero();
    float last_hi_sample_t_ = NAN;
    float last_hi_apply_t_  = NAN;

    Eigen::Vector3f mag_world_ref_uT_ = Eigen::Vector3f::Zero();
    bool mag_world_ref_valid_ = false;

    float last_mag_sample_t_ = NAN;

    float last_mag_tilt_frame_yaw_rad_ = NAN;
    float last_mag_startup_yaw_correction_rad_ = NAN;

    // Yaw gauge acquired while the MEKF was still held, applied at handoff.
    float pending_yaw_abs_rad_ = NAN;

    float mag_north_lock_time_sec_ = NAN;
    float live_time_sec_           = NAN;

    bool  mag_refine_started_  = false;
    bool  mag_refine_done_     = false;
    float mag_refine_time_sec_ = NAN;

    AdaptiveWaveDetrender3D displacement_detrender_{};
    AdaptiveWaveDetrender3D::Output displacement_det_out_{};
    Eigen::Vector3f displacement_up_m_ = Eigen::Vector3f::Zero();

    Vec3LPF gravity_gate_acc_lpf_{};
    float   mag_gravity_good_sec_ = 0.0f;
    bool    mag_gravity_aligned_branch_ = false;
    float   mag_init_eligible_t0_ = NAN;

    StartupTiltObserver bootstrap_tilt_obs_{};
    Vec3LPF             bootstrap_gravity_slow_lpf_{};
    float               bootstrap_gravity_good_sec_ = 0.0f;
};
