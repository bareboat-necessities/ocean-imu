#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy

  Released under the MIT License

  SeaStateFusionFilter_OU_II

  Marine Inertial Navigational System (INS) Filter for IMU

  Combines multiple real-time estimators into a cohesive ocean-state tracker:

    • Quaternion-based attitude and linear motion estimation via Kalman3D_Wave_OU_II

    • Dominant frequency tracking using one of:
          – AranovskiyFreqTracker     (frequency estimator)
          – KalmANFFreqTracker        (adaptive notch / Kalman frequency tracker)
          - PLLFreqTracker            (PLL frequency tracker)
          – SchmittTrigger            (zero-cross event detector)

    • Dual-stage frequency smoothing:
          – Fast 1st-order IIR (≈ few s, ~90% step) for demodulation / direction
          – Slow 1st-order IIR (≈ longer s, ~90% step) for auto-tuning / moments

    • Online auto-tuning of Kalman filter parameters (τ, σₐ, R_p0, R_v0) through
      SeaStateAutoTuner, which estimates acceleration variance and applies the
      σₐ·τ² and σₐ·τ regularization laws to stabilize displacement drift
      correction.  The time scale those laws are built on comes from
      WavePeriodEstimator, not from the acceleration-band frequency tracker,
      at every instant of the run: before that estimator has a value a fixed
      wave-band prior stands in rather than the tracker, and the tracker has no
      path into the adaptation at all.  See "Where" below.  The variance
      channel first applies a period-scaled measurement-only wave band, and
      averages its variance over a horizon of K wave periods.

  Where
  – τ (tau):  OU process time constant = ½ · T_z, half the zero-crossing period
              of the surface elevation as estimated by WavePeriodEstimator.
              It is deliberately *not* half the dominant period of acceleration:
              an ocean acceleration spectrum is the elevation spectrum weighted
              by (2πf)⁴, so its apparent frequency sits well above the spectral
              peak and barely moves as the sea grows, which pins τ to roughly
              one value across every sea state.
  – σₐ:       Stationary acceleration scale from period-scaled wave-band RMS
  – R_p0:     Pseudo-measurement noise controlling p drift suppression
  – R_v0:     Pseudo-measurement noise controlling v drift suppression
  – R_p0_xy:  Anisotropic X/Y weight on the p pseudo-measurement

  Adaptive update: tau/sigma smoothing over 0.40*T_sea (fixed seconds retained for ablation)

  ------------------------------------------------------------------------
  ADAPTATION AND STARTUP POLICY, shared with SeaStateFusionFilter_OU_III
  ------------------------------------------------------------------------

  The tuner and the attitude front end sit outside the estimator and do not
  know which of the two filters they are driving, so what OU-III measured
  about them applies here unchanged.  The two filters differ only in the
  translational state structure -- OU-III carries the extra integral
  displacement state and regularizes it with a single r_S, this one
  regularizes p and v separately -- and none of the following depends on that.

  1. EXOGENEITY IS A TIMING PROPERTY, NOT JUST A SIGNAL CHOICE.  Feeding the
     tuner from the complementary observer keeps its *inputs* independent of
     the filter.  That is necessary and not sufficient: the schedule smoothed
     during step k is committed at the top of the next update(), before
     y_{k+1} reaches the MEKF, so the active schedule at step k+1 is
     measurable with respect to data through k.

  2. THE VARIANCE CHANNEL RUNS IN A JONSWAP-SIMILAR WAVE BAND.  σₐ is
     estimated from the same exogenous levelled acceleration the wave-period
     estimator uses, after a band-pass whose two corners are fixed multiples
     of the tuner's own wave frequency.  Away from the absolute safety clamps
     the transfer shape is therefore fixed in f/f_tune, which is the condition
     the σₐ similarity argument needs.  The bench noise floor is referred
     through that band's own time-varying coefficients before subtraction.

  3. THE PSEUDO-MEASUREMENT CADENCE IS SELF-SIMILAR IN τ, AND THE
     REGULARIZERS FOLLOW IT.  One pseudo update has covariance r²; at one
     update per T_S seconds the continuous-equivalent information rate goes as
     1/(r² T_S).  Scaling T_S with τ while holding r fixed would therefore
     change the regularization strength with sea state as a side effect.
     T_S = (0.015/1.1)·τ, clamped to [5, 250] ms, with both filter inputs
     renormalized by sqrt(T_0/T_S).  That renormalization is deliberately not
     re-clamped, so the small-sea end may go below the base floor once
     T_S > T_0.

  4. STARTUP: THE PROXY OWNS TILT AND MAGNETIC LEARNING.  See
     SeaStateFusion_OU_II::maybeHandOffToMekf_() and
     docs/ou-iii-startup-init.md.

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
#include "tuner/AccelVibrationGuard.h"
#include "tuner/MagAutoTuner.h"
#include "tuner/ContinuousMagHardIronEstimator.h"
#include "kalman_ou_ii/Kalman3D_Wave_OU_II.h"
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

// JONSWAP-similar acceleration-variance band, matching OU-III's.  Away from
// the absolute safety clamps, its transfer shape is fixed in f/f_tune, which
// is the condition the sigma_aw similarity argument needs.  The band is a
// property of the tuner rather than of the estimator behind it, so nothing in
// its derivation depends on which of the two filters consumes the result.
constexpr float SIGMA_BAND_LOW_RATIO_DEFAULT  = 0.5f;
constexpr float SIGMA_BAND_HIGH_RATIO_DEFAULT = 4.0f;
constexpr float SIGMA_BAND_MIN_HZ_DEFAULT     = 0.01f;
constexpr float SIGMA_BAND_MAX_HZ_DEFAULT     = 6.0f;

constexpr float MIN_TAU_S     = 0.02f;  // sec
// tau now scales with the zero-crossing wave period rather than with an
// acceleration-band frequency, so the ceiling has to admit a developed sea:
// T_z reaches 8.6 s at H_s = 8.5 m and a long swell goes further.  The old
// 3.0 s ceiling was reached at H_s = 8.5 m, which clipped the operating point
// exactly where the filter was losing.
constexpr float MAX_TAU_S     = 12.0f;  // sec
constexpr float MAX_SIGMA_A   = 6.0f;
// r_p0 ~ sigma_aw * tau^2 and r_v0 ~ sigma_aw * tau inherit that range.  The
// old 18 m ceiling on r_p0 was the binding constraint at H_s = 8.5 m: the
// calibrated point sat exactly on it.  At the wave-band operating point the
// largest deployed values on the reference records are r_p0 = 20 m and
// r_v0 = 8.6 m/s at H_s = 8.5 m, so these ceilings sit five to seven times
// above the working range and act as saturation safeguards, not as tuning.
//
// regularizer_floor-test asks the *deployed* law whether these still hold.
// Under PhysicalMSE the near-still stress case asks for r_p0 = 0.070 m against
// this 0.05 floor, where the Empirical schedule asked 0.058 m, so the deployed
// law clears the floor by more than the one it replaced and the clamps are
// unchanged.
constexpr float MIN_R_p0_std  = 0.05f;
constexpr float MAX_R_p0_std  = 150.0f;
constexpr float MIN_R_v0_std  = 0.01f;
constexpr float MAX_R_v0_std  = 40.0f;

// Legacy fixed-second horizon remains for explicit ablations.  The deployed
// common tau/sigma EMA follows measured sea time T_sea=T_z/2, matching OU-III.
constexpr float ADAPT_TAU_SEC                  = 1.8f;
constexpr float ADAPT_TAU_SEA_PERIODS          = 0.40f;
constexpr float ADAPT_EVERY_SECS               = 0.1f;
// Smoothing horizons of the two drift-correction channels, in units of
// tau_target.  Measured on the versioned records against synthesized sea-state
// transitions: the error during a transition falls monotonically as these
// shorten, the stationary worst-record vertical error rises monotonically, and
// 3.0 is where the two cross at an acceptable cost for both channels.  See
// docs/ou-ema-adaptation-tuning.md.
constexpr float ADAPT_R_p0_MULT            = 3.0f;   // dimensionless
constexpr float ADAPT_R_v0_MULT            = 3.0f;   // dimensionless
// Discrepancy, in natural-log units of the target-to-applied ratio, above
// which the smoothing horizon of either drift-correction channel shortens.
// Zero keeps the plain proportional horizon; see
// docs/ou-ema-adaptation-tuning.md.
constexpr float ADAPT_R_SLEW_LOG           = 0.0f;   // ln units
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
// two_ki estimates the gyro bias.  The vertical channel ran at zero because
// everything downstream of it is high-passed, but an attitude seed keeps
// whatever static tilt the bias leaves -- about 2b/two_kp, i.e. 0.71 deg at
// 0.05 deg/s -- so the observer that serves both has to estimate it.
constexpr float STARTUP_PROXY_TWO_KP_DEFAULT = 0.2f;
constexpr float STARTUP_PROXY_TWO_KI_DEFAULT = 0.02f;

// Front-end accelerometer vibration guard, armed by default.
//
// docs/engine-noise-degradation.md is the measurement.  Machinery vibration
// rectifies into a standing tilt error and from there into a displacement
// offset, and at a routine cruise condition that costs a factor of eight in
// pooled 3-D error; the guard holds every engine condition swept inside a
// factor of 1.6 of the engine-off baseline.
//
// The study was run on OU-III, but nothing in it is specific to the
// translational state structure, which is the only thing OU-II and OU-III
// disagree about.  The mechanism is the attitude loop: vibration reaches the
// accelerometer, the attitude correction wobbles at the vibration
// frequencies, and the nonlinearity of the measurement model rectifies that
// wobble into a static tilt.  Both filters level from the same private Mahony
// observer with the same gains and feed the same accelerometer measurement to
// a MEKF, so both inherit the defect and the remedy.
//
// Arming it costs nothing when there is no machinery, because the guard is
// gated on its own detector and returns its input unchanged below the lower
// rail.  Across the eight stationary records the clean detector reading is
// 0.00796 to 0.00805 m/s^2 against a 0.03 rail, and the replays are
// bit-identical to the unguarded ones.  That is why this is a default rather
// than an option: there is no quiet-water case to trade away.
//
// Two poles at 14 Hz sits in the gap between the wave band and the lowest
// crank order a small auxiliary diesel puts on the hull, and costs 22.7 ms of
// group delay at full engagement.  Set the cutoff to zero to remove the guard
// entirely and restore the unconditioned measurement path.
constexpr float ACC_VIBRATION_GUARD_HZ_DEFAULT = 14.0f;
constexpr int   ACC_VIBRATION_GUARD_POLES_DEFAULT = 2;

// Vibration-aware accelerometer measurement covariance, on by default.
//
// Conditioning removes the machinery the guard can reach; what survives its
// stopband still arrives as measurement error the MEKF does not know about.
// This raises the commanded accelerometer sigma by the guard's own gated
// excess, so the covariance and the measurement describe the same conditions.
//
// The gain is swept in docs/engine-noise-degradation.md, on OU-III, and 0.75
// sits at the displacement optimum there with margin below the cliff: past
// about 1.25 the accelerometer is de-weighted enough that the wave estimate
// starts leaning on the OU prior instead, and displacement turns back up.
// The trade is between the same two channels here, so the same point is
// carried over rather than re-derived.
//
// Zero disables it.  Like the guard, it is driven by a gated excess that is
// identically zero on a quiet installation, so it is bit-transparent there.
constexpr float ACC_VIBRATION_RACC_GAIN_DEFAULT = 0.75f;

// Frequency smoother dt (SeaStateFusionFilter_OU_II is designed for 200 Hz)
constexpr float FREQ_SMOOTHER_DT = 1.0f / 200.0f;

// Self-similar drift-regularizer pseudo-measurement cadence.
//
// The S=0 pseudo updates of OU-III and the p=0/v=0 pseudo updates here share
// one mechanism and one 15 ms default: a single periodic tick inside the MEKF
// fires them.  One pseudo update has covariance r^2, so at one update per T_S
// seconds the continuous-equivalent information rate goes as 1/(r^2 T_S).
// Scaling T_S with tau while holding r fixed therefore changes the
// regularization strength with sea state as a side effect, which is not what
// either schedule is supposed to say.  The deployed wrapper historically used
// T_S = 15 ms while its initial applied OU time constant is tau = 1.1 s, so
// this ratio preserves that exact operating point and scales the cadence with
// tau thereafter.
constexpr float PSEUDO_UPDATE_PERIOD_NOMINAL_S = 0.015f;
constexpr float PSEUDO_UPDATE_TAU_NOMINAL_S = 1.1f;
constexpr float PSEUDO_UPDATE_TAU_RATIO_DEFAULT =
    PSEUDO_UPDATE_PERIOD_NOMINAL_S / PSEUDO_UPDATE_TAU_NOMINAL_S;
// A pseudo update cannot occur more often than the nominal 200 Hz IMU schedule.
// The upper guard is inactive over the current tau <= 12 s operating envelope.
constexpr float PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT = FREQ_SMOOTHER_DT;
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
// mu_-2 = M_-2 / (sigma_a^2 tau^6).  Evaluating that moment on the eight
// reference spectra at their measured OU-II operating points and solving
// Eq. (wp-opt) record by record gives C_P ~ 0.112, so the default is an
// analytical prediction rather than a fitted number.  See
// tools/ou2_dual_mse_coefficients.py and docs/ou-ii-dual-mse-adaptation.md.
constexpr float R_PSEUDO_MSE_COEFF_DEFAULT = 0.1116f;

// Channel ratio C_P/C_V of the PhysicalMSE law.  This is the corollary the
// note identifies as the more tightly constrained half of the prediction,
// because q_eff and the cadence normalization both cancel out of it:
//     (r_p / r_v)^2 = (3/2) M_-2 / M_0,
// hence r_p/r_v = C_P/C_V * tau for a fixed normalized sea shape.  The eight
// finite-band reference spectra give 0.431..0.485 with a mean of 0.461; an
// ideal infinite-band Pierson-Moskowitz sea gives sqrt(3/(4 pi)) = 0.489 and an
// ideal gamma = 3.3 JONSWAP 0.465.  The Empirical law's c_p/c_v = 0.500 is the
// same relative law with a slightly looser position channel, which is the
// agreement Sec. (ratio) of the note reports.
//
// Applying the ratio rather than a second independent power is exact, not an
// approximation: sigma^(4/5) tau^(7/5) = sigma^(4/5) tau^(12/5) / tau.  It is
// also what keeps the whole schedule down to a single transcendental.
constexpr float R_PSEUDO_MSE_RATIO_DEFAULT = 0.4611f;

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
// OU-II objective is not, and using the bench figure here was a
// misidentification the calibration sweep caught.
//
// The filter already carries a measured estimate of the right quantity:
// ACC_NOISE_FLOOR_SIGMA_DEFAULT, the pre-band vertical acceleration noise
// floor the tuner subtracts as non-wave energy, at 0.12 m/s^2 -- eight times
// the bench 0.0148 m/s^2.  Referring the law to it moves the schedule by
// (0.12/0.0148)^(1/5) = 1.52, and the eight-record scale sweep puts the
// complete-MEKF vertical optimum within 3 % of exactly that; see
// docs/ou-ii-dual-mse-adaptation.md.  So the analytical C_P is deployed as
// derived, and this constant is what changed.
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

// Unified SeaState fusion filter
template<TrackerType trackerT>
class SeaStateFusionFilter_OU_II {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using TrackingPolicy = TrackerPolicy<trackerT>;

    enum class StartupStage {
        Cold,        // just booted or just had a big tilt reset
        TunerWarm,   // front end running, tuner collecting stats
        TunerReady,  // tuner trusted, MEKF still held by an external bootstrap
        Live         // MEKF owns the attitude; full adaptation & extras allowed
    };

    explicit SeaStateFusionFilter_OU_II(bool with_mag = true)
        : with_mag_(with_mag),
          time_(0.0),
          last_adapt_time_sec_(0.0),
          freq_hz_(FREQ_GUESS),
          freq_hz_slow_(FREQ_GUESS)
    {
        // Default cutoff ~max_freq_hz_ Hz: passes waves, kills 8–37 Hz engine band
        freq_input_lpf_.setCutoff(max_freq_hz_);
        setAccelVibrationGuard(ACC_VIBRATION_GUARD_HZ_DEFAULT,
                               ACC_VIBRATION_GUARD_POLES_DEFAULT);
        setAccelVibrationRaccGain(ACC_VIBRATION_RACC_GAIN_DEFAULT);
        freq_stillness_.setTargetFreqHz(min_freq_hz_);
        startup_stage_   = StartupStage::Cold;
        startup_stage_t_ = 0.0f;
    }

    StartupStage getStartupStage() const noexcept { return startup_stage_; }
    bool isAdaptiveLive() const noexcept { return startup_stage_ == StartupStage::Live; }

    // The operating point is trustworthy.  It is reached while the MEKF is
    // still held, and it is one of the conditions the caller waits on before
    // handing the attitude over.
    bool isTunerReady() const noexcept {
        return startup_stage_ == StartupStage::TunerReady ||
               startup_stage_ == StartupStage::Live;
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
    // goLive(), then call this.  Both wrappers in this file do exactly that.
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

    // Out-of-band accelerometer guard, ahead of the proxy and the MEKF.
    //
    // The cutoff sits in the gap between the wave band and the machinery band
    // and stops vibration reaching the attitude loop, where it rectifies into a
    // standing tilt error.  Armed by default; pass zero to remove it and
    // restore the unconditioned measurement path exactly.  The cost is group
    // delay, accelVibrationGuardDelaySec(), which appears in displacement as
    // amplitude * 2 pi f * delay -- so prefer the highest corner that removes
    // the machinery, not the lowest corner that fits above the waves.
    void setAccelVibrationGuard(float cutoff_hz, int poles = 2) {
        accel_guard_.setPoles(poles);
        accel_guard_.setCutoffHz(cutoff_hz);
    }

    // Engagement band of the guard's detector.  Exposed mainly so a study can
    // force the guard on over a quiet input and separate the delay it costs
    // from the vibration it removes.
    void setAccelVibrationEngagement(float lo_mps2, float hi_mps2,
                                     float slew_tau_sec) noexcept {
        accel_guard_.setEngagement(lo_mps2, hi_mps2, slew_tau_sec);
    }

    [[nodiscard]] float accelVibrationGuardDelaySec() const noexcept {
        return accel_guard_.groupDelaySec();
    }

    // Vibration-aware accelerometer measurement covariance.
    //
    // The guard removes the machinery it can, but what survives its stopband
    // still reaches the MEKF as measurement error the filter does not know
    // about.  This tells it: the commanded accelerometer standard deviation
    // becomes sqrt(sigma_base^2 + (gain * excess)^2), where excess is the
    // guard's detector reading above its engagement floor.
    //
    // Zero disables it and leaves the commanded covariance exactly as the
    // startup and stage logic set it.  Because the drive is the guard's own
    // gated excess, it is identically zero on a quiet installation, so an
    // enabled gain is still bit-transparent there.
    void setAccelVibrationRaccGain(float gain) noexcept {
        if (std::isfinite(gain) && gain >= 0.0f) racc_vibration_gain_ = gain;
    }

    [[nodiscard]] float accelVibrationRaccGain() const noexcept {
        return racc_vibration_gain_;
    }

    // Accelerometer sigma currently commanded to the MEKF, m/s^2 per axis.
    [[nodiscard]] Eigen::Vector3f accelVibrationRaccStd() const noexcept {
        return racc_effective_;
    }

    [[nodiscard]] float accelVibrationGuardCutoffHz() const noexcept {
        return accel_guard_.cutoffHz();
    }

    [[nodiscard]] int accelVibrationGuardPoles() const noexcept {
        return accel_guard_.poles();
    }

    // How far the guard is currently engaged, in [0, 1].  Zero means the
    // measurement path is the unconditioned one.
    [[nodiscard]] float accelVibrationGuardEngagement() const noexcept {
        return accel_guard_.engagement();
    }

    // Smoothed RMS of the out-of-band content the guard is removing, m/s^2.
    // Zero when the guard is disabled, so it reads as a health signal only
    // where it is actually measuring something.
    [[nodiscard]] float accelVibrationRms() const noexcept {
        return accel_guard_.removedRms();
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
    void updateCore_(float dt, const Eigen::Vector3f& gyro, const Eigen::Vector3f& acc,
                     float tempC, bool drive_mekf)
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
        const Eigen::Vector3f acc_in = accel_guard_.step(acc, dt);

        // Private Mahony observer for the wave-period estimator, the default
        // sigma channel and the startup attitude.  It is fed gyro and
        // accelerometer before the MEKF sees them, so that the levelling it
        // provides stays a pure function of the measurements.  Stepping it
        // unconditionally keeps its transient off the critical path when the
        // input source is switched at runtime.
        vertical_accel_comp_.update(dt, gyro, acc_in, g_std);

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
                    mekf_->initialize_from_acc_preserve_yaw(acc_in);
                } else {
                    // During startup stages, accel-only re-lock is acceptable.
                    mekf_->initialize_from_acc(acc_in);
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

        // Tuner gets vertical accel and the wave-band frequency.  That single
        // substitution fixed both halves of the old operating point: tau is no
        // longer derived from an acceleration-band frequency that barely moves
        // with the sea state, and the tuner's variance horizon (a few periods)
        // is no longer shorter than one wave period, which was biasing
        // sigma_aw low.  The acceleration-band tracker reaches none of it; it
        // stays where it is needed, as the wave-direction demodulator carrier.
        //
        // The variance channel sees the same exogenous levelled acceleration
        // the wave-period estimator does, but only after a period-scaled
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
        // ablation; tests/kalman_ou_iii/tuner_coupling-test.cpp asserts the
        // default is exogenous bit-for-bit and bounds the Leveled path's gain.
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
    // Magnetometer correction
    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!with_mag_ || !mekf_) return;
        if (time_ < mag_delay_sec_) return;

        mekf_->measurement_update_mag_only(mag_body_ned);
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

    void setAccNoiseFloorSigma(float s) {
        if (std::isfinite(s) && s > 0.0f) acc_noise_floor_sigma_ = s;
    }

    float getAccNoiseFloorSigma() const noexcept { return acc_noise_floor_sigma_; }

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

    // Configure LPF on the levelled vertical acceleration the tracker runs on
    void setFreqInputCutoffHz(float fc) { freq_input_lpf_.setCutoff(fc); }

    void enableClamp(bool flag = true) { enable_clamp_ = flag; }
    void enableTuner(bool flag = true) {
        enable_tuner_ = flag;
        if (!flag) online_tune_apply_pending_ = false;
    }

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

    // Self-similar drift-regularizer pseudo-measurement cadence
    // T_S = c_T * tau_applied.  Enabled by default; disabling restores the
    // historical fixed 15 ms cadence for direct old-versus-new ablation.
    // Whenever the cadence changes while Live, reapply r_p0 and r_v0 so their
    // per-update covariances stay information-rate matched.
    void setTauScaledPseudoUpdateCadence(bool flag) {
        tau_scaled_pseudo_cadence_ = flag;
        if (!mekf_) return;
        if (flag) apply_pseudo_update_cadence_();
        else mekf_->set_pseudo_update_period_s(pseudo_update_fixed_period_s_);
        if (startup_stage_ == StartupStage::Live) {
            apply_R_p0_tune_();
            apply_R_v0_tune_();
        }
    }
    bool tauScaledPseudoUpdateCadence() const noexcept { return tau_scaled_pseudo_cadence_; }
    void setPseudoUpdateTauRatio(float ratio) {
        if (!(std::isfinite(ratio) && ratio > 0.0f)) return;
        pseudo_update_tau_ratio_ = ratio;
        if (tau_scaled_pseudo_cadence_) {
            apply_pseudo_update_cadence_();
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
            apply_pseudo_update_cadence_();
            if (startup_stage_ == StartupStage::Live) {
                apply_R_p0_tune_();
                apply_R_v0_tune_();
            }
        }
    }
    float getPseudoUpdateTauRatio() const noexcept { return pseudo_update_tau_ratio_; }
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

    // Fixed-seconds compatibility mode. Selecting it explicitly disables
    // sea-period scaling for the common tau/sigma EMA.
    void setAdaptationTimeConstants(float tau_sec) {
        if (std::isfinite(tau_sec) && tau_sec > 0.0f) {
            adapt_tau_sec_ = tau_sec;
            adapt_tau_sea_periods_ = 0.0f;
        }
    }

    void setAdaptationSeaPeriods(float periods) {
        if (std::isfinite(periods) && periods > 0.0f) {
            adapt_tau_sea_periods_ = periods;
        }
    }

    // No-ops kept so existing ablation scripts still compile and run.  The
    // tuner's second frequency EMA was removed when WavePeriodEstimator took
    // ownership of the canonical log-period state, and there is nothing left
    // here to smooth.
    void setTunerFreqSmoothingSeaPeriods(float /*periods*/) {}
    void setTunerFreqSmoothingTimeConstant(float /*tau_sec*/) {}

    float getAdaptationSeaPeriods() const noexcept { return adapt_tau_sea_periods_; }
    float getTunerFreqSmoothingSeaPeriods() const noexcept {
        return tuner_.getFrequencySmoothingSeaPeriods();
    }

    // Smoothing-horizon multipliers for the two drift-correction channels.
    // The EMA time constant of each channel is mult * tau_target, so the
    // horizon follows the sea state instead of being pinned to one second
    // count.  r_p0 ~ sigma_aw * tau^2 and r_v0 ~ sigma_aw * tau amplify a tau
    // error by different powers, so the two channels need not share a
    // multiplier; see docs/ou-ema-adaptation-tuning.md.
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

    void setAdaptationUpdatePeriod(float every_sec) {
        if (std::isfinite(every_sec) && every_sec > 0.0f) adapt_every_secs_ = every_sec;
    }

    void setOnlineTuneWarmupSec(float warmup_sec) {
        if (std::isfinite(warmup_sec) && warmup_sec >= 0.0f) online_tune_warmup_sec_ = warmup_sec;
    }

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

    inline float getFreqHz()           const noexcept { return freq_hz_; }
    inline float getFreqSlowHz()       const noexcept { return freq_hz_slow_; }
    inline float getFreqRawHz()        const noexcept { return f_raw; }
    inline float getTauApplied()       const noexcept { return mekf_ ? mekf_->get_aw_time_constant() : NAN; }
    inline float getSigmaApplied()     const noexcept { return mekf_ ? mekf_->get_aw_stationary_std().z() : NAN; }
    inline float getR_p0_std_applied() const noexcept { return mekf_ ? mekf_->get_Rp0_noise_std().z() : NAN; }
    inline float getR_v0_std_applied() const noexcept { return mekf_ ? mekf_->get_Rv0_noise_std().z() : NAN; }
    inline float getTauTarget()        const noexcept { return tau_target_; }
    inline float getSigmaTarget()      const noexcept { return sigma_target_; }
    inline float getR_p0_std_target()  const noexcept { return R_p0_std_target_; }
    inline float getR_v0_std_target()  const noexcept { return R_v0_std_target_; }

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

    inline float getHeaveAbs() const noexcept {
        if (!mekf_) return NAN;
        return std::fabs(mekf_->get_position().z());
    }

    inline float getDisplacementScale(bool smoothed = true) const noexcept {
        const float tau   = smoothed ? tune_.tau_applied   : tau_target_;
        const float sigma = smoothed ? tune_.sigma_applied : sigma_target_;
        if (!std::isfinite(sigma) || !std::isfinite(tau)) return NAN;
        constexpr float C_HS = 2.0f * std::sqrt(2.0f) / (std::numbers::pi_v<float> * std::numbers::pi_v<float>);
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

    using FreqInputLPF = seastate::common::FreqInputLPF;
    using StillnessAdapter = seastate::common::StillnessAdapter;

    // Bench white-noise sigma referred to the current adaptive band.  The band
    // is time varying, so the gain has to come from the filter's own state
    // rather than from a closed-form transfer function.
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
        if (startup_stage_ == StartupStage::Live) {
            apply_R_p0_tune_();
            apply_R_v0_tune_();
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
        // Commit the pseudo-update cadence with the same applied tau so
        // T_S/tau stays constant apart from explicit safety clamps.
        apply_pseudo_update_cadence_();

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
        if (!periodic_aw_cov_sync_ || !mekf_) return;
        if (startup_stage_ != StartupStage::Live) return;
        if (time_ - last_aw_cov_sync_sec_ <= adapt_every_secs_) return;
        mekf_->synchronize_aw_covariance_to_stationary();
        last_aw_cov_sync_sec_ = time_;
    }

    // set_Rp0_noise_std()/set_Rv0_noise_std() accept standard deviations, so one
    // pseudo update has covariance r^2.  With updates every T_S seconds, the
    // continuous-equivalent information rate is proportional to 1/(r^2 T_S).
    // Preserve the historical 15 ms information rate by normalizing the
    // filter-input standard deviations:
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
        const float period = mekf_->get_pseudo_update_period_s();
        if (!(std::isfinite(period) && period > 0.0f)) return 1.0f;
        return std::sqrt(pseudo_update_fixed_period_s_ / period);
    }

    // Pseudo-update period the cadence scheduler would select for a given tau.
    // Used by the PhysicalMSE law so that the target pair and the target
    // cadence refer to the same operating point, including the safety clamps.
    float pseudo_update_period_for_(float tau) const noexcept {
        if (!tau_scaled_pseudo_cadence_) return pseudo_update_fixed_period_s_;
        if (!(std::isfinite(tau) && tau > 0.0f)) return pseudo_update_fixed_period_s_;
        return std::min(std::max(pseudo_update_tau_ratio_ * tau,
                                 pseudo_update_period_min_s_),
                        pseudo_update_period_max_s_);
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
                    // The operating point is trusted, but the attitude is not
                    // this filter's to decide.  Park here and let the bootstrap
                    // call goLive() once it has tilt and north.
                    startup_stage_   = StartupStage::TunerReady;
                    startup_stage_t_ = 0.0f;
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
        if (!std::isfinite(f_tune) || f_tune < f_tune_floor) f_tune = f_tune_floor;
        if (f_tune > f_tune_ceil) f_tune = f_tune_ceil;

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
            tau_target_   = std::min(std::max(tau_raw, min_tau_s_), max_tau_s_);
            sigma_target_ = std::min(sigma_wave * sigma_coeff_, max_sigma_a_);
        } else {
            tau_target_   = tau_raw;
            sigma_target_ = sigma_wave;
        }

        if (!tuner_.isVarReady()) {
            sigma_target_ = std::max(sigma_target_, std::max(0.05f, band_noise_sigma));
        }

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

        const float sea_time_sec = 0.5f / f_tune;  // T_sea = T_z/2
        adapt_mekf(dt, tau_target_, sigma_target_, R_p0_std_target_, R_v0_std_target_,
                   sea_time_sec);
    }

    void adapt_mekf(float dt, float tau_t, float sigma_t, float R_p0_t, float R_v0_t,
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
        tracker_policy_ = TrackingPolicy{};
        wave_period_    = WavePeriodEstimator{};
        vertical_accel_comp_.reset();
        sigma_wave_band_.reset();
        freq_input_lpf_ = FreqInputLPF{};
        freq_stillness_ = StillnessAdapter(g_std, min_freq_hz_, FREQ_GUESS);
        freq_input_lpf_.setCutoff(max_freq_hz_);
        freq_stillness_.setTargetFreqHz(min_freq_hz_);

        tuner_.reset();

        freq_fast_smoother_ = FirstOrderIIRSmoother<float>(FREQ_SMOOTHER_DT, 3.5f);
        freq_slow_smoother_ = FirstOrderIIRSmoother<float>(FREQ_SMOOTHER_DT, 10.0f);

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

    void apply_racc_vibration_inflation_() {
        if (!mekf_ || !(racc_vibration_gain_ > 0.0f)) return;

        const Eigen::Vector3f base = racc_base_std_();
        if (!(base.minCoeff() > 0.0f)) return;

        const float excess = accel_guard_.excessRms();
        if (!(excess > 0.0f)) {
            // Hand the base back once on the way down, then stay quiet, so a
            // dormant guard leaves the stage logic's covariance untouched.
            if (racc_inflated_) {
                mekf_->set_Racc_std(base);
                racc_effective_ = base;
                racc_inflated_ = false;
            }
            return;
        }

        const float added = racc_vibration_gain_ * excess;
        const Eigen::Vector3f effective =
            (base.array().square() + added * added).sqrt().matrix();
        mekf_->set_Racc_std(effective);
        racc_effective_ = effective;
        racc_inflated_ = true;
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

    StartupStage startup_stage_ = StartupStage::Cold;
    float        startup_stage_t_ = 0.0f;

    // Vibration-aware measurement covariance, armed in the constructor at
    // ACC_VIBRATION_RACC_GAIN_DEFAULT and inert until the guard sees machinery.
    float racc_vibration_gain_        = 0.0f;
    bool  racc_inflated_              = false;
    Eigen::Vector3f racc_effective_   = Eigen::Vector3f::Zero();

    Eigen::Vector3f Racc_nominal_std_ = Eigen::Vector3f::Constant(0.0f);

    bool accel_bias_locked_ = true;
    int  mag_updates_applied_ = 0;
    static constexpr int MAG_UPDATES_TO_UNLOCK = 250;
    int  mag_updates_to_unlock_ = MAG_UPDATES_TO_UNLOCK;
    bool acc_bias_hold_ = false;

    bool   with_mag_;
    double time_;
    double last_adapt_time_sec_;

    float first_mag_update_time_ = NAN;

    float tilt_over_limit_sec_ = 0.0f;
    float tilt_reset_cooldown_sec_ = 0.0f;

    float freq_hz_      = FREQ_GUESS;
    float freq_hz_slow_ = FREQ_GUESS;
    float f_raw         = FREQ_GUESS;


    bool enable_clamp_ = true;
    bool enable_tuner_ = true;
    bool online_tune_apply_pending_ = false;

    // Covariance-inflation policy; see setPeriodicAwCovarianceSync.
    bool   periodic_aw_cov_sync_ = true;
    double last_aw_cov_sync_sec_ = 0.0;

    bool  tau_scaled_pseudo_cadence_ = true;
    float pseudo_update_tau_ratio_ = PSEUDO_UPDATE_TAU_RATIO_DEFAULT;
    float pseudo_update_period_min_s_ = PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT;
    float pseudo_update_period_max_s_ = PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT;
    float pseudo_update_fixed_period_s_ = PSEUDO_UPDATE_PERIOD_NOMINAL_S;

    WavePeriodInputSource wave_period_input_ = WavePeriodInputSource::Complementary;
    float tune_freq_prior_hz_     = TUNE_FREQ_PRIOR_HZ;
    float min_tune_freq_hz_       = MIN_TUNE_FREQ_HZ;
    float max_tune_freq_hz_       = MAX_TUNE_FREQ_HZ;
    float min_freq_hz_            = MIN_FREQ_HZ;
    float max_freq_hz_            = MAX_FREQ_HZ;
    float min_tau_s_              = MIN_TAU_S;
    float max_tau_s_              = MAX_TAU_S;
    float max_sigma_a_            = MAX_SIGMA_A;
    float MIN_R_p0_std_           = MIN_R_p0_std;
    float MAX_R_p0_std_           = MAX_R_p0_std;
    float MIN_R_v0_std_           = MIN_R_v0_std;
    float MAX_R_v0_std_           = MAX_R_v0_std;
    float adapt_tau_sec_              = ADAPT_TAU_SEC;
    float adapt_tau_sea_periods_      = ADAPT_TAU_SEA_PERIODS;
    float adapt_R_p0_mult_            = ADAPT_R_p0_MULT;
    float adapt_R_v0_mult_        = ADAPT_R_v0_MULT;
    float adapt_R_slew_log_       = ADAPT_R_SLEW_LOG;
    float adapt_every_secs_       = ADAPT_EVERY_SECS;
    float online_tune_warmup_sec_ = ONLINE_TUNE_WARMUP_SEC;
    float mag_delay_sec_          = MAG_DELAY_SEC;

    // Per-axis horizontal p-regularization scale, against the vertical one.
    // These were a single scalar until the split; the history below is that
    // scalar's, and both axes start from the same value.
    //
    // 0.31 was fitted against the acceleration-band operating point, where it
    // made the horizontal high-pass 3.2x stronger than the vertical one.  That
    // was a small-sea optimum applied to every sea state; with tau tied to the
    // wave band, 1.0 took about 10 percent off the mean 3D RMS across the four
    // stationary records while leaving normalized vertical error flat to within
    // 0.03 percentage points, and 1.0 is what shipped.
    //
    // Nothing between 0.31 and 1.0 was scored at that time.  Swept there over
    // the eight scored records and three IMU seed triplets
    // (tools/ou_low_sea_error_study.py xy --family OU_II), 3D displacement RMS
    // has an interior minimum at 0.65 -- -3.56 percent against 1, with the same
    // sign in all 24 record x seed cells -- and 0.8 and 0.55 bracket it at
    // -2.95 and -2.80.  Vertical is unchanged at -0.01 percent.  That is the
    // same shape and nearly the same optimum OU-III's own regularizer has, on a
    // different pseudo-measurement, which is the strongest evidence either
    // family offers that the effect is structural rather than a fit to one
    // wrapper.
    //
    // The deployed 0.72 is OU-III's measured optimum, carried here so the three
    // families share one horizontal-anisotropy constant.  It is not this
    // family's own argmin -- 0.65 is -- but the basin is flat enough that the
    // difference does not matter: measured on the same 24 cells, 0.72 scores
    // -3.50 percent of 3D RMS against 1 where 0.65 scores -3.56, so the shared
    // constant gives up six hundredths of a point.
    // See docs/ou-horizontal-anisotropy-per-axis-split.md.
    float R_p0_x_factor_ = 0.72f;
    float R_p0_y_factor_ = 0.72f;
    float P_factor_       = 1.5f;

    TrackingPolicy               tracker_policy_{};
    FirstOrderIIRSmoother<float> freq_fast_smoother_{FREQ_SMOOTHER_DT, 3.5f};
    FirstOrderIIRSmoother<float> freq_slow_smoother_{FREQ_SMOOTHER_DT, 10.0f};
    SeaStateAutoTuner            tuner_;
    WavePeriodEstimator          wave_period_;
    // One private Mahony observer, serving both the vertical channel and the
    // startup attitude.
    //
    // The two jobs disagree about the integral term.  The vertical channel had
    // always run at two_ki = 0 and accepted the ~2b/two_kp static tilt a gyro
    // bias leaves, on the grounds that its two high-pass stages reject
    // anything static.  Nothing high-passes an *attitude seed*: the same error
    // is a standing roll and pitch bias for the whole run.  Turning the
    // integral term on for both settles that in the only way that leaves one
    // observer, and it is also the better answer for the vertical channel on
    // its own terms -- the static tilt it used to tolerate was leaking gravity
    // into the levelled acceleration, it was simply being high-passed away
    // afterwards.
    //
    // two_kp stays at 0.2 for the reason it always did: the correction corner
    // must sit an order of magnitude below the wave band, or the observer
    // levels itself against the orbital specific force instead of gravity.
    VerticalAccelComplementary   vertical_accel_comp_{
        STARTUP_PROXY_TWO_KP_DEFAULT,
        STARTUP_PROXY_TWO_KI_DEFAULT};

    // Armed in the constructor at ACC_VIBRATION_GUARD_HZ_DEFAULT, and dormant
    // until its own detector sees machinery, so an unconditioned replay is
    // bit-identical to a guarded one.
    seastate::tuner::AccelVibrationGuard accel_guard_{};

    AdaptiveWaveBandPass         sigma_wave_band_{
        SIGMA_BAND_LOW_RATIO_DEFAULT,
        SIGMA_BAND_HIGH_RATIO_DEFAULT,
        SIGMA_BAND_MIN_HZ_DEFAULT,
        SIGMA_BAND_MAX_HZ_DEFAULT};
    TuneState                    tune_;

    float tau_target_      = NAN;
    float sigma_target_    = NAN;
    float R_p0_std_target_ = NAN;
    float R_v0_std_target_ = NAN;

    float acc_noise_floor_sigma_ = ACC_NOISE_FLOOR_SIGMA_DEFAULT;

    // r_p0 = R_p0_coeff * sigma_aw * tau^2, r_v0 = R_v0_coeff * sigma_aw * tau,
    // and tau = tau_coeff * T_z / 2.  All are re-fitted for the wave-band period
    // on the four stationary JONSWAP records, jointly with the corrected
    // accelerometer-bias prior in Kalman3D_Wave_OU_II.h: tau_coeff = 1 is both
    // the documented intent, tau equal to half the zero-crossing period, and the
    // optimum of the scan, while R_p0_coeff fell from 1.6 to 0.6 because the
    // same law now sees a tau two to three times longer.  Along the good ridge
    // the conserved quantity is R_p0_coeff * tau_coeff^2, so the two must be
    // re-fitted together rather than one at a time.
    //
    // The two regularizer coefficients were then re-fitted again, because that
    // fit predates the parity change: the sigma channel moved behind
    // AdaptiveWaveBandPass and now reads lower, so r = c * sigma_aw * tau^k came
    // out below what the records want, and the estimator was over-regularized.
    // R_v0_coeff 1.1 -> 1.3 and R_p0_coeff 0.6 -> 0.65, measured on all eight
    // stationary records at six IMU seeds and on synthesized sea-state
    // transitions.  0.65 is where it stops rather than at the 3D optimum near
    // 0.70: it is the largest position coefficient at which no record's mean
    // vertical error degrades -- seven of the eight improve and the eighth is
    // 1.0003 of the shipped one -- and at 0.70 three records lose, two of them
    // by half a percent.  See docs/ou-ii-pseudo-variance-tuning.md.
    // c_p and c_v of the Empirical law.  Retained under PhysicalMSE so that
    // switching laws at runtime restores the calibrated empirical schedule.
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
    float tau_coeff_   = 1.0f;
    float sigma_coeff_ = 0.85f;

    std::unique_ptr<Kalman3D_Wave_OU_II<float>> mekf_;
    KalmanWaveDirection dir_filter_{2.0f * static_cast<float>(M_PI) * FREQ_GUESS};

    FreqInputLPF     freq_input_lpf_;
    StillnessAdapter freq_stillness_;

    WaveDirectionDetector<float> dir_sign_{0.002f, 0.005f};
    WaveDirection                dir_sign_state_ = UNCERTAIN;
};

template<TrackerType trackerT>
class SeaStateFusion_OU_II {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    struct Config {
        bool with_mag = true;

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
        //
        // Waiting for the good one before reporting anything would put first
        // heading around 105 s, which is not a usable device.  So the first
        // stage locks a provisional reference as soon as the gravity gate
        // allows -- heading and a live filter in roughly 20 s, as before --
        // and the second stage re-learns the reference once the observer has
        // actually converged.  The correction lands long before the scored
        // window opens.
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

        // Magnetometer updates that must land after the filter goes live
        // before the accelerometer-bias gate opens.
        //
        // Accelerometer bias and a tilt error are only weakly separable in
        // waves -- a roll error tips gravity into body Y and reads as a Y bias
        // -- so opening this gate while the attitude is still settling lets the
        // bias absorb the error and hold it.  The proxy bootstrap reaches live
        // far earlier than a staged warmup would, which moves this gate earlier
        // in absolute terms unless it is set to account for that.
        int acc_bias_unlock_mag_updates = 250;

        Eigen::Vector3f sigma_a = Eigen::Vector3f(0.2f, 0.2f, 0.2f);
        Eigen::Vector3f sigma_g = Eigen::Vector3f(0.01f, 0.01f, 0.01f);
        Eigen::Vector3f sigma_m = Eigen::Vector3f(0.3f, 0.3f, 0.3f);

        // The remaining MEKF variances the Kalman3D_Wave_OU_II constructor
        // takes.  They were reachable only through initialize_ext(), which
        // this wrapper never called, so every deployment ran on the header
        // defaults; docs/ou-iii-qmekf-variances.md is the sweep that gauged
        // the same seven for OU-III, and this family carries r_p0 and r_v0
        // where that one carries r_S.  The values here reproduce those
        // defaults exactly.
        //
        //   Pq0      initial attitude-error variance, rad^2.  The proxy
        //            handoff overwrites the attitude block, so this only ever
        //            seeds the covariance the handoff replaces.
        //   Pb0      initial gyro-bias variance, (rad/s)^2.
        //   b0       gyro-bias random-walk variance density, (rad/s)^2/s.
        //   R_p0_noise, R_v0_noise  initial position and velocity pseudo-
        //            measurement variances.  The tuner overwrites both at
        //            Live; they set the pre-Live values.
        float Pq0        = 5e-4f;
        float Pb0        = 1e-6f;
        float b0         = 1e-11f;
        float R_p0_noise = 1.5f;
        float R_v0_noise = 0.3f;
        float gravity_magnitude = g_std;

        // Period-scaled sigma-band shape.  These are dimensionless wave-band
        // ratios except for the absolute safety clamps.  Keeping the ratios
        // fixed is what gives the JONSWAP sigma channel its similarity law.
        float sigma_band_low_ratio  = SIGMA_BAND_LOW_RATIO_DEFAULT;
        float sigma_band_high_ratio = SIGMA_BAND_HIGH_RATIO_DEFAULT;
        float sigma_band_min_hz     = SIGMA_BAND_MIN_HZ_DEFAULT;
        float sigma_band_max_hz     = SIGMA_BAND_MAX_HZ_DEFAULT;

        // Mag-start gate: gravity-direction agreement using current tilt.
        //
        // The residual is measured on the world-frame specific force, averaged
        // over the wave band; see gravityAlignResidualSinWorld() and
        // mag_gravity_align_world_tau_sec below.
        float mag_gravity_align_max_sin   = 0.075f; // sin(deg)
        float mag_gravity_align_hold_sec  = 2.0f;

        // Horizon of the world-frame average the gate is judged on, and how
        // long that average must have been running before its verdict counts.
        //
        // The average has to span whole wave periods for orbital acceleration
        // to cancel out of it, so the horizon is set against the longest swell
        // the device is expected to start up in rather than against the sea it
        // happens to be in -- the frequency tracker is not converged this
        // early, and being conservative here costs settling time rather than
        // accuracy.  12 s covers the band these filters work in.
        //
        // The warmup exists because the average and the observer are seeded
        // from the *same* accelerometer sample.  Until the average has moved
        // off that seed, a small residual only says the two agree about the
        // instant they both started from, which they do by construction even
        // when the boat was mid-wave and both are wrong.
        //
        // It is set so the gate's earliest possible verdict lands with the
        // magnetometer's first eligible sample -- mag_delay_sec less the hold
        // the gate has to serve anyway -- which is the last moment at which it
        // is free.  Beyond that it delays a calm start for nothing: the
        // magnetometer cannot begin averaging before mag_delay_sec however
        // early the gate closes, so any warmup shorter than this buys no time
        // and any warmup longer than it costs time one-for-one.
        float mag_gravity_align_world_tau_sec    = 12.0f;
        float mag_gravity_align_world_warmup_sec = 5.0f;

        float mag_tilt_fallback_sec       = 30.0f;
        float mag_extreme_gyro_dps        = 30.0f; // veto only truly violent motion
        float mag_init_min_mag_norm       = 1e-3f;

        // Mag reference acquisition.
        // Samples are accumulated in the current MEKF world frame:
        //   mag_world = q_mekf_body_to_world * mag_body
        //
        // MagAutoTuner then estimates the yaw gauge of that world frame and returns
        // a gauge-fixed reference:
        //   B_ref = [horizontal_magnitude, 0, vertical]
        //
        // The wrapper removes that same yaw gauge once from the MEKF:
        //   q_new = Rz(-yaw_gauge) * q_old
        //
        // Then normal 3D mag EKF updates run.
        //
        // That average is taken in the estimator's own tilt frame, so whatever
        // tilt error survives the window survives in the reference.  In waves
        // the error is periodic, so what the window has to buy is whole wave
        // periods, not samples: 128 samples is 5.1 s at a 25 Hz mag ODR, short
        // enough to lock in the phase it started on rather than cancel it.
        // 15 s covers a couple of periods across the band these filters work
        // in and captures most of what a much longer window would, at a
        // startup cost of 15 s rather than 40 s.  Held in seconds so it does
        // not silently shorten at a higher ODR.
        int   mag_min_samples    = 128;
        float mag_min_window_sec = 15.0f;
        float mag_max_window_sec = 0.0f;          // no forced timeout
        float mag_sample_dt_sec  = 1.0f / 200.0f;
        
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
        // Ported from OU-III unchanged, because the thing it corrects is not a
        // property of the translational model.  The startup estimate above is a
        // single window, and a single window is where the offset is least
        // identifiable: the excitation is whatever tilt the hull happened to
        // take in fifteen seconds.  This one never closes its accumulation, so
        // the question stops being "can the offset be read out of these fifteen
        // seconds" and becomes "keep watching, and correct when the data
        // finally say something".
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

        bool enable_displacement_detrend = false;
        bool use_custom_displacement_detrend_cfg = false;
        AdaptiveWaveDetrender3D::Config displacement_detrend_cfg{};
    };

    void begin(const Config& cfg) {
        cfg_ = cfg;

        begun_ = true;
        stage_ = Stage::Bootstrap;
        t_ = 0.0f;

        gravity_gate_acc_world_lpf_.reset();
        gravity_gate_world_elapsed_sec_ = 0.0f;
        mag_gravity_good_sec_ = 0.0f;
        mag_gravity_aligned_branch_ = false;
        mag_init_eligible_t0_ = NAN;
        last_mag_sample_t_ = NAN;
      
        mag_ref_set_ = false;

        last_mag_tilt_frame_yaw_rad_ = NAN;
        last_mag_startup_yaw_correction_rad_ = NAN;

        MagAutoTuner::Config mag_cfg;
        mag_cfg.mag_norm_min = cfg_.mag_init_min_mag_norm;
        mag_cfg.min_samples = cfg_.mag_min_samples;
        mag_cfg.min_window_sec = cfg_.mag_min_window_sec;
        mag_cfg.max_window_sec = cfg_.mag_max_window_sec;
        mag_cfg.sample_dt_sec = cfg_.mag_sample_dt_sec;
        mag_cfg.gravity_ref = g_std;
        mag_cfg.enable_quality_weighting = cfg_.mag_enable_quality_weighting;
        mag_cfg.estimate_hard_iron       = cfg_.mag_estimate_hard_iron;
        mag_cfg.min_effective_weight = cfg_.mag_min_effective_weight;
        mag_cfg.acc_norm_rel_soft = cfg_.mag_acc_norm_rel_soft;
        mag_cfg.gyro_soft_dps = cfg_.mag_gyro_soft_dps;
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
        impl_.setMagUpdatesToUnlockAccBias(cfg_.acc_bias_unlock_mag_updates);

        // The provisional reference is deliberately cheap and early, so the
        // accelerometer bias must not be allowed to fit itself to it; see
        // SeaStateFusionFilter_OU_II::setAccBiasHold().
        impl_.setAccBiasHold(cfg_.with_mag && cfg_.mag_refine_enabled);
        impl_.setMagDelaySec(0.0f); // outer wrapper owns startup delay
        impl_.setOnlineTuneWarmupSec(cfg_.online_tune_warmup_sec);
        impl_.setSigmaWaveBandRatios(cfg_.sigma_band_low_ratio,
                                     cfg_.sigma_band_high_ratio);
        impl_.setSigmaWaveBandLimitsHz(cfg_.sigma_band_min_hz,
                                       cfg_.sigma_band_max_hz);

        impl_.initialize_ext(cfg_.sigma_a, cfg_.sigma_g, cfg_.sigma_m,
                             cfg_.Pq0, cfg_.Pb0, cfg_.b0,
                             cfg_.R_p0_noise, cfg_.R_v0_noise,
                             cfg_.gravity_magnitude);
        last_impl_startup_stage_ = impl_.getStartupStage();

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

        last_acc_body_ned_  = acc_body_ned;
        last_gyro_body_ned_ = gyro_body_ned;
        have_last_imu_      = true;

        // The gravity-lock bootstrap is the Mahony observer inside the front
        // end, which runs from the first sample, so there is no phase in which
        // the wrapper withholds IMU data.
        if (stage_ != Stage::Live) {
            // Bootstrap: front end only, MEKF held.
            impl_.updateFrontEnd(dt, gyro_body_ned, acc_body_ned);
        } else {
            impl_.updateTime(dt, gyro_body_ned, acc_body_ned, tempC);
        }

        // Whose tilt the magnetometer gate is judged against.  Before
        // handoff this is the observer's, so the gate measures the
        // attitude that will actually frame the magnetic reference rather
        // than one the MEKF is still converging toward.
        //
        // The residual is taken in that attitude's own world frame rather
        // than in the body frame.  A body-frame average of the specific
        // force is not gravity under way: the hull rolls and pitches
        // through the window, so the orbital term the average is there to
        // remove is smeared across it instead of cancelling.  What the
        // body-frame gate then reports is the sea state, not the levelling
        // error -- on the 8.5 m reference record its residual sits between
        // 0.03 and 0.45 for the whole run against a 0.075 threshold, so
        // the gate simply never closes and startup falls through to its
        // timeout.  Rotating first fixes the frame the average is taken
        // in, orbital acceleration is zero mean there, and the residual
        // becomes the tilt error it was always meant to be: on that same
        // record it settles below 0.05 within about twenty seconds.
        gravity_gate_acc_world_lpf_.step(
            seastate::common::accWorldFromBody(attitudeReferenceQuat_(), acc_body_ned),
            dt,
            cfg_.mag_gravity_align_world_tau_sec);

        const Eigen::Vector3f acc_gate_world_lp = gravity_gate_acc_world_lpf_.state;

        gravity_gate_world_elapsed_sec_ += dt;

        const bool gate_average_warm =
            gravity_gate_world_elapsed_sec_ >= cfg_.mag_gravity_align_world_warmup_sec;

        const float align_sin =
            gate_average_warm
                ? seastate::common::gravityAlignResidualSinWorld(acc_gate_world_lp)
                : 1.0f;

        // The sine residual is the same at an angle and at its supplement,
        // so it accepts an attitude flipped through 180 deg just as readily
        // as the right one.  The branch is the sign of the world down
        // component, and the gate has to carry it: this gate is what
        // certifies the tilt that frames the magnetic reference and that is
        // handed to the MEKF.
        //
        // The branch is deliberately not held behind the warmup.  It is the
        // one part of the certificate an unaveraged sample can answer -- a
        // specific force pointing down in the world frame is a filter that
        // has been seeded upside down, not a wave -- and the handoff
        // timeout is gated on it, so withholding it early would let a
        // stalled startup sit unbranched rather than fail closed.
        const bool aligned_branch =
            seastate::common::gravityAlignedBranchWorld(acc_gate_world_lp);

        mag_gravity_aligned_branch_ = aligned_branch;

        const float gyro_dps = gyro_body_ned.norm() * 57.295779513f;

        // Main gate: gravity-direction agreement only.
        // Gyro only vetoes truly violent motion.
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
            if (mag_gravity_good_sec_ > 10.0f) mag_gravity_good_sec_ = 10.0f;
        } else {
            mag_gravity_good_sec_ = std::max(0.0f, mag_gravity_good_sec_ - 2.0f * dt);
        }

        const Eigen::Vector3f pos_ned_m = impl_.mekf().get_position();
        displacement_up_m_ = Eigen::Vector3f(pos_ned_m.x(), pos_ned_m.y(), -pos_ned_m.z());

        if (cfg_.enable_displacement_detrend) {
            const float wave_hz = impl_.getFreqHz();
            const bool ext_freq_valid =
                isLive() &&
                std::isfinite(wave_hz) &&
                (wave_hz >= displacement_detrender_.config().min_wave_freq_hz) &&
                (wave_hz <= displacement_detrender_.config().max_wave_freq_hz);

            displacement_det_out_ = displacement_detrender_.update(
                displacement_up_m_, dt, wave_hz, ext_freq_valid);
        } else {
            displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};
            displacement_det_out_.input = displacement_up_m_;
            displacement_det_out_.baseline_slow = Eigen::Vector3f::Zero();
            displacement_det_out_.wave_raw = displacement_up_m_;
            displacement_det_out_.wave_clean = displacement_up_m_;
        }

        const auto cur_stage = impl_.getStartupStage();
        if (cur_stage != last_impl_startup_stage_) {
            if (cur_stage == SeaStateFusionFilter_OU_II<trackerT>::StartupStage::Cold) {
                mag_ref_set_ = false;
                mag_auto_tuner_.reset();
                gravity_gate_acc_world_lpf_.reset();
                gravity_gate_world_elapsed_sec_ = 0.0f;
                mag_gravity_good_sec_ = 0.0f;
                mag_gravity_aligned_branch_ = false;
                mag_init_eligible_t0_ = NAN;
                last_mag_sample_t_ = NAN;

                last_mag_tilt_frame_yaw_rad_ = NAN;
                last_mag_startup_yaw_correction_rad_ = NAN;

                if (stage_ != Stage::Live) {
                    // Inner filter already re-locked tilt internally.
                    displacement_up_m_.setZero();
                    displacement_det_out_ = AdaptiveWaveDetrender3D::Output{};
                    if (cfg_.enable_displacement_detrend) {
                        displacement_detrender_.reset(0.0f, 0.0f, 0.0f);
                    }
                }
            }
        
            last_impl_startup_stage_ = cur_stage;
        }

        if (stage_ != Stage::Live) {
            maybeHandOffToMekf_();
        }
    }

    void updateMag(const Eigen::Vector3f& mag_body_ned) {
        if (!begun_ || !cfg_.with_mag) return;
        // There is no withheld stage to wait out: the observer has been
        // levelling since the first sample, and learning north before the MEKF
        // starts is the entire point.  The quality gate below still decides
        // when accumulation may begin.
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
        if (!mag_ref_set_ && t_ < cfg_.proxy_mag_settle_sec) {
            return;
        }

        if (!std::isfinite(mag_init_eligible_t0_)) {
            mag_init_eligible_t0_ = t_;
        }
    
        const bool gravity_trusted = mag_gravity_good_sec_ >= cfg_.mag_gravity_align_hold_sec;
        const bool fallback_ok = (t_ - mag_init_eligible_t0_) >= cfg_.mag_tilt_fallback_sec;
    
        if (!mag_ref_set_) {
            if (!gravity_trusted && !fallback_ok) return;
            if (!have_last_imu_) return;
    
            const float dt_mag =
                std::isfinite(last_mag_sample_t_) && t_ > last_mag_sample_t_
                    ? t_ - last_mag_sample_t_
                    : cfg_.mag_sample_dt_sec;
    
            last_mag_sample_t_ = t_;

            // Accumulate in a tilt frame with yaw removed.
            //
            // Stripping yaw makes the frame invariant to the estimator's
            // arbitrary startup heading, so this leaks no yaw into the learned
            // reference -- q_bw and Rz(psi) q_bw give the same tilt.
            //
            // The tilt comes from the private Mahony observer, which is
            // gyro-propagated through the wave band and corrected below it.
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
                    // This reference was learned in a yaw-stripped tilt frame,
                    // so it carries no estimator heading.  It is a model
                    // parameter, so writing it before the MEKF has been handed
                    // the attitude is safe and leaves it ready to use the
                    // magnetometer from its first live sample.
                    setMagWorldRef_(mag_world_ref_uT);

                    const float yaw_gauge_rad =
                        mag_auto_tuner_.getYawGaugeCorrectionRad();

                    if (std::isfinite(yaw_gauge_rad)) {
                        const float yaw_abs_rad =
                            wrapPi_(-yaw_gauge_rad);

                        // Before handoff there is no MEKF attitude to rewrite;
                        // the gauge is carried to the handoff instead and
                        // composed with the proxy tilt there, so the filter's
                        // very first attitude already has north in it.
                        if (stage_ != Stage::Live) {
                            pending_yaw_abs_rad_ = yaw_abs_rad;

                            last_mag_tilt_frame_yaw_rad_ = wrapPi_(yaw_gauge_rad);
                            last_mag_startup_yaw_correction_rad_ = yaw_abs_rad;
                        } else {
                            Eigen::Quaternionf q_bw =
                                impl_.mekf().quaternion_boat();
                            q_bw.normalize();

                            const Eigen::Quaternionf q_new =
                                boatQuatWithAbsoluteYaw_(
                                    q_bw,
                                    yaw_abs_rad);

                            if (q_new.coeffs().allFinite()) {
                                impl_.mekf().set_quaternion_boat(q_new);

                                last_mag_tilt_frame_yaw_rad_ = wrapPi_(yaw_gauge_rad);
                                last_mag_startup_yaw_correction_rad_ = yaw_abs_rad;
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
                }
            }
        }

        maybeRefineMagReference_(mag_body_ned);
        maybeApplyContinuousHardIron_();

        // Magnetometer corrections go to the MEKF only once it owns the
        // attitude.  Before handoff its state is not the one being solved.
        if (mag_ref_set_ && stage_ == Stage::Live) {
            impl_.updateMag(mag_body_ned - mag_hard_iron_body_uT_);
        }
    }

    bool hasMagNorthLock() const noexcept { return mag_ref_set_; }

    // True once the second-stage acquisition has replaced the provisional
    // reference; see maybeRefineMagReference_().
    bool hasRefinedMagReference() const noexcept { return mag_refine_done_; }
    float magRefineTimeSec() const noexcept { return mag_refine_time_sec_; }

    // Wall-clock marks for the startup sequence, for anyone measuring
    // time-to-first-fix rather than steady-state accuracy.
    float magNorthLockTimeSec() const noexcept { return mag_north_lock_time_sec_; }
    float liveTimeSec() const noexcept { return live_time_sec_; }

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


    bool  isLive() const { return stage_ == Stage::Live; }
    float freqHz() const { return impl_.getFreqHz(); }
    float waveDirectionDeg() const { return impl_.getWaveDirectionDeg(); }

    // Best available boat attitude, BODY -> WORLD (NED).
    //
    // The MEKF holds its initial quaternion until the handoff, so reading it
    // directly during the bootstrap would report a level identity attitude
    // rather than the platform's.  The bootstrap observer's solution is what
    // is actually known at that point, and it is what this returns; after
    // handoff this returns the MEKF's.
    //
    // Heading is only meaningful once hasMagNorthLock() is true (or the build
    // has no magnetometer); before that the bootstrap yaw is arbitrary.  The
    // linear outputs -- displacement, velocity -- stay gated on isLive().
    Eigen::Quaternionf attitudeQuat() const {
        return attitudeReferenceQuat_();
    }

    const Eigen::Vector3f& displacementUpMeters() const { return displacement_up_m_; }
    const AdaptiveWaveDetrender3D::Output& displacementDetrend() const { return displacement_det_out_; }

    SeaStateFusionFilter_OU_II<trackerT>& raw() { return impl_; }
    const SeaStateFusionFilter_OU_II<trackerT>& raw() const { return impl_; }

private:
    enum class Stage {
        Bootstrap,   // proxy owns the attitude, MEKF held
        Live         // MEKF owns the attitude
    };

    struct Vec3LPF {
        Eigen::Vector3f state = Eigen::Vector3f::Zero();
        bool initialized = false;

        void reset() {
            state.setZero();
            initialized = false;
        }

        Eigen::Vector3f step(const Eigen::Vector3f& x, float dt, float tau_sec) {
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

    // The attitude the startup machinery judges and frames things against.
    //
    // Once the MEKF is live it is the answer -- it has the magnetometer, the
    // linear block and the bias states, and the proxy has none of them.
    // Before that the MEKF has nothing to say and the observer does.
    Eigen::Quaternionf attitudeReferenceQuat_() const {
        if (stage_ != Stage::Live) {
            return impl_.startupProxyQuat();
        }
        return impl_.mekf().quaternion_boat();
    }

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
    // adopting its magnitude and dip costs roll accuracy that the offset
    // correction itself does not.
    void maybeApplyContinuousHardIron_() {
        if (!cfg_.mag_continuous_hard_iron) return;
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
    // is live and settled.
    //
    // Both the reference vector and the heading gauge are replaced.  The yaw
    // write is a step, and deliberately so: it is the coarse-to-fine alignment
    // correction, it happens once, and it lands well before the scored window
    // opens.
    void maybeRefineMagReference_(const Eigen::Vector3f& mag_body_ned) {
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
        // the field in it re-derives the error it was meant to remove.
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

    // Hand the bootstrap attitude to the MEKF and start it live.
    //
    // The quality exit needs three things at once: a tilt that has held
    // agreement with gravity long enough to be trusted, a magnetic north gauge
    // (when a magnetometer is fitted), and an operating point the tuner
    // stands behind.  Waiting for all three is what lets the MEKF skip the
    // staged warmup entirely -- there is nothing left for it to converge.
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
        // handoff onto the antipodal branch.  Seeding the MEKF with an attitude
        // that disagrees with measured gravity by more than a right angle is
        // worse than waiting out the extra samples it takes the observer to
        // leave a set it is not attracted to in the first place.
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
    }

    static float wrapPi_(float a) {
        constexpr float PI_F = 3.14159265358979323846f;
        constexpr float TWO_PI_F = 2.0f * PI_F;
        if (!std::isfinite(a)) {
            return NAN;
        }
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
    SeaStateFusionFilter_OU_II<trackerT> impl_{false};

    bool begun_ = false;

    Stage stage_ = Stage::Bootstrap;
    float t_ = 0.0f;

    typename SeaStateFusionFilter_OU_II<trackerT>::StartupStage last_impl_startup_stage_ =
        SeaStateFusionFilter_OU_II<trackerT>::StartupStage::Cold;

    // Last IMU sample for mag-init gating.
    Eigen::Vector3f last_acc_body_ned_  = Eigen::Vector3f::Zero();
    Eigen::Vector3f last_gyro_body_ned_ = Eigen::Vector3f::Zero();
    bool have_last_imu_ = false;

    // Mag-init state.
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

    Vec3LPF gravity_gate_acc_world_lpf_{};
    float   gravity_gate_world_elapsed_sec_ = 0.0f;
    float   mag_gravity_good_sec_ = 0.0f;
    bool    mag_gravity_aligned_branch_ = false;
    float   mag_init_eligible_t0_ = NAN;

};
