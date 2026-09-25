#pragma once

// Defaults shared by the OU-II, OU-III and TFG sea-state orchestrators.
//
// Only one kind of constant belongs here: a property of the marine front end,
// the sensor, the startup/magnetic machinery or the adaptation *mechanics*
// that is the same physical or procedural quantity for every estimator
// behind it.  Each value below had been written out separately in two or
// three wrappers with identical numbers; it is defined once here so the
// wrappers cannot drift apart silently.
//
// Constants that merely happen to share a value today but belong to an
// estimator's model (OU prior bounds, regularizer bounds and coefficients,
// horizontal anisotropy, the r_S/r_p0/r_v0 smoothing multipliers) stay in the
// estimator's own wrapper, and so do the ones that are intentionally
// different; see the "Intentionally different" block at the end, which names
// them so a reader can tell a deliberate difference from drift.

namespace seastate::common::defaults {

// ---------------------------------------------------------------------------
// Sample schedule
// ---------------------------------------------------------------------------

// Nominal IMU sample interval.  The orchestrators and every constant expressed
// per sample below are designed for 200 Hz.
inline constexpr float NOMINAL_IMU_DT_S = 1.0f / 200.0f;

// ---------------------------------------------------------------------------
// Measurement-only front end
// ---------------------------------------------------------------------------

// Estimated pre-band vertical accel noise floor (1 sigma), m/s^2.
// The adaptive sigma band propagates this white-noise variance through its
// actual time-varying coefficients before subtraction.
inline constexpr float ACC_NOISE_FLOOR_SIGMA = 0.12f;

// Bounds of the acceleration-band frequency tracker (wave-direction carrier).
// They do not reach the OU operating point; see the tuning-frequency bounds.
inline constexpr float ACCEL_TRACKER_MIN_FREQ_HZ = 0.2f;
inline constexpr float ACCEL_TRACKER_MAX_FREQ_HZ = 6.0f;

// Floor for the wave-band tuning frequency.  The tracker floor above is far
// too high for a zero-crossing period: 0.03 Hz admits a 33 s swell.  The
// ceiling is estimator-specific and stays with each wrapper.
inline constexpr float MIN_TUNE_FREQ_HZ = 0.03f;

// Wave-band tuning frequency used before WavePeriodEstimator has a value.
//
// It has to be a constant: the whole point of the wave-band source is that no
// part of the adaptation path reads the acceleration-band tracker, and a
// constant is trivially exogenous.  0.2 Hz is a 5 s zero-crossing period, and
// the estimator reports 2.3-8.4 s across the reference family, so the prior is
// never worse than a factor of two off; the estimator replaces it after about
// 50 s in any case.  It is the constant SeaStateFusionFilter_TFG has always
// used for this, which is where the value comes from rather than a fit.
// Sensitivity to it was measured by sweeping it over 0.1-0.4 Hz, which leaves
// every scored 900 s metric unchanged to four decimal places in both OU
// families, because the estimator has replaced it 250 s before the window
// opens.
inline constexpr float TUNE_FREQ_PRIOR_HZ = 0.2f;

// JONSWAP-similar acceleration-variance band.  Away from the absolute safety
// clamps its transfer shape is fixed in f/f_tune, which is the condition the
// sigma_aw similarity argument needs.  The band is a property of the tuner
// rather than of the estimator behind it.
inline constexpr float SIGMA_BAND_LOW_RATIO  = 0.5f;
inline constexpr float SIGMA_BAND_HIGH_RATIO = 4.0f;
inline constexpr float SIGMA_BAND_MIN_HZ     = 0.01f;
inline constexpr float SIGMA_BAND_MAX_HZ     = 6.0f;

// Gains of the private Mahony observer, which levels the vertical channel and
// solves the startup attitude.
//
// two_kp is the accelerometer-to-gyro correction corner: it must stay an order
// of magnitude below the wave band, or the observer levels itself against the
// orbital specific force instead of gravity.
//
// two_ki estimates the gyro bias.  The vertical channel ran at zero for years
// because everything downstream of it is high-passed, but an attitude seed
// keeps whatever static tilt the bias leaves -- about 2b/two_kp, i.e. 0.71 deg
// at 0.05 deg/s -- so the observer that serves both has to estimate it.
inline constexpr float STARTUP_PROXY_TWO_KP = 0.2f;
inline constexpr float STARTUP_PROXY_TWO_KI = 0.02f;

// Front-end accelerometer vibration guard, armed by default.
//
// docs/engine-noise-degradation.md is the measurement.  Machinery vibration
// rectifies into a standing tilt error and from there into a displacement
// offset, and at a routine cruise condition that costs a factor of eight in
// pooled 3-D error; the guard holds every engine condition swept inside a
// factor of 1.6 of the engine-off baseline.
//
// The study was run on OU-III, but the defect is in the attitude loop, not in
// the translational state: vibration reaches the accelerometer, the attitude
// correction wobbles at the vibration frequencies, and the nonlinearity of
// the measurement model rectifies that wobble into a static tilt.  All three
// orchestrators level from the same private Mahony observer at the same gains
// and feed the same accelerometer to their own MEKF, so they inherit the
// defect and the remedy alike.
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
// group delay at full engagement.  A zero cutoff removes the guard entirely
// and restores the unconditioned measurement path.
inline constexpr float ACC_VIBRATION_GUARD_HZ    = 14.0f;
inline constexpr int   ACC_VIBRATION_GUARD_POLES = 2;

// Vibration-aware accelerometer measurement covariance, on by default.
//
// Conditioning removes the machinery the guard can reach; what survives its
// stopband still arrives as measurement error the MEKF does not know about.
// This raises the commanded accelerometer sigma by the guard's own gated
// excess, so the covariance and the measurement describe the same conditions.
//
// The gain is swept in docs/engine-noise-degradation.md, on OU-III.  Pooled
// over the eight stationary records it takes the residual against the
// engine-off baseline from 1.141x to 1.104x at cruise, 1.291x to 1.125x near
// maximum revs, and 1.597x to 1.170x with the sensor at the engine bed, while
// the standing tilt offset falls 1.180 to 0.444 deg and yaw 4.17 to 1.84 deg.
// 0.75 sits at the displacement optimum with margin below the cliff: past
// about 1.25 the accelerometer is de-weighted enough that the wave estimate
// starts leaning on the OU prior instead, and displacement turns back up.
// The trade is between the same two channels in every orchestrator, so the
// same point is carried rather than re-derived.
//
// Zero disables it.  Like the guard, it is driven by a gated excess that is
// identically zero on a quiet installation, so it is bit-transparent there.
inline constexpr float ACC_VIBRATION_RACC_GAIN = 0.75f;

// ---------------------------------------------------------------------------
// Adaptation mechanics (timing and smoothing, not the laws)
// ---------------------------------------------------------------------------

// Common tau/sigma_aw EMA.  A positive sea-period multiplier makes it
// self-similar:
//     tau_ema = ADAPT_TAU_SEA_PERIODS * T_sea,  T_sea = T_z/2.
// The measured default 0.40 is the 0.20*T_z winner from the paired OU-III
// sweep.  Zero selects the legacy fixed-second ADAPT_TAU_SEC path, retained
// for controlled ablations.
inline constexpr float ADAPT_TAU_SEC         = 1.8f;
inline constexpr float ADAPT_TAU_SEA_PERIODS = 0.40f;

// Parameter activation cadence and posterior a_w covariance maintenance
// period.  Candidates are smoothed on every sample; this only throttles when
// the smoothed candidate is committed to the MEKF.
inline constexpr float ADAPT_EVERY_SECS = 0.1f;

// Self-similar drift-regularizer pseudo-measurement cadence
// T_S = (T_0/tau_0) tau.  Every orchestrator historically ran its zero
// pseudo-measurements at T_0 = 15 ms with an initial applied OU time constant
// of tau_0 = 1.1 s, so this ratio preserves that exact operating point and
// scales the cadence with tau thereafter.  A pseudo update cannot occur more
// often than the nominal 200 Hz IMU schedule.  The upper cadence clamp is
// estimator-specific and stays with each wrapper.
inline constexpr float PSEUDO_UPDATE_PERIOD_NOMINAL_S = 0.015f;
inline constexpr float PSEUDO_UPDATE_TAU_NOMINAL_S    = 1.1f;
inline constexpr float PSEUDO_UPDATE_TAU_RATIO =
    PSEUDO_UPDATE_PERIOD_NOMINAL_S / PSEUDO_UPDATE_TAU_NOMINAL_S;
inline constexpr float PSEUDO_UPDATE_PERIOD_MIN_S = NOMINAL_IMU_DT_S;

// ---------------------------------------------------------------------------
// Startup, handoff and magnetic acquisition
// ---------------------------------------------------------------------------

// Magnetometer path held off after power-on.
inline constexpr float MAG_DELAY_SEC = 7.0f;

// Online-tuner warmup the deployed proxy-startup front ends use.  (The inner
// OU filters keep a 5 s compatibility default for direct callers; see
// SeaStateOUFamilyDefaults.h.)
inline constexpr float STARTUP_ONLINE_TUNE_WARMUP_SEC = 10.0f;

// Earliest and latest the proxy bootstrap may hand over.
inline constexpr float PROXY_STARTUP_MIN_SEC     = 8.0f;
inline constexpr float PROXY_STARTUP_TIMEOUT_SEC = 150.0f;

// Covariance seeded at handoff: tilt, yaw with a magnetic north gauge, yaw
// without one.
inline constexpr float PROXY_HANDOFF_TILT_SIGMA_RAD     = 0.035f;   // ~2 deg
inline constexpr float PROXY_HANDOFF_YAW_SIGMA_RAD      = 0.087f;   // ~5 deg
inline constexpr float PROXY_HANDOFF_YAW_SIGMA_FREE_RAD = 1.5708f;  // ~90 deg

// World-frame gravity gate that certifies the startup tilt.
inline constexpr float GRAVITY_GATE_MAX_SIN    = 0.075f;
inline constexpr float GRAVITY_GATE_HOLD_SEC   = 2.0f;
inline constexpr float GRAVITY_GATE_LPF_SEC    = 12.0f;
inline constexpr float GRAVITY_GATE_WARMUP_SEC = 5.0f;
inline constexpr float MAG_TILT_FALLBACK_SEC   = 30.0f;
inline constexpr float MAG_EXTREME_GYRO_DPS    = 30.0f;
inline constexpr float MAG_INIT_MIN_MAG_NORM   = 1e-3f;

// Magnetic reference acquisition window.
inline constexpr int   MAG_MIN_SAMPLES    = 128;
inline constexpr float MAG_MIN_WINDOW_SEC = 15.0f;
inline constexpr float MAG_MAX_WINDOW_SEC = 0.0f;  // no forced timeout
inline constexpr float MAG_SAMPLE_DT_SEC  = 1.0f / 200.0f;
inline constexpr float PROXY_MAG_SETTLE_SEC  = 0.0f;
inline constexpr float MAG_REFINE_WINDOW_SEC = 30.0f;

// MagAutoTuner quality weighting (off: accel/gyro weighting can phase-select
// wave motion).
inline constexpr float MAG_ACC_NORM_REL_SOFT = 0.22f;
inline constexpr float MAG_GYRO_SOFT_DPS     = 45.0f;

// Continuous hard-iron estimation and its slew.
inline constexpr float MAG_HI_MEMORY_SEC           = 600.0f;
inline constexpr float MAG_HI_MODEL_RIDGE          = 5.0e-4f;
inline constexpr float MAG_HI_MODEL_RIDGE_RELATIVE = 0.25f;
inline constexpr float MAG_HI_MIN_INFORMATION      = 0.1f;
inline constexpr float MAG_HI_MIN_EFFECTIVE_WEIGHT = 500.0f;
inline constexpr float MAG_HI_MAX_RESIDUAL_RMS_UT  = 3.0f;
inline constexpr float MAG_HI_MAX_BIAS_FRACTION    = 0.35f;
inline constexpr float MAG_HI_APPLY_FRACTION       = 1.0f;
inline constexpr float MAG_HI_SLEW_TAU_SEC         = 45.0f;

// ---------------------------------------------------------------------------
// Intentionally different (NOT shared; listed so a difference reads as a
// decision rather than as drift)
// ---------------------------------------------------------------------------
//
//                              OU-II        OU-III       TFG
//   max tuning frequency       1.5 Hz       1.2 Hz       1.5 Hz
//   max sigma_aw               6.0          4.0          6.0
//   max pseudo-update period   250 ms       150 ms       250 ms
//   still-water sigma decay    5 s (knob)   1 s          none
//   tau_coeff                  0.95         1.0          1.0
//   sigma_coeff                0.85         0.9          0.8
//   horizontal a_w factor      P = 1.5      S = 1.0      S = 1.0
//   drift-channel EMA mult     3.0 (p, v)   1.5 (r_S)    1.5 (r_S)
//   regularizer law            PhysicalMSE  SpectralMSE  SpectralMSE
//   mag refine start           90 s         90 s         30 s
//
// Each is documented, with its measurement, where it is defined.

}  // namespace seastate::common::defaults
