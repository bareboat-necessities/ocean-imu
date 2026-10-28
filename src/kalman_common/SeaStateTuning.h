#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy
*/

// Tuning envelope, integral-regularizer (S = 0) law and committed tune of the
// deployed OU-III orchestrator.  These are OU-III's estimator-model constants;
// the OU-II and TFG values differ (see "Intentionally different" in
// SeaStateFusionDefaults.h).  They are defined here rather than in
// SeaStateFusionFilter_OU_III.h, which re-exports them into ocean_imu::ou3,
// so that the CARMA research wrappers, which keep OU-III's schedule and
// replace only its motion prior, use the same values without including the
// OU-III estimator.

#include <cstdint>

#include "kalman_common/SeaStateFusionDefaults.h"

namespace seastate::common::tuning {

// Ceiling for the wave-band tuning frequency: a safety limit, not a tuning
// surface.  MAX_FREQ_HZ is the tracker's bound and does not belong in the
// adaptation path; setFreqBounds() is a wave-direction knob and must not move
// the OU operating point.  1.2 Hz is a 0.83 s zero-crossing period, shorter
// than any sea a hull responds to.  OU-II and TFG keep 1.5 Hz.
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

// Smoothing horizon of the r_S channel, in units of tau_target.  Shorter
// tracks transitions faster; longer filters tau noise, which r_S amplifies by
// at least the third power (it grows at least as fast as tau^3).  The common
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


// Integral-regularizer adaptation laws.  Derivations and the calibration
// evidence are in the OU-III paper (Sec. "Implemented SpectralMSE schedule").
// r_a ~ R_a*dt is the acceleration-noise density of the reduced scalar model
// and q_eff = 2 r_a (1 - 1/sqrt(1+zeta)), zeta = 2 sigma_aw^2 tau / r_a, its
// posterior acceleration-noise intensity.
//
//   Cubic            r_S,base = C_R sqrt(R_a) tau^3, renormalized by
//                    sqrt(T_0/T_S).  C_R places the normalized regularizer
//                    corner kappa = omega_R tau: C_R = sqrt(2h/T_0)/kappa^3.
//                    No transcendental per tuner update, so it is the
//                    supported low-cost configuration for targets without
//                    hardware powf; it needs C_R rather than C_J.
//   StrongRiccati    the strong-observation branch q_eff = 2 r_a:
//                    r_S = sqrt(2 r_a) tau^3 / (sqrt(T_S) kappa^3), the same
//                    schedule as Cubic at the analytical C_R.
//   PosteriorRiccati the full q_eff law; r_S ~ sigma_aw tau^(7/2)/sqrt(T_S)
//                    as zeta -> 0 and StrongRiccati as zeta -> infinity.
//   SpectralMSE      the pole minimizing drift plus wave-distortion error,
//                        r_S ~ C_J q_eff^(1/14) sigma_a^(6/7) tau^(24/7)
//                              / sqrt(T_S).
//                    DEPLOYED DEFAULT.
//
// The Riccati laws take a one-parameter amplitude tilt (p = 0: none; p = 1:
// a normalized sigma_aw factor).  It does not apply to SpectralMSE.
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
constexpr float R_S_ACCEL_NOISE_DENSITY_DEFAULT = 0.0148f * 0.0148f * defaults::NOMINAL_IMU_DT_S;
// Analytical pole-placement value of C_R, Eq. (adapt-cR-kappa) of the OU-III
// paper: C_R = sqrt(2h/T_S,0) / kappa^3 with h = 5 ms (200 Hz), T_S,0 = 15 ms
// and kappa = 0.3627.  At this value the Cubic base, after cadence
// renormalization, is identical to the StrongRiccati law.
constexpr float R_S_COEFF_ANALYTICAL_REFERENCE = 17.112f;

// Coefficient C_J of the SpectralMSE law,
//     r_S = C_J q_eff^(1/14) sigma_a,B^(6/7) tau^(24/7) / sqrt(T_S).
// It absorbs the dimensionless spectral moment int x^-8 Phi_a(x) dx of the
// self-similar acceleration spectrum, raised to 3/7; the default is that
// analytical value, not a fit.
constexpr float R_S_MSE_COEFF_DEFAULT = 0.0538f;

struct TuneState {
    float tau_applied   = 1.1f;    // s
    float sigma_applied = 1e-2f;   // m/s²
    float RS_applied    = 0.5f;    // m*s
};

}  // namespace seastate::common::tuning
