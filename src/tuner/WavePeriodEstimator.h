#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy

  WavePeriodEstimator - zero-crossing wave period from vertical acceleration.

  Why this exists
  ---------------
  The OU time constant and the integral regularization scale are properties of
  the *wave* band: they have to follow the sea state.  A frequency tracker
  running on acceleration cannot supply that.  An ocean acceleration spectrum is
  the elevation spectrum weighted by (2*pi*f)^4, so its apparent frequency sits
  well above the spectral peak and, because the weighting fights the spectral
  roll-off, barely moves as the sea grows: across a JONSWAP family spanning
  H_s = 0.27 m to 8.5 m the elevation peak moves 0.335 -> 0.085 Hz while the
  acceleration-band frequency stays near 0.42-0.55 Hz.

  What it estimates
  -----------------
  The zero-crossing period T_z = 2*pi*sqrt(m0/m2) of the surface elevation,
  where m_n are the elevation spectral moments.  The elevation proxy is one
  leaky integration H(s) = 1/(s + lambda) further from the acceleration than the
  velocity proxy is, so for a narrow band at omega

      sigma_v / sigma_eta = sqrt(omega^2 + lambda^2),

  and omega^2 = (sigma_v/sigma_eta)^2 - lambda^2 removes the leak exactly.  The
  relation holds for *any* filtering the two proxies share, because they differ
  by exactly one integrator; a broadband input gives the mean of the same
  relation weighted by whatever the shared response passes.

  That is what the two first-order high-pass stages in front of the integrators
  are for.  Double integration weights a spectrum by 1/omega^4, so the sub-band
  energy a real accelerometer carries - bias random walk, gravity leaking
  through attitude error - dominates the elevation proxy and drags the period
  estimate toward the leak corner.  On the reference records the leak alone
  reported 12-16 s for every sea state, against a true 2.5-8.6 s, and in the
  wrong order.  Two high-pass stages at the same corner make the shared response
  s^2/(s+lambda)^2, which cancels the 1/omega^4 growth below the band while
  costing under 6 percent in gain at the lowest wave frequency of interest.  The
  ratio relation above is untouched, because both proxies see the same stages.

  Independence
  ------------
  The integrators are driven by the accelerometer only, never by the filter's
  own displacement.  Deriving the period from the filter's displacement instead
  would close a positive feedback loop: the integral regularizer high-passes
  displacement, which raises the apparent frequency, which shortens tau, which
  strengthens the regularizer.
*/

#include <algorithm>
#include <cmath>

class WavePeriodEstimator {
public:
    // high_pass_hz : corner shared by the two high-pass stages and the two
    //                integrator leaks.  It must sit below the lowest wave
    //                frequency of interest; 0.03 Hz (33 s) covers any swell a
    //                small vessel meets and still rejects instrument drift.
    // horizon_periods / min_horizon_sec / max_horizon_sec : averaging window
    //                for the variances, expressed in wave periods and clamped.
    explicit WavePeriodEstimator(float high_pass_hz = 0.02f,
                                 float horizon_periods = 8.0f,
                                 float min_horizon_sec = 20.0f,
                                 float max_horizon_sec = 180.0f)
        : lambda_(2.0f * 3.14159265358979323846f * std::max(1e-4f, high_pass_hz)),
          horizon_periods_(std::max(1.0f, horizon_periods)),
          min_horizon_sec_(std::max(1.0f, min_horizon_sec)),
          max_horizon_sec_(std::max(min_horizon_sec_, max_horizon_sec))
    {
        reset();
    }

    void reset() {
        accel_prev_ = 0.0f;
        high_pass_1_ = 0.0f;
        high_pass_1_prev_ = 0.0f;
        high_pass_2_ = 0.0f;
        velocity_ = 0.0f;
        elevation_ = 0.0f;
        velocity_mean_ = 0.0f;
        velocity_sq_ = 0.0f;
        elevation_mean_ = 0.0f;
        elevation_sq_ = 0.0f;
        weight_ = 0.0f;
        elapsed_sec_ = 0.0f;
        period_sec_ = NAN;
    }

    // vertical_accel_ms2 : up-positive vertical acceleration proxy, gravity
    // already removed.  Offsets and drift are rejected by the high-pass stages.
    void update(float dt_sec, float vertical_accel_ms2) {
        if (!(dt_sec > 0.0f) || !std::isfinite(dt_sec) ||
            !std::isfinite(vertical_accel_ms2)) {
            return;
        }

        // Exponential steps: exact for a step-held input, so nothing here
        // depends on the sample rate.
        const float decay = std::exp(-lambda_ * dt_sec);
        const float gain = (lambda_ > 1e-9f) ? ((1.0f - decay) / lambda_) : dt_sec;

        // Two shared high-pass stages, y <- decay*(y + x - x_prev).
        const float stage1 = decay * (high_pass_1_ + vertical_accel_ms2 - accel_prev_);
        accel_prev_ = vertical_accel_ms2;
        const float stage2 = decay * (high_pass_2_ + stage1 - high_pass_1_prev_);
        high_pass_1_prev_ = stage1;
        high_pass_1_ = stage1;
        high_pass_2_ = stage2;

        velocity_ = decay * velocity_ + gain * stage2;
        elevation_ = decay * elevation_ + gain * velocity_;

        elapsed_sec_ += dt_sec;

        // Hold the integrators for a few leak time constants before believing
        // their statistics; the startup transient is pure bias otherwise.
        const float settle_sec = 6.0f / lambda_;
        if (elapsed_sec_ < settle_sec) return;

        const float horizon = averaging_horizon_sec();
        const float alpha = 1.0f - std::exp(-dt_sec / horizon);

        weight_ = (1.0f - alpha) * weight_ + alpha;
        velocity_mean_ = (1.0f - alpha) * velocity_mean_ + alpha * velocity_;
        velocity_sq_ = (1.0f - alpha) * velocity_sq_ + alpha * velocity_ * velocity_;
        elevation_mean_ = (1.0f - alpha) * elevation_mean_ + alpha * elevation_;
        elevation_sq_ = (1.0f - alpha) * elevation_sq_ + alpha * elevation_ * elevation_;

        if (!(weight_ > 1e-3f)) return;

        const float velocity_mean = velocity_mean_ / weight_;
        const float elevation_mean = elevation_mean_ / weight_;
        const float velocity_var =
            std::max(0.0f, velocity_sq_ / weight_ - velocity_mean * velocity_mean);
        const float elevation_var =
            std::max(0.0f, elevation_sq_ / weight_ - elevation_mean * elevation_mean);

        if (!(elevation_var > 1e-12f) || !(velocity_var > 1e-12f)) return;

        // omega^2 = (sigma_v/sigma_eta)^2 - lambda^2 inverts the leak exactly
        // for a narrow band and is the moment-weighted mean otherwise.
        const float ratio_sq = velocity_var / elevation_var;
        const float omega_sq = ratio_sq - lambda_ * lambda_;
        if (!(omega_sq > 1e-8f)) return;

        const float period = 2.0f * 3.14159265358979323846f / std::sqrt(omega_sq);
        if (std::isfinite(period) && period > 0.0f) {
            period_sec_ = period;
        }
    }

    // Zero-crossing period of the elevation [s]; NaN until settled.
    float getPeriodSec() const { return period_sec_; }

    float getFrequencyHz() const {
        return (period_sec_ > 1e-6f) ? (1.0f / period_sec_) : NAN;
    }

    bool isReady() const {
        return std::isfinite(period_sec_) && weight_ > 0.5f;
    }

    // Elevation standard deviation of the band-limited proxy [m].  Not used for
    // tuning, but it makes the estimator inspectable from logs.
    float getElevationStd() const {
        if (!(weight_ > 1e-3f)) return NAN;
        const float mean = elevation_mean_ / weight_;
        return std::sqrt(std::max(0.0f, elevation_sq_ / weight_ - mean * mean));
    }

private:
    float averaging_horizon_sec() const {
        const float period = std::isfinite(period_sec_) ? period_sec_ : 6.0f;
        return std::min(max_horizon_sec_,
                        std::max(min_horizon_sec_, horizon_periods_ * period));
    }

    float lambda_;
    float horizon_periods_;
    float min_horizon_sec_;
    float max_horizon_sec_;

    float accel_prev_ = 0.0f;
    float high_pass_1_ = 0.0f;
    float high_pass_1_prev_ = 0.0f;
    float high_pass_2_ = 0.0f;
    float velocity_ = 0.0f;
    float elevation_ = 0.0f;
    float velocity_mean_ = 0.0f;
    float velocity_sq_ = 0.0f;
    float elevation_mean_ = 0.0f;
    float elevation_sq_ = 0.0f;
    float weight_ = 0.0f;
    float elapsed_sec_ = 0.0f;
    float period_sec_ = NAN;
};
