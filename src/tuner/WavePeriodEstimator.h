#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy

  WavePeriodEstimator - zero-crossing wave period from vertical acceleration.

  The physical estimator is the spectral-moment identity

      T_z = 2*pi*sqrt(m0/m2),

  implemented through the variance ratio of two proxies that differ by one
  leaky integration.  The exponentially weighted moment estimates remain the
  primary statistic; period/frequency are derived from those moments rather
  than averaged as raw reciprocal samples.

  Two time scales are deliberately separated:
    * the moment horizon is expressed in wave periods, so that as the sea grows
      the estimator keeps averaging over a comparable number of cycles rather
      than a fixed wall-clock window.  The absolute floor below still wins on
      short seas: at the retuned default of 4 periods, min_horizon_sec = 20 s
      binds for every T_z under 5 s, so across the eight reference seas
      (T_z ~= 2.3..8.4 s) the realized horizon runs about 8 periods on the
      shortest chop down to 4 on the longest swell.  The retuning sweep scored
      the period counts with this floor in place, so moving the floor is itself
      a retune and needs a re-gauge, not a tweak;
    * a canonical log-period state smooths the nonlinear moment-ratio output.
      The moment horizon follows this slowly moving positive state, not the raw
      instantaneous ratio, avoiding the old direct self-modulation of estimator
      bandwidth by its noisiest output.

  Because log(f) = -log(T), the canonical state is invariant to whether a
  consumer wants period or frequency: getPeriodSec() and getFrequencyHz() are
  exact reciprocals of one another.
*/

#include <algorithm>
#include <cmath>

#include "tuner/SeaStateAdaptationLimits.h"

class WavePeriodEstimator {
public:
    // high_pass_hz: corner shared by the two high-pass stages and integrator
    // leaks.  It must sit below the wave band of interest.
    // moment_horizon_periods: EW moment time constant in T_z units.
    // log_smoothing_periods: canonical log(T_z) EMA time constant in T_z units.
    // min/max_horizon_sec: absolute guards on the moment horizon only.
    explicit WavePeriodEstimator(float high_pass_hz = 0.02f,
                                 float moment_horizon_periods = 4.0f,
                                 float log_smoothing_periods = 0.05f,
                                 float min_horizon_sec = 20.0f,
                                 float max_horizon_sec = 180.0f)
        : lambda_(2.0f * 3.14159265358979323846f * std::max(1e-4f, high_pass_hz)),
          moment_horizon_periods_(std::max(1.0f, moment_horizon_periods)),
          log_smoothing_periods_(std::max(0.0f, log_smoothing_periods)),
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
        raw_period_sec_ = NAN;
        log_period_sec_ = NAN;
        last_moment_horizon_sec_ = min_horizon_sec_;
        last_log_horizon_sec_ = 0.0f;
    }

    void setMomentHorizonPeriods(float periods) {
        if (std::isfinite(periods) && periods > 0.0f) {
            moment_horizon_periods_ = std::max(1.0f, periods);
        }
    }
    float getMomentHorizonPeriods() const { return moment_horizon_periods_; }

    void setLogSmoothingPeriods(float periods) {
        if (std::isfinite(periods) && periods >= 0.0f) {
            log_smoothing_periods_ = periods;
        }
    }
    float getLogSmoothingPeriods() const { return log_smoothing_periods_; }

    void setMomentHorizonBounds(float min_sec, float max_sec) {
        if (!(std::isfinite(min_sec) && std::isfinite(max_sec))) return;
        if (!(min_sec > 0.0f && max_sec >= min_sec)) return;
        min_horizon_sec_ = min_sec;
        max_horizon_sec_ = max_sec;
    }

    float getMomentHorizonSec() const { return last_moment_horizon_sec_; }
    float getLogSmoothingHorizonSec() const { return last_log_horizon_sec_; }

    // vertical_accel_ms2: up-positive vertical acceleration proxy, gravity
    // already removed.  Offsets and drift are rejected by the high-pass stages.
    void update(float dt_sec, float vertical_accel_ms2) {
        if (!(dt_sec > 0.0f) || !std::isfinite(dt_sec) ||
            !std::isfinite(vertical_accel_ms2)) {
            return;
        }

        // Exact first-order discretization for a step-held input.
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

        const float horizon = moment_horizon_sec_();
        last_moment_horizon_sec_ = horizon;
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

        const float raw_period =
            2.0f * 3.14159265358979323846f / std::sqrt(omega_sq);
        if (!(std::isfinite(raw_period) && raw_period > 0.0f)) return;

        raw_period_sec_ = raw_period;
        update_log_period_(dt_sec, raw_period);
    }

    // Canonical zero-crossing period [s].  It is positive by construction and
    // NaN until the moment ratio has produced its first valid value.
    float getPeriodSec() const {
        return std::isfinite(log_period_sec_) ? std::exp(log_period_sec_) : NAN;
    }

    // Exact reciprocal of the same canonical log state.  No second frequency
    // smoother exists in this estimator chain.
    float getFrequencyHz() const {
        return std::isfinite(log_period_sec_) ? std::exp(-log_period_sec_) : NAN;
    }

    // Raw moment-ratio period before log-domain smoothing, for diagnostics and
    // controlled ablations only.
    float getRawPeriodSec() const { return raw_period_sec_; }

    bool isReady() const {
        return std::isfinite(log_period_sec_) && weight_ > 0.5f;
    }

    // Elevation standard deviation of the band-limited proxy [m].  Not used for
    // tuning, but it makes the estimator inspectable from logs.
    float getElevationStd() const {
        if (!(weight_ > 1e-3f)) return NAN;
        const float mean = elevation_mean_ / weight_;
        return std::sqrt(std::max(0.0f, elevation_sq_ / weight_ - mean * mean));
    }

private:
    float period_for_horizon_sec_() const {
        const float period = getPeriodSec();
        return (std::isfinite(period) && period > 0.0f) ? period : 6.0f;
    }

    float moment_horizon_sec_() const {
        const float requested = moment_horizon_periods_ * period_for_horizon_sec_();
        return std::min(max_horizon_sec_, std::max(min_horizon_sec_, requested));
    }

    void update_log_period_(float dt_sec, float raw_period_sec) {
        const float log_raw = std::log(raw_period_sec);
        if (!std::isfinite(log_period_sec_)) {
            // Initialize from the first valid physical moment ratio; debiasing a
            // second smoother here would only extend startup without adding
            // information.
            log_period_sec_ = log_raw;
            last_log_horizon_sec_ = 0.0f;
            return;
        }

        if (!(log_smoothing_periods_ > 0.0f)) {
            log_period_sec_ = log_raw;
            last_log_horizon_sec_ = dt_sec;
            return;
        }

        const float sea_period = std::exp(log_period_sec_);
        const float requested = log_smoothing_periods_ * sea_period;
        const float horizon =
            seastate::tuner::limits::clampDynamicEmaHorizonSec(requested, dt_sec);
        last_log_horizon_sec_ = horizon;
        const float alpha = 1.0f - std::exp(-dt_sec / horizon);
        log_period_sec_ += alpha * (log_raw - log_period_sec_);
    }

    float lambda_;
    float moment_horizon_periods_;
    float log_smoothing_periods_;
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
    float raw_period_sec_ = NAN;
    float log_period_sec_ = NAN;
    float last_moment_horizon_sec_ = 0.0f;
    float last_log_horizon_sec_ = 0.0f;
};
