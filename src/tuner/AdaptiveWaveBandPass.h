#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy
  Released under the MIT License

  AdaptiveWaveBandPass

  Measurement-only, period-scaled band-pass used by the OU-III sea-state
  tuner.  Its two corners are fixed multiples of an externally supplied wave
  frequency, so away from practical absolute clamps the transfer function has
  a fixed shape in normalized frequency f/f_ref.  That is the condition used
  by the JONSWAP similarity theorem for sigma_aw.

  The filter is a one-pole high-pass (implemented as x - LP_low{x}) followed
  by a one-pole low-pass.  In addition to the filtered sample it propagates the
  exact two-state covariance driven by unit-variance discrete white noise.
  Therefore a pre-band bench noise variance can be referred to the current
  adaptive band even while the corners move.
*/

#include <algorithm>
#include <cmath>

class AdaptiveWaveBandPass {
public:
    explicit AdaptiveWaveBandPass(float low_ratio = 0.5f,
                                  float high_ratio = 4.0f,
                                  float min_hz = 0.01f,
                                  float max_hz = 6.0f)
        : low_ratio_(low_ratio),
          high_ratio_(high_ratio),
          min_hz_(min_hz),
          max_hz_(max_hz)
    {
        sanitizeConfig_();
        reset();
    }

    void reset() noexcept {
        lowpass_low_ = 0.0f;
        band_ = 0.0f;
        p00_ = 0.0f;
        p01_ = 0.0f;
        p11_ = 0.0f;
        low_hz_ = NAN;
        high_hz_ = NAN;
        ready_ = false;
    }

    void setRatios(float low_ratio, float high_ratio) noexcept {
        if (!(std::isfinite(low_ratio) && std::isfinite(high_ratio))) return;
        if (!(low_ratio > 0.0f && high_ratio > low_ratio)) return;
        low_ratio_ = low_ratio;
        high_ratio_ = high_ratio;
    }

    void setLimitsHz(float min_hz, float max_hz) noexcept {
        if (!(std::isfinite(min_hz) && std::isfinite(max_hz))) return;
        if (!(min_hz > 0.0f && max_hz > min_hz)) return;
        min_hz_ = min_hz;
        max_hz_ = max_hz;
    }

    float step(float x, float dt, float f_ref_hz) noexcept {
        if (!(std::isfinite(x) && std::isfinite(dt) && std::isfinite(f_ref_hz))) {
            return band_;
        }
        if (!(dt > 0.0f && f_ref_hz > 0.0f)) return band_;

        const float nyquist_guard_hz = 0.45f / dt;
        const float upper_limit = std::min(max_hz_, nyquist_guard_hz);
        if (!(upper_limit > min_hz_)) return band_;

        float low_hz = std::max(min_hz_, low_ratio_ * f_ref_hz);
        low_hz = std::min(low_hz, upper_limit / 1.05f);

        float high_hz = std::min(upper_limit, high_ratio_ * f_ref_hz);
        high_hz = std::max(high_hz, low_hz * 1.05f);
        high_hz = std::min(high_hz, upper_limit);

        if (!(high_hz > low_hz)) return band_;

        const float alpha_low = 1.0f - std::exp(-(2.0f * 3.14159265358979323846f) * low_hz * dt);
        const float alpha_high = 1.0f - std::exp(-(2.0f * 3.14159265358979323846f) * high_hz * dt);
        const float q_low = 1.0f - alpha_low;
        const float q_high = 1.0f - alpha_high;

        // State equations, with z=[LP_low, band]^T and input x:
        //   LP_low' = q_l LP_low + alpha_l x
        //   band'   = -alpha_h q_l LP_low + q_h band
        //             + alpha_h q_l x
        const float low_prev = lowpass_low_;
        const float band_prev = band_;
        lowpass_low_ = q_low * low_prev + alpha_low * x;
        const float high_passed = q_low * (x - low_prev);
        band_ = q_high * band_prev + alpha_high * high_passed;

        // Exact covariance recursion for unit-variance discrete white input.
        // This is intentionally time varying: the same recurrence remains
        // correct when f_ref, and therefore the two coefficients, move.
        const float a00 = q_low;
        const float a10 = -alpha_high * q_low;
        const float a11 = q_high;
        const float b0 = alpha_low;
        const float b1 = alpha_high * q_low;

        const float p00_new = a00 * a00 * p00_ + b0 * b0;
        const float p01_new = a00 * (a10 * p00_ + a11 * p01_) + b0 * b1;
        const float p11_new = a10 * a10 * p00_
                            + 2.0f * a10 * a11 * p01_
                            + a11 * a11 * p11_
                            + b1 * b1;

        p00_ = std::max(0.0f, p00_new);
        p01_ = p01_new;
        p11_ = std::max(0.0f, p11_new);

        low_hz_ = low_hz;
        high_hz_ = high_hz;
        ready_ = true;
        return band_;
    }

    float output() const noexcept { return band_; }
    float lowHz() const noexcept { return low_hz_; }
    float highHz() const noexcept { return high_hz_; }
    float lowRatio() const noexcept { return low_ratio_; }
    float highRatio() const noexcept { return high_ratio_; }
    bool isReady() const noexcept { return ready_; }

    // Output variance produced by unit-variance, zero-mean, discrete white
    // noise at the input, including the current filter transient.
    float whiteNoiseVarianceGain() const noexcept {
        return std::max(0.0f, p11_);
    }

private:
    void sanitizeConfig_() noexcept {
        if (!(std::isfinite(low_ratio_) && low_ratio_ > 0.0f)) low_ratio_ = 0.5f;
        if (!(std::isfinite(high_ratio_) && high_ratio_ > low_ratio_)) high_ratio_ = 4.0f;
        if (!(std::isfinite(min_hz_) && min_hz_ > 0.0f)) min_hz_ = 0.01f;
        if (!(std::isfinite(max_hz_) && max_hz_ > min_hz_)) max_hz_ = 6.0f;
    }

    float low_ratio_ = 0.5f;
    float high_ratio_ = 4.0f;
    float min_hz_ = 0.01f;
    float max_hz_ = 6.0f;

    float lowpass_low_ = 0.0f;
    float band_ = 0.0f;

    // Covariance of [LP_low, band] for unit white input.
    float p00_ = 0.0f;
    float p01_ = 0.0f;
    float p11_ = 0.0f;

    float low_hz_ = NAN;
    float high_hz_ = NAN;
    bool ready_ = false;
};
