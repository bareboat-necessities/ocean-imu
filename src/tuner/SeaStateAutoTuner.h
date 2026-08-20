#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy

  SeaStateAutoTuner — minimal online acceleration statistics.

  The wave-frequency input is already the canonical, measurement-only output of
  WavePeriodEstimator.  This class therefore does not estimate or smooth that
  time scale again.  It stores the latest supplied frequency only so the
  period-scaled acceleration band and variance horizon can use the same
  one-sample-predictable value.

  Acceleration variance is the ordinary exponentially weighted second central
  moment

      sigma_a^2 = E_w[a^2] - E_w[a]^2,

  with identical weights for the first and second raw moments.  The previous
  second EMA of this already-smoothed variance has been removed; any desired
  downstream parameter slew belongs to the filter scheduler, not to the
  statistical variance estimator.
*/

#include <algorithm>
#include <cmath>

#include "tuner/SeaStateAdaptationLimits.h"

struct DebiasedEMA {
    float value  = 0.0f;
    float weight = 0.0f;
    inline void reset() { value = 0.0f; weight = 0.0f; }
    inline void update(float x, float alpha) {
        value  = (1.0f - alpha) * value + alpha * x;
        weight = (1.0f - alpha) * weight + alpha;
    }
    inline float get() const {
        return (weight > 1e-12f) ? value / weight : 0.0f;
    }
    inline bool isReady() const { return weight > 1e-6f; }
};

class SeaStateAutoTuner {
public:
    // K_periods: dimensionless time constant of the acceleration-moment EMA in
    // periods of the canonical wave-frequency input.
    //
    // The second constructor argument is retained only for source compatibility
    // with older call sites; frequency smoothing no longer occurs here.
    explicit SeaStateAutoTuner(float K_periods_ = 2.0f,
                               float /*legacy_tau_freq_sec*/ = 1.0f)
        : K_periods(std::max(1e-3f, K_periods_)) {
        reset();
    }

    inline void reset() {
        tau_var_sec = 0.0f;
        frequency_hz = NAN;
        A_mean.reset();
        A_sq.reset();
    }

    inline void update(float dt_s, float accel, float f_input_hz) {
        if (!(dt_s > 0.0f) || !std::isfinite(accel) ||
            !std::isfinite(f_input_hz)) {
            return;
        }

        float f_eff = f_input_hz;
        f_eff = std::max(f_min_hz, std::min(f_max_hz, f_eff));
        frequency_hz = f_eff;

        const float sea_time_sec =
            seastate::tuner::limits::clampDynamicEmaTimeScaleSec(0.5f / f_eff);
        const float T_eff = 2.0f * sea_time_sec;  // T_z = 2 T_sea

        const float tau_var_requested = std::max(
            tau_var_min_sec, std::min(tau_var_max_sec, K_periods * T_eff));
        tau_var_sec = seastate::tuner::limits::clampDynamicEmaHorizonSec(
            tau_var_requested, dt_s);

        const float alpha_var = 1.0f - std::exp(-dt_s / tau_var_sec);
        A_mean.update(accel, alpha_var);
        A_sq.update(accel * accel, alpha_var);
    }

    inline float getAccelVariance() const {
        if (!(A_mean.isReady() && A_sq.isReady())) return 0.0f;
        const float mu = A_mean.get();
        return std::max(0.0f, A_sq.get() - mu * mu);
    }

    inline float getAccelStd() const {
        return std::sqrt(getAccelVariance());
    }

    inline float getFrequencyHz() const { return frequency_hz; }

    inline float getPeriodSec() const {
        return (std::isfinite(frequency_hz) && frequency_hz > 1e-9f)
            ? (1.0f / frequency_hz)
            : 0.0f;
    }

    inline bool isReady() const { return isVarReady() && isFreqReady(); }
    inline bool isFreqReady() const {
        return std::isfinite(frequency_hz) && frequency_hz > 0.0f;
    }
    inline bool isVarReady() const { return A_mean.isReady() && A_sq.isReady(); }

    // The requested horizon is now the actual first/second-moment time constant;
    // there is no hidden second variance EMA doubling its memory.
    inline float getVarianceHorizonSec() const { return tau_var_sec; }

    inline void setKPeriods(float k) {
        if (std::isfinite(k) && k > 0.0f) K_periods = k;
    }
    inline float getKPeriods() const { return K_periods; }

    inline void setVarianceHorizonBounds(float min_s, float max_s) {
        if (!(std::isfinite(min_s) && std::isfinite(max_s))) return;
        if (!(min_s > 0.0f && max_s >= min_s)) return;
        tau_var_min_sec = min_s;
        tau_var_max_sec = max_s;
    }
    inline float getVarianceHorizonMinSec() const { return tau_var_min_sec; }
    inline float getVarianceHorizonMaxSec() const { return tau_var_max_sec; }

    inline void setFrequencyBounds(float min_hz, float max_hz) {
        if (!(std::isfinite(min_hz) && std::isfinite(max_hz))) return;
        if (!(min_hz > 0.0f && max_hz >= min_hz)) return;
        f_min_hz = min_hz;
        f_max_hz = max_hz;
    }

    // Compatibility shims for older wrappers/ablation code.  Frequency
    // smoothing has intentionally moved upstream into WavePeriodEstimator's
    // canonical log-period state, so these knobs no longer alter the tuner.
    inline void setTauFreq(float) {}
    inline void setFrequencySmoothingSeaPeriods(float) {}
    inline float getFrequencySmoothingSeaPeriods() const { return 0.0f; }
    inline float getFrequencySmoothingHorizonSec() const { return 0.0f; }

private:
    float K_periods = 2.0f;
    float tau_var_min_sec = 0.3f;
    float tau_var_max_sec = 60.0f;
    float f_min_hz = 0.05f;
    float f_max_hz = 5.0f;
    float tau_var_sec = 0.0f;
    float frequency_hz = NAN;

    DebiasedEMA A_mean;
    DebiasedEMA A_sq;
};

#ifdef SEA_STATE_TUNER_TEST
#include <iostream>

static inline void SeaStateAutoTuner_test() {
    constexpr float Fs   = 200.0f;
    constexpr float DT   = 1.0f / Fs;
    constexpr float F_HZ = 0.5f;
    const float omega = 2.0f * 3.14159265358979323846f * F_HZ;

    SeaStateAutoTuner tuner(1.5f);

    float t = 0.0f;
    for (int n = 0; n < int(20.0f / DT); ++n) {
        t += DT;
        const float a = -std::sin(omega * t) * omega * omega;
        tuner.update(DT, a, F_HZ);
    }

    std::cerr << "[AutoTuner] sigma_a=" << tuner.getAccelStd()
              << " m/s2, f=" << tuner.getFrequencyHz()
              << " Hz\n";
}
#endif
