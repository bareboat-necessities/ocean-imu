#pragma once

/*
  Copyright (c) 2025-2026  Mikhail Grushinskiy

  SeaStateAutoTuner — minimal online acceleration stats

  Purpose:
    • Estimate acceleration variance σ_a² (time-domain EWMA)
    • Smooth externally provided frequency f_in (Hz) via EMA

  Notes:
    • K_periods is a dimensionless factor controlling the variance horizon:
        τ_var_dyn ≈ K_periods * T_eff,
      where T_eff = 1 / f_eff and f_eff is the smoothed frequency.
    • The horizon is expressed in *periods of the frequency it is fed*, so what
      that frequency means decides what the horizon means.  Fed the
      acceleration-band tracker it was a few seconds whatever the sea did,
      because that band barely moves with sea state; fed the wave-band
      zero-crossing period it becomes a genuine multiple of the wave period and
      stretches from about 5 s on a short chop to 17 s on a developed swell.
      The clamps below therefore have to admit that whole span, and the
      defaults are documented in docs/ou-sigma-horizon.md.
*/

#include <cmath>
#include <algorithm>

// Where the tuning frequency the OU adaptation runs on comes from.
//
// WaveBand is the deployed choice: the operating point is a property of the
// wave band, so it is read from WavePeriodEstimator alone and, before that
// estimator has settled, from a fixed wave-band prior.  Nothing in the
// adaptation path then reads the acceleration-band frequency tracker, which
// stays where it is needed -- as the carrier of the wave-direction
// demodulator.
//
// TrackerFallback is the previous behaviour, kept as an ablation: the same
// wave-band frequency once the period estimator is ready, but the tracker's
// output for the first minute or so of every run.
//
// WaveBandGated separates the two things WaveBand changes at once.  It drops
// the tracker exactly like WaveBand but keeps the legacy readiness gate, so
// comparing the three attributes the effect to "no tracker" or to "adopt the
// period estimator as soon as it has a value" rather than to their sum.
enum class TunerFrequencySource {
    WaveBand,         // wave-period estimator as soon as it has a value, then a prior
    WaveBandGated,    // ablation: same, but only once the estimator reports ready
    TrackerFallback,  // legacy: acceleration-band tracker until the estimator is ready
};

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
    // K_periods: dimensionless horizon in periods.
    //   τ_var_dyn ≈ K_periods * T_eff,  T_eff = 1 / f_eff
    // Default K_periods ≈ 2.0 → ~6.0 periods for ~95% response.
    explicit SeaStateAutoTuner(float K_periods_   = 2.0f,
                               float tau_freq_sec = 1.0f)   // frequency smoothing horizon (seconds)
    : K_periods(K_periods_), tau_freq(tau_freq_sec) {
        reset();
    }

    inline void reset() {
        last_dt_freq = -1.0f;
        alpha_freq = 0.0f;
        tau_var_sec = 0.0f;
        A_mean.reset(); A_sq.reset(); A_var.reset();
        Freq_smoothed.reset();
    }

    // main update
    inline void update(float dt_s, float accel, float f_input_hz) {
        if (!(dt_s > 0.0f) || !std::isfinite(accel) || !std::isfinite(f_input_hz))
            return;

        // Update alpha for frequency smoothing (fixed horizon in seconds)
        updateAlphaFreq(dt_s);

        // Smooth incoming frequency first
        Freq_smoothed.update(f_input_hz, alpha_freq);

        // Effective frequency for variance horizon (use smoothed if ready)
        float f_eff = Freq_smoothed.isReady() ? Freq_smoothed.get() : f_input_hz;
        // Clamp to avoid insane horizons at tiny/zero frequency
        if (!std::isfinite(f_eff)) f_eff = f_min_hz;
        f_eff = std::max(f_min_hz, std::min(f_max_hz, f_eff));

        const float T_eff = 1.0f / f_eff;

        // Dynamic variance time constant: a few periods
        tau_var_sec = std::max(tau_var_min_sec,
                               std::min(tau_var_max_sec, K_periods * T_eff));

        // Compute alpha for variance based on dynamic tau
        const float alpha_var = 1.0f - std::exp(-dt_s / tau_var_sec);

        // Time-domain EWMA variance
        A_mean.update(accel, alpha_var);
        A_sq.update(accel * accel, alpha_var);
        const float mu       = A_mean.get();
        const float var_inst = std::max(0.0f, A_sq.get() - mu * mu);
        A_var.update(var_inst, alpha_var);
    }

    // accessors
    inline float getAccelVariance() const { return A_var.get(); }       // σ_a²
    inline float getAccelStd()      const { return std::sqrt(std::max(0.0f, A_var.get())); }
    inline float getFrequencyHz()   const { return Freq_smoothed.get(); }
    inline float getPeriodSec()     const {
        const float f = Freq_smoothed.get();
        return (f > 1e-9f) ? (1.0f / f) : 0.0f;
    }

    inline bool isReady() const { return A_var.isReady() && Freq_smoothed.isReady(); }
    inline bool isFreqReady() const { return Freq_smoothed.isReady(); }
    inline bool isVarReady()  const { return A_var.isReady(); }

    // Horizon of the σ_a EWMA, in seconds, as last used.  The variance is a
    // two-stage EWMA at this horizon (mean/square, then the variance itself),
    // so the effective memory is about twice it.
    inline float getVarianceHorizonSec() const { return tau_var_sec; }

    inline void setTauFreq(float t) {
        tau_freq = std::max(1e-3f, t);
        last_dt_freq = -1.0f;  // force recompute on next update
    }

    // σ_a averaging horizon, in periods of the tuning frequency.
    inline void setKPeriods(float k) {
        if (std::isfinite(k) && k > 0.0f) K_periods = k;
    }
    inline float getKPeriods() const { return K_periods; }

    // Absolute clamps on that horizon [s].  They exist so a degenerate tuning
    // frequency cannot make the estimator either a pass-through or a constant;
    // in the wave band the product K_periods * T_z should decide, not these.
    inline void setVarianceHorizonBounds(float min_s, float max_s) {
        if (!(std::isfinite(min_s) && std::isfinite(max_s))) return;
        if (!(min_s > 0.0f && max_s >= min_s)) return;
        tau_var_min_sec = min_s;
        tau_var_max_sec = max_s;
    }
    inline float getVarianceHorizonMinSec() const { return tau_var_min_sec; }
    inline float getVarianceHorizonMaxSec() const { return tau_var_max_sec; }

    // Clamps on the tuning frequency the horizon is derived from [Hz].  The
    // 0.05 Hz floor is a 20 s period, which is past the longest swell the wave
    // band admits.
    inline void setFrequencyBounds(float min_hz, float max_hz) {
        if (!(std::isfinite(min_hz) && std::isfinite(max_hz))) return;
        if (!(min_hz > 0.0f && max_hz >= min_hz)) return;
        f_min_hz = min_hz;
        f_max_hz = max_hz;
    }

private:
    // K_periods: dimensionless factor such that τ_var_dyn ≈ K_periods * T_eff
    float K_periods = 2.0f;
    float tau_var_min_sec = 0.3f;   // seconds: don't go *too* twitchy
    float tau_var_max_sec = 60.0f;  // seconds: don't be glacial
    float f_min_hz = 0.05f;         // 20 s period max
    float f_max_hz = 5.0f;          // avoid crazy high freq
    float tau_var_sec = 0.0f;       // last horizon used, for inspection
    float tau_freq  = 1.0f;   // seconds
    float last_dt_freq  = -1.0f;
    float alpha_freq = 0.0f;

    DebiasedEMA A_mean, A_sq, A_var;
    DebiasedEMA Freq_smoothed;

    inline void updateAlphaFreq(float dt_s) {
        if (dt_s == last_dt_freq) return;
        last_dt_freq = dt_s;
        alpha_freq = 1.0f - std::exp(-dt_s / tau_freq);
    }
};

#ifdef SEA_STATE_TUNER_TEST
#include <iostream>

static inline void SeaStateAutoTuner_test() {
    constexpr float Fs   = 200.0f;
    constexpr float DT   = 1.0f / Fs;
    constexpr float F_HZ = 0.5f;
    const float omega = 2.0f * 3.14159265358979323846f * F_HZ;

    // K_periods ≈ 1.5 → τ_var_dyn ≈ 1.5 * T = 3 s at 0.5 Hz
    // ~9 s (~4.5 periods) for 95% response
    SeaStateAutoTuner tuner(1.5f, 3.0f);

    float t = 0.0f;
    for (int n = 0; n < int(20.0f / DT); ++n) {
        t += DT;
        const float a = -std::sin(omega * t) * omega * omega;
        tuner.update(DT, a, F_HZ);
    }

    std::cerr << "[AutoTuner] σ_a=" << tuner.getAccelStd()
              << " m/s², f=" << tuner.getFrequencyHz()
              << " Hz\n";
}
#endif
