/* Copyright (c) 2026 Mikhail Grushinskiy */
#define EIGEN_NON_ARDUINO
#include "detrend/AdaptiveWaveDetrender3D.h"
#include "detrend/AdaptiveWaveDetrender.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
using D = AdaptiveWaveDetrender3D;
using V = Eigen::Vector3f;
void require(bool ok, const std::string& what) {
    if (!ok) throw std::runtime_error(what);
}
D::Config config(float period) {
    D::Config c;
    c.init_wave_freq_hz = 1.0f / period;
    c.freq_smooth_tau_s = 12.0f;
    c.startup_hold_s = 0.0f;
    c.freq_learn_axis = 2;
    // Fixed threshold and an effectively transparent derivative LP make the
    // public-input crossing fixture's measured periods deterministic.
    c.slope_lpf_tau_s = 1e-6f;
    c.min_slope_threshold_abs = 0.25f;
    c.max_slope_threshold_abs = 0.25f;
    return c;
}
double normalized_error(float period, float initial, float target) {
    return std::log(double(period) / target) / std::log(double(initial) / target);
}
void external_reciprocal_steps() {
    D longer(config(4.0f)), shorter(config(16.0f));
    longer.reset(V::Zero()); shorter.reset(V::Zero());
    constexpr float dt = 0.005f;
    double max_symmetry_error = 0.0, max_law_error = 0.0;
    double t90_longer = -1, t90_shorter = -1;
    for (int k = 1; k <= 12000; ++k) {
        auto a = longer.update(V::Zero(), dt, 1.0f / 16.0f, true);
        auto b = shorter.update(V::Zero(), dt, 1.0f / 4.0f, true);
        const double ea = normalized_error(a.wave_period_s, 4, 16);
        const double eb = normalized_error(b.wave_period_s, 16, 4);
        const double elapsed = double(k) * dt;
        const double expected = std::exp(-elapsed / 12.0);
        max_symmetry_error = std::max(max_symmetry_error, std::abs(ea-eb));
        max_law_error = std::max(max_law_error, std::max(std::abs(ea-expected), std::abs(eb-expected)));
        require(a.freq_valid && b.freq_valid, "external frequency freshness lost");
        require(a.wave_clean.isZero() && b.wave_clean.isZero(), "frequency guidance invented heave");
        if (ea <= 0.1 && t90_longer < 0) t90_longer = elapsed;
        if (eb <= 0.1 && t90_shorter < 0) t90_shorter = elapsed;
    }
    require(max_symmetry_error < 0.002, "reciprocal external steps are not symmetric in log space");
    require(max_law_error < 0.002, "external EMA does not follow exp(-t/tau) in log space");
    require(std::abs(t90_longer-t90_shorter) < 0.08, "long period makes the return disproportionately slow");
    std::cout << "external 4->16/16->4: max_log_symmetry_error=" << max_symmetry_error
              << " max_log_law_error=" << max_law_error << " t90_s=" << t90_longer << ',' << t90_shorter << '\n';
}

// The learning channel is x - previous_baseline. Feeding the public baseline
// back ONLY in this test lets us specify a triangular learning signal without
// depending on the separate adaptive baseline's phase lag. Production internals
// remain private; both the derivative and Schmitt detection still execute.
struct CrossingInput {
    D d;
    float wave = 0.0f;
    bool first = true;
    explicit CrossingInput(float initial) : d(config(initial)) { d.reset(V::Zero()); }
    D::Output step(int index, int samples_per_period) {
        constexpr float dt = 0.0625f; // exactly represented timestamps/crossings
        const float slope = index % samples_per_period < samples_per_period/2 ? 1.0f : -1.0f;
        // 0 -> 0.4 crosses +0.25 at 5/8 dt, exactly the same
        // fraction as every later -1 -> +1 crossing. Without this seed,
        // the first measured period includes a known partial-sample offset.
        wave += dt * (first ? 0.4f : slope);
        first = false;
        V input = d.currentBaselineSlow(); input.z() += wave;
        return d.update(input, dt);
    }
};
void internal_measurement_law(float initial, float target) {
    CrossingInput input(initial);
    constexpr float dt = 0.0625f;
    const int period_samples = static_cast<int>(target/dt);
    double measurement_time = 0.0, max_error = 0.0;
    int updates = 0;
    float previous = input.d.currentWaveFreqHz();
    for (int k=0; k < 8*period_samples; ++k) {
        const auto out = input.step(k, period_samples);
        if (out.wave_freq_hz != previous) {
            // At constant forced period every accepted same-polarity crossing
            // supplies target seconds of EMA weight (not the old held period).
            measurement_time += target;
            ++updates;
            const double expected = std::exp(-measurement_time/12.0);
            const double got = normalized_error(out.wave_period_s, initial, target);
            max_error = std::max(max_error, std::abs(got-expected));
            previous = out.wave_freq_hz;
        }
    }
    require(updates >= 6, "internal fixture did not exercise enough accepted periods");
    require(input.d.frequencyValid(), "internal accepted periods are not fresh");
    require(max_error < 0.0005, "internal frequency does not use the measured-interval log EMA");
    std::cout << "internal " << initial << "->" << target << ": accepted=" << updates
              << " max_log_law_error=" << max_error << '\n';
}
double learned_period_return(float initial, float target) {
    CrossingInput input(initial);
    constexpr float dt = 0.0625f;
    const int old_samples=static_cast<int>(initial/dt), new_samples=static_cast<int>(target/dt);
    for(int k=0;k<12*old_samples;++k) input.step(k,old_samples);
    require(input.d.frequencyValid(), "initial period was not actually learned");
    require(std::abs(input.d.currentWavePeriodS()/initial-1) < 0.001f, "initial period did not settle");
    double t90=-1;
    for(int k=0;k<1600;++k) {
        const auto out=input.step(k,new_samples);
        if (std::abs(normalized_error(out.wave_period_s,initial,target)) <= 0.1 && t90<0)
            t90=(k+1)*double(dt);
    }
    require(t90>0 && t90<50, "learned-period transition did not settle promptly");
    return t90;
}
// The scalar fixture must not be the only guard of the 1D frequency law:
// check that implementation independently against the closed-form log EMA.
AdaptiveWaveDetrender::Config scalar_config(float period) {
    AdaptiveWaveDetrender::Config c;
    c.init_wave_freq_hz = 1.0f / period;
    c.freq_smooth_tau_s = 12.0f;
    c.startup_hold_s = 0.0f;
    c.slope_lpf_tau_s = 1e-6f;
    c.min_slope_threshold_abs = c.max_slope_threshold_abs = 0.25f;
    return c;
}
void scalar_frequency_laws() {
    AdaptiveWaveDetrender longer(scalar_config(4)), shorter(scalar_config(16));
    longer.reset(0); shorter.reset(0);
    constexpr float dt = 0.0625f;
    double max_external_error = 0;
    for (int k = 1; k <= 960; ++k) {
        const auto a = longer.update(0, dt, 1.0f / 16, true);
        const auto b = shorter.update(0, dt, 1.0f / 4, true);
        const double expected = std::exp(-double(k) * dt / 12.0);
        max_external_error = std::max(max_external_error,
            std::max(std::abs(normalized_error(a.wave_period_s, 4, 16) - expected),
                     std::abs(normalized_error(b.wave_period_s, 16, 4) - expected)));
        require(a.wave_clean == 0 && b.wave_clean == 0, "scalar guidance invented heave");
        if (k >= 2) require(a.freq_valid && b.freq_valid, "scalar guidance is not fresh");
    }
    require(max_external_error < 0.002, "scalar external frequency violates the log EMA");
    for (float initial : {4.0f, 16.0f}) {
        const float target = initial == 4.0f ? 16.0f : 4.0f;
        AdaptiveWaveDetrender d(scalar_config(initial));
        d.reset(0);
        const int samples = static_cast<int>(target / dt);
        float input = 0, previous = d.currentWaveFreqHz();
        double measurement_time = 0, max_error = 0;
        int accepted = 0;
        for (int k = 0; k < 8 * samples; ++k) {
            const float slope = k % samples < samples / 2 ? 1.0f : -1.0f;
            // 1D learns from raw input slope, not the 3D baseline residual.
            input += dt * (k == 0 ? 0.4f : slope);
            const auto out = d.update(input, dt);
            if (out.wave_freq_hz != previous) {
                measurement_time += target;
                ++accepted;
                max_error = std::max(max_error,
                    std::abs(normalized_error(out.wave_period_s, initial, target)
                             - std::exp(-measurement_time / 12.0)));
                previous = out.wave_freq_hz;
            }
        }
        require(accepted >= 6 && d.frequencyValid(), "scalar accepted-period fixture is inactive");
        require(max_error < 0.0005, "scalar internal frequency violates measured-interval log EMA");
        std::cout << "scalar internal " << initial << "->" << target
                  << ": accepted=" << accepted << " max_log_law_error=" << max_error << '\n';
    }
    std::cout << "scalar external max_log_law_error=" << max_external_error << '\n';
}
void invalid_guidance() {
    D d(config(8)); d.reset(V::Zero());
    for(float bad : {0.f,-1.f,NAN,INFINITY}) {
        d.update(V::Zero(),.005f,bad,true);
        require(d.currentWavePeriodS()==8.f && !d.frequencyValid(), "invalid external guidance contaminated frequency");
    }
    d.update(V::Zero(),.005f,.25f,false);
    require(d.currentWavePeriodS()==8.f, "invalid-flag guidance was applied");
}
}
int main() {
    try {
        scalar_frequency_laws();
        external_reciprocal_steps();
        internal_measurement_law(4,16);
        internal_measurement_law(16,4);
        const double longer=learned_period_return(4,16), shorter=learned_period_return(16,4);
        require(shorter<=longer+0.125, "a previously learned long period delays the shorter-period return");
        std::cout << "learned internal 4->16/16->4 t90_s=" << longer << ',' << shorter << '\n';
        invalid_guidance();
        std::cout << "detrend-log-frequency-test: all checks passed\n";
        return 0;
    } catch(const std::exception& e) { std::cerr << "FAIL: " << e.what() << '\n'; return 1; }
}
