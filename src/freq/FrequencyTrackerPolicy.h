#pragma once

/*
  Copyright 2025-2026, Mikhail Grushinskiy
*/

#include <cmath>
#include <type_traits>

#include "util/WaveFilesSupport.h"
#include "freq/AranovskiyFreqTracker.h"
#include "freq/KalmANFFreqTracker.h"
#include "freq/PLLFreqTracker.h"
#include "freq/SchmittTriggerZCFreqTracker.h"

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

#ifndef ZERO_CROSSINGS_HYSTERESIS
#define ZERO_CROSSINGS_HYSTERESIS 0.04f
#endif

#ifndef ZERO_CROSSINGS_PERIODS
#define ZERO_CROSSINGS_PERIODS 1
#endif

// Tracker policy traits
template<TrackerType>
struct TrackerPolicy; // primary template (undefined)

// Aranovskiy
template<>
struct TrackerPolicy<TrackerType::ARANOVSKIY> {
    using Tracker = AranovskiyFreqTracker<double>;
    Tracker t;

    TrackerPolicy() : t() {
        const double omega_up   = (FREQ_GUESS * 2.0) * (2.0 * M_PI);
        const double k_gain     = 20.0;
        const double x1_0       = 0.0;
        const double omega_init = (FREQ_GUESS / 1.5) * 2.0 * M_PI;
        const double theta_0    = -(omega_init * omega_init);
        const double sigma_0    = theta_0;

        t.setParams(omega_up, k_gain);
        t.setState(x1_0, theta_0, sigma_0);
    }

    double run(float a, float dt) {
        t.update(static_cast<double>(a) / static_cast<double>(g_std), static_cast<double>(dt));
        return getFrequencyHz();
    }

    double getFrequencyHz() const { return t.getFrequencyHz(); }
    double getRawFrequencyHz() const { return t.getRawFrequencyHz(); }
    double getConfidence() const { return t.getConfidence(); }
    bool isLocked() const { return t.isLocked(); }
    bool hasCoarseEstimate() const { return t.hasCoarseEstimate(); }
    double getCoarseFrequencyHz() const { return t.getCoarseFrequencyHz(); }
};

// KalmANF
template<>
struct TrackerPolicy<TrackerType::KALMANF> {
    using Tracker = KalmANFFreqTracker<double>;
    Tracker t{};

    double run(float a, float dt) {
        double e = 0.0;
        t.process(static_cast<double>(a) / static_cast<double>(g_std),
                  static_cast<double>(dt), &e);
        return getFrequencyHz();
    }

    double getFrequencyHz() const { return t.getFrequencyHz(); }
    double getRawFrequencyHz() const { return t.getRawFrequencyHz(); }
    double getConfidence() const { return t.getConfidence(); }
    bool isLocked() const { return t.isLocked(); }
    bool hasCoarseEstimate() const { return t.hasCoarseEstimate(); }
    double getCoarseFrequencyHz() const { return t.getCoarseFrequencyHz(); }
};

// PLL
template<>
struct TrackerPolicy<TrackerType::PLL> {
    using Tracker = PLLFreqTracker<double>;
    using Config = typename Tracker::Config;

    Tracker t{};

    void configure(const Config& cfg) {
        t.configure(cfg);
    }

    void reset(double f_init_hz) {
        t.reset(f_init_hz);
    }

    void update(float a, float dt) {
        t.update(static_cast<double>(a) / static_cast<double>(g_std), static_cast<double>(dt));
    }

    double run(float a, float dt) {
        update(a, dt);
        return getFrequencyHz();
    }

    double getFrequencyHz() const { return t.getFrequencyHz(); }
    double getRawFrequencyHz() const { return t.getRawFrequencyHz(); }
    double getConfidence() const { return t.getConfidence(); }
    bool isLocked() const { return t.isLocked(); }
    bool hasCoarseEstimate() const { return t.hasCoarseEstimate(); }
    double getCoarseFrequencyHz() const { return t.getCoarseFrequencyHz(); }
};

// ZeroCross
template<>
struct TrackerPolicy<TrackerType::ZEROCROSS> {
    using Tracker = SchmittTriggerZCFreqTracker;
    Tracker t{ZERO_CROSSINGS_HYSTERESIS, ZERO_CROSSINGS_PERIODS};

    double run(float a, float dt) {
        t.update(a / g_std,
                 ZERO_CROSSINGS_SCALE /* max g */,
                 ZERO_CROSSINGS_DEBOUNCE_TIME,
                 ZERO_CROSSINGS_STEEPNESS_TIME, dt);
        return getFrequencyHz();
    }

    double getFrequencyHz() const {
        const float raw = t.getFrequencyHz();
        return isZeroCrossFallback_(raw) ? static_cast<double>(FREQ_GUESS) : static_cast<double>(raw);
    }
    double getRawFrequencyHz() const { return static_cast<double>(t.getRawFrequencyHz()); }
    double getConfidence() const { return static_cast<double>(t.getConfidence()); }
    bool isLocked() const { return t.isLocked(); }
    bool hasCoarseEstimate() const { return t.hasCoarseEstimate(); }
    double getCoarseFrequencyHz() const { return static_cast<double>(t.getCoarseFrequencyHz()); }

private:
    static bool isZeroCrossFallback_(float f) {
        return (f == SCHMITT_TRIGGER_FREQ_INIT || f == SCHMITT_TRIGGER_FALLBACK_FREQ);
    }
};

namespace marine_obs {
namespace detail {

template<typename T>
struct EmptyAccelFreqTrackerConfig {
    T f_init_hz = T(0.12);
};

template<typename Tracker, typename T, typename = void>
struct tracker_config_type {
    using type = EmptyAccelFreqTrackerConfig<T>;
};

template<typename Tracker, typename T>
struct tracker_config_type<Tracker, T, std::void_t<typename Tracker::Config>> {
    using type = typename Tracker::Config;
};

template<typename Tracker, typename T>
using tracker_config_t = typename tracker_config_type<Tracker, T>::type;

template<typename>
inline constexpr bool always_false_v = false;

template<typename Config, typename T>
Config make_default_tracker_config() {
    Config cfg{};

    if constexpr (requires { cfg.f_min_hz; })               cfg.f_min_hz = T(0.045);
    if constexpr (requires { cfg.f_max_hz; })               cfg.f_max_hz = T(0.35);
    if constexpr (requires { cfg.f_init_hz; })              cfg.f_init_hz = T(0.12);

    if constexpr (requires { cfg.pre_hp_hz; })              cfg.pre_hp_hz = T(0.015);
    if constexpr (requires { cfg.pre_lp_hz; })              cfg.pre_lp_hz = T(0.45);
    if constexpr (requires { cfg.demod_lp_hz; })            cfg.demod_lp_hz = T(0.05);
    if constexpr (requires { cfg.loop_bandwidth_hz; })      cfg.loop_bandwidth_hz = T(0.018);
    if constexpr (requires { cfg.loop_damping; })           cfg.loop_damping = T(1.0);
    if constexpr (requires { cfg.max_dfdt_hz_per_s; })      cfg.max_dfdt_hz_per_s = T(0.04);
    if constexpr (requires { cfg.recenter_tau_s; })         cfg.recenter_tau_s = T(12.0);
    if constexpr (requires { cfg.output_smooth_tau_s; })    cfg.output_smooth_tau_s = T(4.0);
    if constexpr (requires { cfg.power_tau_s; })            cfg.power_tau_s = T(14.0);
    if constexpr (requires { cfg.confidence_tau_s; })       cfg.confidence_tau_s = T(10.0);
    if constexpr (requires { cfg.lock_rms_min; })           cfg.lock_rms_min = T(0.012);
    if constexpr (requires { cfg.enable_coarse_assist; })   cfg.enable_coarse_assist = true;
    if constexpr (requires { cfg.coarse_hysteresis_frac; }) cfg.coarse_hysteresis_frac = T(0.20);
    if constexpr (requires { cfg.coarse_smooth_tau_s; })    cfg.coarse_smooth_tau_s = T(4.5);
    if constexpr (requires { cfg.coarse_pull_tau_s; })      cfg.coarse_pull_tau_s = T(3.5);
    if constexpr (requires { cfg.coarse_timeout_s; })       cfg.coarse_timeout_s = T(18.0);

    return cfg;
}

template<typename Config, typename T>
T tracker_init_frequency_hz(const Config& cfg, T fallback = T(0.12)) {
    if constexpr (requires { cfg.f_init_hz; }) {
        return static_cast<T>(cfg.f_init_hz);
    } else {
        return fallback;
    }
}

template<typename Tracker, typename Config>
void tracker_configure(Tracker& tracker, const Config& cfg) {
    if constexpr (requires { tracker.configure(cfg); }) {
        tracker.configure(cfg);
    } else {
        (void)tracker;
        (void)cfg;
    }
}

template<typename Tracker, typename T>
void tracker_reset(Tracker& tracker, T f_init_hz) {
    if constexpr (requires { tracker.reset(f_init_hz); }) {
        tracker.reset(f_init_hz);
    } else if constexpr (std::is_default_constructible_v<Tracker>) {
        (void)f_init_hz;
        tracker = Tracker{};
    } else {
        static_assert(always_false_v<Tracker>,
                      "Tracker must provide reset(f_init_hz) or be default-constructible.");
    }
}

template<typename Tracker, typename T>
void tracker_step(Tracker& tracker, T a_meas, T dt) {
    if constexpr (requires { tracker.update(a_meas, dt); }) {
        tracker.update(a_meas, dt);
    } else if constexpr (requires {
        tracker.run(static_cast<float>(a_meas), static_cast<float>(dt));
    }) {
        (void)tracker.run(static_cast<float>(a_meas), static_cast<float>(dt));
    } else {
        static_assert(always_false_v<Tracker>,
                      "Tracker must provide update(a, dt) or run(a, dt).");
    }
}

template<typename Tracker, typename T>
T tracker_get_frequency_hz(const Tracker& tracker) {
    return static_cast<T>(tracker.getFrequencyHz());
}

template<typename Tracker, typename T>
T tracker_get_raw_frequency_hz(const Tracker& tracker) {
    return static_cast<T>(tracker.getRawFrequencyHz());
}

template<typename Tracker, typename T>
T tracker_get_confidence(const Tracker& tracker) {
    return static_cast<T>(tracker.getConfidence());
}

template<typename Tracker>
bool tracker_is_locked(const Tracker& tracker) {
    return tracker.isLocked();
}

template<typename Tracker>
bool tracker_has_coarse(const Tracker& tracker) {
    return tracker.hasCoarseEstimate();
}

template<typename Tracker, typename T>
T tracker_get_coarse_frequency_hz(const Tracker& tracker) {
    return static_cast<T>(tracker.getCoarseFrequencyHz());
}

} // namespace detail
} // namespace marine_obs
