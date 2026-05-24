#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy
*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

#include "pii_observer/VerticalPIIObserver.h"
#include "freq/FrequencyTrackerPolicy.h"

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

template<typename T = float, bool WithBias = true, TrackerType TT = TrackerType::PLL>
class AdaptiveVerticalPII {
    static_assert(std::is_floating_point<T>::value,
                  "AdaptiveVerticalPII<T>: T must be a floating-point type.");

public:
    using Observer = VerticalPIIObserver<T, WithBias>;
    using ObserverConfig = typename Observer::Config;
    using ObserverAdaptConfig = typename Observer::AdaptConfig;
    using ObserverSnapshot = typename Observer::Snapshot;

    using AccelFreqTracker = TrackerPolicy<TT>;
    using AccelFreqTrackerConfig = detail::tracker_config_t<AccelFreqTracker, T>;

    struct Config {
        ObserverConfig observer = [] {
            ObserverConfig cfg{};
            cfg.r = T(0.150);
            cfg.tau_a = T(0.68);
            cfg.tau_d = T(49.0);
            cfg.kb = T(2.5e-5);
            cfg.lambda_b = T(3.0e-3);
            cfg.bias_limit = T(0.12);
            cfg.a_f_limit = T(50.0);
            cfg.v_limit = T(50.0);
            cfg.p_limit = T(20.0);
            cfg.S_limit = T(200.0);
            cfg.d_limit = T(20.0);
            return cfg;
        }();

        ObserverAdaptConfig adaptation = [] {
            ObserverAdaptConfig cfg{};
            cfg.enabled = true;
            cfg.min_confidence = T(0.22);
            cfg.f_disp_ref_hz = T(0.12);
            cfg.sigma_a_ref = T(0.95);
            cfg.input_smooth_tau = T(4.5);
            cfg.param_smooth_tau = T(7.5);
            cfg.r_freq_exp = T(0.28);
            cfg.r_sigma_exp = T(0.02);
            cfg.tau_a_freq_exp = T(-0.40);
            cfg.tau_a_sigma_exp = T(-0.03);
            cfg.tau_d_freq_exp = T(-0.03);
            cfg.tau_d_sigma_exp = T(-0.01);
            cfg.kb_freq_exp = T(0.02);
            cfg.kb_sigma_exp = T(0.08);
            cfg.r_min = T(0.145);
            cfg.r_max = T(0.225);
            cfg.tau_a_min = T(0.50);
            cfg.tau_a_max = T(0.90);
            cfg.tau_d_min = T(44.0);
            cfg.tau_d_max = T(58.0);
            cfg.kb_min = T(5e-6);
            cfg.kb_max = T(6e-5);
            return cfg;
        }();

        AccelFreqTrackerConfig accel_freq_tracker =
            detail::make_default_tracker_config<AccelFreqTrackerConfig, T>();

        T sigma_mean_tau_s = T(20.0);
        T sigma_var_tau_s  = T(6.0);
        T sigma_floor      = T(1e-4);

        bool auto_schedule_from_accel_freq = true;
        T auto_schedule_period_s = T(0.50);

        bool force_enable_adaptation_when_auto_schedule = true;
        T default_external_confidence = T(1.0);

        T fallback_confidence_floor = T(0.52);
        T fallback_confidence_when_locked = T(0.82);

        T coarse_schedule_blend = T(0.48);
        T coarse_schedule_confidence_floor = T(0.62);

        struct EnvelopeConfig {
            bool enabled = true;

            // Same concept as OU ACC_NOISE_FLOOR_SIGMA_DEFAULT.
            T acc_noise_floor_sigma = T(0.12);

            // SeaStateAutoTuner-style variance horizon:
            // tau_var ~= K_periods / f.
            T tuner_K_periods = T(2.0);

            // SeaStateAutoTuner frequency EMA horizon.
            T tuner_tau_freq_s = T(1.0);

            // OU-style mapping:
            // tau = tau_coeff * 0.5 / f
            // sigma = sigma_coeff * sqrt(max(var - noise_floor^2, 0)).
            T tau_coeff = T(1.38);
            T sigma_coeff = T(0.90);

            // Smooth tau/sigma toward targets.
            T adapt_tau_s = T(1.80);

            T tau_min_s = T(0.02);
            T tau_max_s = T(3.00);
            T sigma_max = T(6.00);
        } envelope{};
    };

    struct Snapshot {
        ObserverSnapshot observer{};

        T accel_mean = T(0);
        T accel_var  = T(0);
        T accel_sigma = T(0);

        T last_sched_dt = T(0);
        T sched_accum_s = T(0);

        bool auto_schedule_from_accel_freq = true;
        bool observer_adaptation_enabled = false;

        T accel_freq_hz = T(0);
        T accel_freq_confidence_raw = T(0);
        T accel_freq_confidence_used = T(0);
        T accel_freq_sched_hz = T(0);
        T accel_freq_coarse_hz = T(0);

        bool accel_freq_locked = false;
        bool accel_freq_has_coarse = false;

        bool envelope_enabled = false;
        bool envelope_ready = false;
        bool adaptation_used_envelope = false;

        T envelope_freq_hz = T(0);
        T envelope_accel_var = T(0);

        T envelope_tau_target = T(0);
        T envelope_tau_applied = T(0);
        T envelope_sigma_target = T(0);
        T envelope_sigma_applied = T(0);

        T displacement_scale_m = T(0);
        T vertical_speed_envelope_mps = T(0);
    };

public:
    explicit AdaptiveVerticalPII(const Config& cfg = Config())
        : observer_(cfg.observer, cfg.adaptation), accel_freq_tracker_() {
        configure(cfg);
        reset();
    }

    void configure(const Config& cfg) {
        cfg_ = sanitizeConfig_(cfg);

        observer_.configure(cfg_.observer);
        observer_.configure_adaptation(cfg_.adaptation);
        detail::tracker_configure(accel_freq_tracker_, cfg_.accel_freq_tracker);

        if (cfg_.auto_schedule_from_accel_freq &&
            cfg_.force_enable_adaptation_when_auto_schedule) {
            observer_.set_adaptation_enabled(true);
        }

        resetSchedulerState_();
        resetEnvelopeState_();
    }

    void reset(T p0 = T(0),
               T v0 = T(0),
               T a_f0 = T(0),
               T S0 = T(0),
               T d0 = T(0),
               T b0 = T(0))
    {
        observer_.reset(p0, v0, a_f0, S0, d0, b0);

        accel_mean_  = T(0);
        accel_var_   = T(0);
        accel_sigma_ = cfg_.sigma_floor;

        detail::tracker_reset(
            accel_freq_tracker_,
            detail::tracker_init_frequency_hz(cfg_.accel_freq_tracker, T(0.12))
        );

        resetSchedulerState_();
        resetEnvelopeState_();
    }

    T update(T a_meas, T dt) {
        if (!(std::isfinite(dt) && dt > T(0))) {
            return observer_.displacement();
        }

        updateAccelSigma_(a_meas, dt);

        detail::tracker_step(accel_freq_tracker_, a_meas, dt);

        updateEnvelope_(a_meas, dt);

        if (cfg_.auto_schedule_from_accel_freq) {
            sched_accum_s_ += dt;

            if (sched_accum_s_ >= cfg_.auto_schedule_period_s) {
                const T dt_sched = sched_accum_s_;
                sched_accum_s_ = T(0);
                last_sched_dt_ = dt_sched;

                if (observer_.adaptation_enabled()) {
                    updateAdaptationFromAccelFrequencyProxy(dt_sched);
                }
            }
        }

        return observer_.update(a_meas, dt);
    }

    void updateAdaptationFromDisplacementFrequency(T f_disp_hz, T dt_est,
                                                   T confidence = std::numeric_limits<T>::quiet_NaN())
    {
        if (!std::isfinite(confidence)) {
            confidence = cfg_.default_external_confidence;
        }
        confidence = clamp01_(confidence);

        last_sched_dt_ = dt_est;
        last_auto_confidence_used_ = confidence;
        last_adaptation_used_envelope_ = false;

        observer_.update_adaptation(
            f_disp_hz,
            accel_sigma_,
            confidence,
            dt_est
        );
    }

    void updateAdaptationExternal(T f_disp_hz, T sigma_a, T dt_est,
                                  T confidence = std::numeric_limits<T>::quiet_NaN())
    {
        if (!std::isfinite(confidence)) {
            confidence = cfg_.default_external_confidence;
        }
        confidence = clamp01_(confidence);

        last_sched_dt_ = dt_est;
        last_auto_confidence_used_ = confidence;
        last_adaptation_used_envelope_ = false;

        observer_.update_adaptation(
            f_disp_hz,
            sigma_a,
            confidence,
            dt_est
        );
    }

    /*
      Auto adaptation source.

      Old behavior:
          f = accel-frequency tracker
          sigma = accel_sigma_ from local EWMA

      New behavior:
          when OU-style envelope tuner is ready:
              f = envelope_freq_hz_
              sigma = envelope_sigma_applied_
          otherwise fallback to old behavior.

      Important:
          The displacement envelope itself is NOT fed back.
          Only the envelope tuner's frequency and acceleration sigma are used
          as cleaner scheduling signals for VerticalPIIObserver.
    */
    void updateAdaptationFromAccelFrequencyProxy(T dt_est) {
        T f_sched_hz = computeScheduledFrequencyHz_();
        T sigma_for_adapt = accel_sigma_;
        bool used_envelope = false;

        if (cfg_.envelope.enabled && envelope_ready_) {
            if (std::isfinite(envelope_freq_hz_) && envelope_freq_hz_ > T(0)) {
                f_sched_hz = envelope_freq_hz_;
                used_envelope = true;
            }

            if (std::isfinite(envelope_sigma_applied_) && envelope_sigma_applied_ >= T(0)) {
                sigma_for_adapt = std::max(envelope_sigma_applied_, cfg_.sigma_floor);
                used_envelope = true;
            }
        }

        if (!(std::isfinite(f_sched_hz) && f_sched_hz > T(0))) {
            return;
        }

        const T conf_used = computeFallbackConfidence_();

        last_auto_confidence_used_ = conf_used;
        last_auto_tracker_locked_ = detail::tracker_is_locked(accel_freq_tracker_);
        last_auto_tracker_has_coarse_ = detail::tracker_has_coarse(accel_freq_tracker_);
        last_auto_sched_freq_hz_ = f_sched_hz;
        last_auto_coarse_freq_hz_ = last_auto_tracker_has_coarse_
            ? detail::tracker_get_coarse_frequency_hz<AccelFreqTracker, T>(accel_freq_tracker_)
            : T(0);

        last_adaptation_used_envelope_ = used_envelope;

        observer_.update_adaptation(
            f_sched_hz,
            sigma_for_adapt,
            conf_used,
            dt_est
        );
    }

    void setAutoScheduleFromAccelFreq(bool on) {
        cfg_.auto_schedule_from_accel_freq = on;
        if (on && cfg_.force_enable_adaptation_when_auto_schedule) {
            observer_.set_adaptation_enabled(true);
        }
    }

    void setAutoSchedulePeriod(T period_s) {
        if (std::isfinite(period_s) && period_s > T(0)) {
            cfg_.auto_schedule_period_s = period_s;
        }
    }

    void setFallbackConfidenceFloor(T c) {
        if (std::isfinite(c)) {
            cfg_.fallback_confidence_floor = clamp01_(c);
        }
    }

    void setFallbackConfidenceWhenLocked(T c) {
        if (std::isfinite(c)) {
            cfg_.fallback_confidence_when_locked = clamp01_(c);
        }
    }

    Observer& observer() { return observer_; }
    const Observer& observer() const { return observer_; }

    AccelFreqTracker& accelFreqTracker() { return accel_freq_tracker_; }
    const AccelFreqTracker& accelFreqTracker() const { return accel_freq_tracker_; }

    T displacement() const { return observer_.displacement(); }
    T velocity() const     { return observer_.velocity(); }
    T accelFiltered() const { return observer_.accel_filtered(); }

    T accelMean() const  { return accel_mean_; }
    T accelVar() const   { return accel_var_; }
    T accelSigma() const { return accel_sigma_; }

    T accelFrequencyHz() const {
        return detail::tracker_get_frequency_hz<AccelFreqTracker, T>(accel_freq_tracker_);
    }

    T accelFrequencyConfidence() const {
        return detail::tracker_get_confidence<AccelFreqTracker, T>(accel_freq_tracker_);
    }

    T lastSchedulerDt() const { return last_sched_dt_; }
    T autoScheduledConfidenceUsed() const { return last_auto_confidence_used_; }

    bool adaptationUsedEnvelope() const { return last_adaptation_used_envelope_; }

    bool envelopeReady() const { return envelope_ready_; }

    T envelopeFrequencyHz() const { return envelope_freq_hz_; }
    T envelopeAccelVariance() const { return envelope_accel_var_; }

    T envelopeTauTarget() const { return envelope_tau_target_; }
    T envelopeTauApplied() const { return envelope_tau_applied_; }

    T envelopeSigmaTarget() const { return envelope_sigma_target_; }
    T envelopeSigmaApplied() const { return envelope_sigma_applied_; }

    T displacementScale(bool smoothed = true) const {
        const T tau = smoothed ? envelope_tau_applied_ : envelope_tau_target_;
        const T sigma = smoothed ? envelope_sigma_applied_ : envelope_sigma_target_;

        if (!(std::isfinite(tau) && tau > T(0)) ||
            !(std::isfinite(sigma) && sigma >= T(0))) {
            return std::numeric_limits<T>::quiet_NaN();
        }

        const T pi = T(3.1415926535897932384626433832795L);
        const T C_HS = T(2) * std::sqrt(T(2)) / (pi * pi);

        return T(0.5) * C_HS * sigma * tau * tau;
    }

    T verticalSpeedEnvelope(bool smoothed = true) const {
        const T tau = smoothed ? envelope_tau_applied_ : envelope_tau_target_;
        const T sigma = smoothed ? envelope_sigma_applied_ : envelope_sigma_target_;

        if (!(std::isfinite(tau) && tau > T(0)) ||
            !(std::isfinite(sigma) && sigma >= T(0))) {
            return std::numeric_limits<T>::quiet_NaN();
        }

        const T pi = T(3.1415926535897932384626433832795L);
        const T K = std::sqrt(T(2)) / pi;

        return K * sigma * tau;
    }

    Snapshot snapshot() const {
        Snapshot s;
        s.observer = observer_.snapshot();

        s.accel_mean = accel_mean_;
        s.accel_var = accel_var_;
        s.accel_sigma = accel_sigma_;

        s.last_sched_dt = last_sched_dt_;
        s.sched_accum_s = sched_accum_s_;

        s.auto_schedule_from_accel_freq = cfg_.auto_schedule_from_accel_freq;
        s.observer_adaptation_enabled = observer_.adaptation_enabled();

        s.accel_freq_hz =
            detail::tracker_get_frequency_hz<AccelFreqTracker, T>(accel_freq_tracker_);
        s.accel_freq_confidence_raw =
            detail::tracker_get_confidence<AccelFreqTracker, T>(accel_freq_tracker_);
        s.accel_freq_confidence_used = last_auto_confidence_used_;
        s.accel_freq_sched_hz = last_auto_sched_freq_hz_;
        s.accel_freq_coarse_hz = last_auto_coarse_freq_hz_;

        s.accel_freq_locked = detail::tracker_is_locked(accel_freq_tracker_);
        s.accel_freq_has_coarse = detail::tracker_has_coarse(accel_freq_tracker_);

        s.envelope_enabled = cfg_.envelope.enabled;
        s.envelope_ready = envelope_ready_;
        s.adaptation_used_envelope = last_adaptation_used_envelope_;

        s.envelope_freq_hz = envelope_freq_hz_;
        s.envelope_accel_var = envelope_accel_var_;

        s.envelope_tau_target = envelope_tau_target_;
        s.envelope_tau_applied = envelope_tau_applied_;
        s.envelope_sigma_target = envelope_sigma_target_;
        s.envelope_sigma_applied = envelope_sigma_applied_;

        s.displacement_scale_m = displacementScale(true);
        s.vertical_speed_envelope_mps = verticalSpeedEnvelope(true);

        return s;
    }

private:
    struct DebiasedEMA_ {
        T value = T(0);
        T weight = T(0);

        void reset() {
            value = T(0);
            weight = T(0);
        }

        void update(T x, T alpha) {
            alpha = std::clamp(alpha, T(0), T(1));
            value  = (T(1) - alpha) * value + alpha * x;
            weight = (T(1) - alpha) * weight + alpha;
        }

        T get() const {
            return (weight > T(1e-12)) ? (value / weight) : T(0);
        }

        bool isReady() const {
            return weight > T(1e-6);
        }
    };

private:
    static constexpr T eps_() {
        return T(1e-9);
    }

    static T clamp01_(T x) {
        return std::clamp(x, T(0), T(1));
    }

    static T finite_or_default_(T x, T def) {
        return std::isfinite(x) ? x : def;
    }

    static T one_pole_alpha_(T dt, T tau) {
        if (!(std::isfinite(dt) && dt > T(0))) return T(0);
        if (!(std::isfinite(tau) && tau > T(0))) return T(1);
        const T a = dt / (tau + dt);
        return std::clamp(a, T(0), T(1));
    }

    static T exp_alpha_(T dt, T tau) {
        if (!(std::isfinite(dt) && dt > T(0))) return T(0);
        if (!(std::isfinite(tau) && tau > T(0))) return T(1);

        T a = T(1) - std::exp(-dt / tau);
        if (!std::isfinite(a)) a = T(1);

        return std::clamp(a, T(0), T(1));
    }

    static Config sanitizeConfig_(Config cfg) {
        cfg.sigma_mean_tau_s =
            std::max(finite_or_default_(cfg.sigma_mean_tau_s, T(20)), eps_());
        cfg.sigma_var_tau_s =
            std::max(finite_or_default_(cfg.sigma_var_tau_s, T(6)), eps_());
        cfg.sigma_floor =
            std::max(finite_or_default_(cfg.sigma_floor, T(1e-4)), T(0));

        cfg.auto_schedule_period_s =
            std::max(finite_or_default_(cfg.auto_schedule_period_s, T(0.50)), eps_());

        cfg.default_external_confidence =
            clamp01_(finite_or_default_(cfg.default_external_confidence, T(1)));
        cfg.fallback_confidence_floor =
            clamp01_(finite_or_default_(cfg.fallback_confidence_floor, T(0.52)));
        cfg.fallback_confidence_when_locked =
            clamp01_(finite_or_default_(cfg.fallback_confidence_when_locked, T(0.82)));
        cfg.coarse_schedule_blend =
            clamp01_(finite_or_default_(cfg.coarse_schedule_blend, T(0.48)));
        cfg.coarse_schedule_confidence_floor =
            clamp01_(finite_or_default_(cfg.coarse_schedule_confidence_floor, T(0.62)));

        cfg.envelope.acc_noise_floor_sigma =
            std::max(finite_or_default_(cfg.envelope.acc_noise_floor_sigma, T(0.12)), T(0));

        cfg.envelope.tuner_K_periods =
            std::max(finite_or_default_(cfg.envelope.tuner_K_periods, T(2.0)), eps_());

        cfg.envelope.tuner_tau_freq_s =
            std::max(finite_or_default_(cfg.envelope.tuner_tau_freq_s, T(1.0)), eps_());

        cfg.envelope.tau_coeff =
            std::max(finite_or_default_(cfg.envelope.tau_coeff, T(1.38)), eps_());

        cfg.envelope.sigma_coeff =
            std::max(finite_or_default_(cfg.envelope.sigma_coeff, T(0.90)), T(0));

        cfg.envelope.adapt_tau_s =
            std::max(finite_or_default_(cfg.envelope.adapt_tau_s, T(1.80)), eps_());

        cfg.envelope.tau_min_s =
            std::max(finite_or_default_(cfg.envelope.tau_min_s, T(0.02)), eps_());

        cfg.envelope.tau_max_s =
            std::max(finite_or_default_(cfg.envelope.tau_max_s, T(3.0)),
                     cfg.envelope.tau_min_s + eps_());

        cfg.envelope.sigma_max =
            std::max(finite_or_default_(cfg.envelope.sigma_max, T(6.0)), T(0));

        if (cfg.auto_schedule_from_accel_freq &&
            cfg.force_enable_adaptation_when_auto_schedule) {
            cfg.adaptation.enabled = true;
        }

        return cfg;
    }

    void resetSchedulerState_() {
        sched_accum_s_ = T(0);
        last_sched_dt_ = T(0);
        last_auto_confidence_used_ = T(0);
        last_auto_sched_freq_hz_ = T(0);
        last_auto_coarse_freq_hz_ = T(0);
        last_auto_tracker_locked_ = false;
        last_auto_tracker_has_coarse_ = false;
        last_adaptation_used_envelope_ = false;
    }

    void resetEnvelopeState_() {
        envelope_A_mean_.reset();
        envelope_A_sq_.reset();
        envelope_A_var_.reset();
        envelope_Freq_smoothed_.reset();

        envelope_ready_ = false;

        envelope_freq_hz_ = detail::tracker_init_frequency_hz(cfg_.accel_freq_tracker, T(0.12));
        envelope_accel_var_ = T(0);

        envelope_tau_target_ = std::clamp(
            cfg_.envelope.tau_coeff * T(0.5) / std::max(envelope_freq_hz_, eps_()),
            cfg_.envelope.tau_min_s,
            cfg_.envelope.tau_max_s);

        envelope_tau_applied_ = envelope_tau_target_;

        envelope_sigma_target_ = T(0);
        envelope_sigma_applied_ = T(0);
    }

    void updateAccelSigma_(T a_meas, T dt) {
        if (!std::isfinite(a_meas) || !(std::isfinite(dt) && dt > T(0))) {
            return;
        }

        const T alpha_mean = one_pole_alpha_(dt, cfg_.sigma_mean_tau_s);
        const T alpha_var  = one_pole_alpha_(dt, cfg_.sigma_var_tau_s);

        accel_mean_ += alpha_mean * (a_meas - accel_mean_);
        const T e = a_meas - accel_mean_;

        accel_var_ += alpha_var * (e * e - accel_var_);
        if (!(accel_var_ >= T(0)) || !std::isfinite(accel_var_)) {
            accel_var_ = T(0);
        }

        accel_sigma_ = std::sqrt(accel_var_);
        if (!(std::isfinite(accel_sigma_) && accel_sigma_ >= cfg_.sigma_floor)) {
            accel_sigma_ = cfg_.sigma_floor;
        }
    }

    void updateEnvelope_(T a_meas, T dt) {
        if (!cfg_.envelope.enabled) {
            envelope_ready_ = false;
            return;
        }

        if (!std::isfinite(a_meas) || !(std::isfinite(dt) && dt > T(0))) {
            return;
        }

        T f_input_hz = computeScheduledFrequencyHz_();

        if (!(std::isfinite(f_input_hz) && f_input_hz > T(0))) {
            f_input_hz = envelope_freq_hz_;
        }

        if (!(std::isfinite(f_input_hz) && f_input_hz > T(0))) {
            f_input_hz = T(0.12);
        }

        constexpr T F_MIN = T(0.05);
        constexpr T F_MAX = T(5.00);

        f_input_hz = std::clamp(f_input_hz, F_MIN, F_MAX);

        const T alpha_freq = exp_alpha_(dt, cfg_.envelope.tuner_tau_freq_s);
        envelope_Freq_smoothed_.update(f_input_hz, alpha_freq);

        T f_eff = envelope_Freq_smoothed_.isReady()
            ? envelope_Freq_smoothed_.get()
            : f_input_hz;

        if (!std::isfinite(f_eff)) {
            f_eff = f_input_hz;
        }

        f_eff = std::clamp(f_eff, F_MIN, F_MAX);
        envelope_freq_hz_ = f_eff;

        constexpr T TAU_VAR_MIN = T(0.30);
        constexpr T TAU_VAR_MAX = T(60.0);

        const T T_eff = T(1) / std::max(f_eff, eps_());
        const T tau_var_dyn =
            std::clamp(cfg_.envelope.tuner_K_periods * T_eff,
                       TAU_VAR_MIN,
                       TAU_VAR_MAX);

        const T alpha_var = exp_alpha_(dt, tau_var_dyn);

        envelope_A_mean_.update(a_meas, alpha_var);
        envelope_A_sq_.update(a_meas * a_meas, alpha_var);

        const T mu = envelope_A_mean_.get();
        const T var_inst = std::max(T(0), envelope_A_sq_.get() - mu * mu);

        envelope_A_var_.update(var_inst, alpha_var);

        envelope_accel_var_ = envelope_A_var_.get();
        if (!(std::isfinite(envelope_accel_var_) && envelope_accel_var_ >= T(0))) {
            envelope_accel_var_ = T(0);
        }

        const T noise_var =
            cfg_.envelope.acc_noise_floor_sigma *
            cfg_.envelope.acc_noise_floor_sigma;

        const T wave_var = std::max(T(0), envelope_accel_var_ - noise_var);

        envelope_sigma_target_ =
            std::sqrt(wave_var) * cfg_.envelope.sigma_coeff;

        envelope_sigma_target_ =
            std::min(std::max(envelope_sigma_target_, T(0)), cfg_.envelope.sigma_max);

        envelope_tau_target_ =
            cfg_.envelope.tau_coeff * T(0.5) / std::max(f_eff, eps_());

        envelope_tau_target_ =
            std::clamp(envelope_tau_target_,
                       cfg_.envelope.tau_min_s,
                       cfg_.envelope.tau_max_s);

        const T alpha_adapt = exp_alpha_(dt, cfg_.envelope.adapt_tau_s);

        if (!(std::isfinite(envelope_tau_applied_) && envelope_tau_applied_ > T(0))) {
            envelope_tau_applied_ = envelope_tau_target_;
        } else {
            envelope_tau_applied_ +=
                alpha_adapt * (envelope_tau_target_ - envelope_tau_applied_);
        }

        if (!(std::isfinite(envelope_sigma_applied_) && envelope_sigma_applied_ >= T(0))) {
            envelope_sigma_applied_ = envelope_sigma_target_;
        } else {
            envelope_sigma_applied_ +=
                alpha_adapt * (envelope_sigma_target_ - envelope_sigma_applied_);
        }

        envelope_ready_ =
            envelope_Freq_smoothed_.isReady() &&
            envelope_A_var_.isReady() &&
            std::isfinite(envelope_freq_hz_) &&
            std::isfinite(envelope_tau_applied_) &&
            std::isfinite(envelope_sigma_applied_);
    }

    T computeScheduledFrequencyHz_() const {
        T f_sched =
            detail::tracker_get_frequency_hz<AccelFreqTracker, T>(accel_freq_tracker_);

        if (!(std::isfinite(f_sched) && f_sched > T(0))) {
            f_sched =
                detail::tracker_get_raw_frequency_hz<AccelFreqTracker, T>(accel_freq_tracker_);
        }

        if (detail::tracker_has_coarse(accel_freq_tracker_)) {
            const T coarse_hz =
                detail::tracker_get_coarse_frequency_hz<AccelFreqTracker, T>(accel_freq_tracker_);

            if (std::isfinite(coarse_hz) && coarse_hz > T(0)) {
                T blend = cfg_.coarse_schedule_blend;

                if (!detail::tracker_is_locked(accel_freq_tracker_)) {
                    blend = std::max(blend, T(0.75));
                }

                if (!(std::isfinite(f_sched) && f_sched > T(0))) {
                    f_sched = coarse_hz;
                } else {
                    f_sched += blend * (coarse_hz - f_sched);
                }
            }
        }

        return f_sched;
    }

    T computeFallbackConfidence_() const {
        T conf =
            detail::tracker_get_confidence<AccelFreqTracker, T>(accel_freq_tracker_);

        if (!std::isfinite(conf)) {
            conf = T(0);
        }

        if (detail::tracker_has_coarse(accel_freq_tracker_)) {
            conf = std::max(
                conf,
                std::max(cfg_.fallback_confidence_floor,
                         cfg_.coarse_schedule_confidence_floor));
        }

        if (detail::tracker_is_locked(accel_freq_tracker_)) {
            conf = std::max(conf, cfg_.fallback_confidence_when_locked);
        }

        return clamp01_(conf);
    }

private:
    Config cfg_{};

    Observer observer_;
    AccelFreqTracker accel_freq_tracker_;

    T accel_mean_  = T(0);
    T accel_var_   = T(0);
    T accel_sigma_ = T(1e-4);

    T sched_accum_s_ = T(0);
    T last_sched_dt_ = T(0);

    T last_auto_confidence_used_ = T(0);
    T last_auto_sched_freq_hz_ = T(0);
    T last_auto_coarse_freq_hz_ = T(0);
    bool last_auto_tracker_locked_ = false;
    bool last_auto_tracker_has_coarse_ = false;
    bool last_adaptation_used_envelope_ = false;

    DebiasedEMA_ envelope_A_mean_{};
    DebiasedEMA_ envelope_A_sq_{};
    DebiasedEMA_ envelope_A_var_{};
    DebiasedEMA_ envelope_Freq_smoothed_{};

    bool envelope_ready_ = false;

    T envelope_freq_hz_ = T(0.12);
    T envelope_accel_var_ = T(0);

    T envelope_tau_target_ = T(0);
    T envelope_tau_applied_ = T(0);

    T envelope_sigma_target_ = T(0);
    T envelope_sigma_applied_ = T(0);
};

template<typename T = float, bool WithBias = true>
using AdaptiveVerticalPII_PLL = AdaptiveVerticalPII<T, WithBias, TrackerType::PLL>;

template<typename T = float, bool WithBias = true>
using AdaptiveVerticalPII_Aranovskiy = AdaptiveVerticalPII<T, WithBias, TrackerType::ARANOVSKIY>;

template<typename T = float, bool WithBias = true>
using AdaptiveVerticalPII_KalmANF = AdaptiveVerticalPII<T, WithBias, TrackerType::KALMANF>;

template<typename T = float, bool WithBias = true>
using AdaptiveVerticalPII_ZeroCross = AdaptiveVerticalPII<T, WithBias, TrackerType::ZEROCROSS>;

} // namespace marine_obs
