#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <utility>

namespace ocean_imu::kalman::filter_detail {
// IntervalRate includes a rotation-equivalent delta-angle / dt. It is NOT an
// instantaneous endpoint sample. Existing time_update callers retain this mode.
enum class GyroSampleMode { IntervalRate, Instantaneous };

template<class Vector> class GyroInterval {
public:
    using T = std::remove_cv_t<std::remove_reference_t<decltype(std::declval<const Vector&>()[0])>>;
    struct Step { Vector rate; T noise_scale; bool predict; };
    void set_mode(GyroSampleMode m) { mode_ = m; }
    void reset() { have_ = false; }
    // Fixed sensor exposure/ODR, distinct from host delivery interval. Zero
    // selects Var(n_k)=Q_g/dt_k; a positive period selects Q_g/period.
    void set_sample_period(T period) { sample_period_ = std::isfinite(period) && period > T(0) ? period : T(0); }
    void prime(const Vector& rate, T dt) {
        if (finite_(rate) && std::isfinite(dt) && dt > T(0)) {
            previous_ = rate; previous_dt_ = dt; have_ = true;
        }
    }
    Step prepare(const Vector& rate, T dt) {
        Step out{rate, T(1), true};
        if (mode_ == GyroSampleMode::Instantaneous) {
            out.predict = have_;
            for (int i = 0; i < 3; ++i)
                out.rate[i] = have_ ? T(.5) * (previous_[i] + rate[i]) : T(0);
            // Do not halve process noise. Adjacent intervals share a sample.
            // For independent endpoint noise Var(n_k)=Q_g/dt_k, Young's
            // inequality bounds the covariance of the entire linearized
            // history by assigning twice each endpoint's contribution:
            // .5 dt^2 (Q_g/dt_prev + Q_g/dt). Keep at least the existing Q_g.
            // This is a conservative correlated-noise bound, not an exact
            // independent-white-noise likelihood. See docs/live-bias-consistency.md.
            if (have_) out.noise_scale = sample_period_ > T(0)
                ? std::max(T(1), dt / sample_period_)
                : std::max(T(1), T(.5) * (T(1) + dt / previous_dt_));
        }
        prime(rate, dt);
        return out;
    }
private:
    GyroSampleMode mode_ = GyroSampleMode::IntervalRate;
    static bool finite_(const Vector& v) {
        return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
    }
    static Vector zero_() { Vector v; for (int i=0; i<3; ++i) v[i]=T(0); return v; }
    Vector previous_ = zero_();
    T previous_dt_ = T(0), sample_period_ = T(0);
    bool have_ = false;
};
} // namespace filter_detail
