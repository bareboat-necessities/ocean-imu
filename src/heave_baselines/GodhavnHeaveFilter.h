#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Literature baseline: the adaptive heave filter of

    J. M. Godhavn, "Adaptive tuning of heave filter in motion sensor",
    OCEANS'98, pp. 174-178, 1998,

  in the form analysed by Richter, Schneider, Walser and Sawodny,
  "Real-time heave motion estimation using adaptive filtering techniques",
  IFAC World Congress 2014 (their Eqs. 3-12).

  Heave is the levelled vertical acceleration through

      H(s) = s^2 / (s^2 + 2 zeta w_c s + w_c^2)^2,   zeta = 1/sqrt(2),

  i.e. a double integrator in series with a fourth-order high pass.  The double
  zero at s = 0 rejects a constant acceleration bias; the price is phase lead
  in the wave band.  It is realised as two identical second-order sections

      x'' + 2 zeta w_c x' + w_c^2 x = u,   y = x',

  whose product is H, integrated with RK4 at the sample rate.  The states carry
  across retuning, so a cutoff change never restarts the filter.

  Adaptation (Godhavn: wave height, wave frequency, sensor noise and bias are
  estimated online and the cutoff follows the sea state).  Every
  adapt_period_s the cutoff target minimises the error budget

      J(w_c) = (A^2/2) |1 - w_p^4 / D(j w_p)^2|^2         sinusoid error
             + q_w / (2^(7/2) w_c^3)                      white accel noise
             + q_b / (2^(7/2) w_c^5)                      accel bias random walk

  with D(s) = s^2 + sqrt(2) w_c s + w_c^2.  The two noise integrals are exact
  for this H.  A and w_p come from running second moments of the filter's
  own heave and heave rate (A^2 = 2 E[p^2], w_p^2 = E[v^2]/E[p^2]); q_w is
  estimated from the first difference of the input; q_b is the sensor
  bias-instability figure, a datasheet input.  The applied cutoff relaxes
  toward the target in log w_c so retuning does not excite transients.  A
  slow mean of the input (Godhavn's bias estimate) is subtracted ahead of the
  filter so the internal DC state stays near zero when w_c moves.
*/

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <type_traits>

namespace heave_baselines {

template <typename T = double>
class GodhavnHeaveFilter {
    static_assert(std::is_floating_point<T>::value,
                  "GodhavnHeaveFilter<T>: T must be a floating-point type.");

public:
    struct Config {
        bool adaptive = true;
        T wc_init = T(0.10);        // rad/s, also the fixed cutoff
        T wc_min = T(0.02);         // rad/s
        T wc_max = T(0.50);         // rad/s
        T warmup_s = T(60);         // fixed cutoff until moments settle
        T adapt_period_s = T(5);
        T wc_relax_tau_s = T(60);   // log-domain relaxation toward target
        T stats_tau_s = T(120);     // second-moment averaging
        T bias_init_s = T(2);       // plain mean before the slow EMA starts
        T bias_tau_s = T(300);
        // Accelerometer bias random-walk intensity, (m/s^2)^2/s.
        T q_bias_rw = T(2.5e-7);
        T wp_min = T(0.2);          // rad/s
        T wp_max = T(6.0);          // rad/s
    };

    GodhavnHeaveFilter() : GodhavnHeaveFilter(Config{}) {}
    explicit GodhavnHeaveFilter(const Config& cfg) : cfg_(cfg) { reset(); }

    void reset() {
        s_.fill(T(0));
        wc_ = cfg_.wc_init;
        wc_target_ = cfg_.wc_init;
        t_ = T(0);
        since_adapt_ = T(0);
        bias_ = T(0);
        bias_n_ = 0;
        m_pp_ = m_vv_ = T(0);
        m_dd_ = T(0);
        prev_u_ = T(0);
        have_prev_ = false;
        heave_ = vel_ = T(0);
    }

    // a_up: levelled vertical specific force minus gravity, m/s^2, z up.
    T update(T a_up, T dt) {
        if (!(dt > T(0)) || !std::isfinite(a_up)) return heave_;
        t_ += dt;

        // Bias estimate: plain mean first, slow EMA afterwards.
        if (t_ <= cfg_.bias_init_s) {
            ++bias_n_;
            bias_ += (a_up - bias_) / T(bias_n_);
        } else {
            bias_ += (dt / cfg_.bias_tau_s) * (a_up - bias_);
        }
        const T u = a_up - bias_;

        step_(u, dt);

        // Sensor white-noise level from the first difference of the input.
        if (have_prev_) {
            const T d = a_up - prev_u_;
            const T k = std::min(T(1), dt / cfg_.stats_tau_s);
            m_dd_ += k * (d * d - m_dd_);
        }
        prev_u_ = a_up;
        have_prev_ = true;

        const T k = std::min(T(1), dt / cfg_.stats_tau_s);
        m_pp_ += k * (heave_ * heave_ - m_pp_);
        m_vv_ += k * (vel_ * vel_ - m_vv_);

        if (cfg_.adaptive && t_ > cfg_.warmup_s) {
            since_adapt_ += dt;
            if (since_adapt_ >= cfg_.adapt_period_s) {
                wc_target_ = optimalCutoff_(dt);
                since_adapt_ = T(0);
            }
            const T a = std::min(T(1), dt / cfg_.wc_relax_tau_s);
            wc_ = std::exp(std::log(wc_) + a * (std::log(wc_target_) - std::log(wc_)));
        }
        return heave_;
    }

    T heave() const { return heave_; }
    T heaveRate() const { return vel_; }
    T cutoff() const { return wc_; }
    T cutoffTarget() const { return wc_target_; }
    T biasEstimate() const { return bias_; }
    T waveFrequencyEstimate() const { return wavePeak_(); }
    T amplitudeEstimate() const { return std::sqrt(T(2) * m_pp_); }
    // Sensor white-noise intensity, (m/s^2)^2/Hz.
    T whiteNoiseIntensity(T dt) const { return T(0.5) * m_dd_ * dt; }

    // Error budget of the class doc; exposed for tests.
    static T errorBudget(T wc, T A, T wp, T q_w, T q_b) {
        const T wc2 = wc * wc;
        const T re = wc2 - wp * wp;
        const T im = std::sqrt(T(2)) * wc * wp;
        // D^2 = (re + j im)^2
        const T d2r = re * re - im * im;
        const T d2i = T(2) * re * im;
        const T den = d2r * d2r + d2i * d2i;
        const T wp4 = wp * wp * wp * wp;
        // E = 1 - wp^4 / D^2
        const T er = T(1) - wp4 * d2r / den;
        const T ei = wp4 * d2i / den;
        const T k = std::pow(T(2), T(3.5));
        return T(0.5) * A * A * (er * er + ei * ei)
            + q_w / (k * wc2 * wc)
            + q_b / (k * wc2 * wc2 * wc);
    }

private:
    T wavePeak_() const {
        if (!(m_pp_ > T(0))) return T(NAN);
        return std::clamp(std::sqrt(m_vv_ / m_pp_), cfg_.wp_min, cfg_.wp_max);
    }

    T optimalCutoff_(T dt) const {
        const T wp = wavePeak_();
        const T A = amplitudeEstimate();
        if (!std::isfinite(wp) || !(A > T(0))) return wc_target_;
        const T q_w = whiteNoiseIntensity(dt);
        constexpr int N = 96;
        const T lo = std::log(cfg_.wc_min);
        const T hi = std::log(cfg_.wc_max);
        T best_w = wc_target_;
        T best_j = std::numeric_limits<T>::infinity();
        for (int i = 0; i < N; ++i) {
            const T w = std::exp(lo + (hi - lo) * T(i) / T(N - 1));
            const T j = errorBudget(w, A, wp, q_w, cfg_.q_bias_rw);
            if (j < best_j) { best_j = j; best_w = w; }
        }
        return best_w;
    }

    // s_ = [x1, x1', x2, x2'].  Output heave = x2', heave rate = x2''.
    void deriv_(const std::array<T, 4>& s, T u, std::array<T, 4>& ds) const {
        const T a = std::sqrt(T(2)) * wc_;
        const T b = wc_ * wc_;
        ds[0] = s[1];
        ds[1] = u - a * s[1] - b * s[0];
        ds[2] = s[3];
        ds[3] = s[1] - a * s[3] - b * s[2];
    }

    void step_(T u, T dt) {
        std::array<T, 4> k1, k2, k3, k4, tmp;
        deriv_(s_, u, k1);
        for (int i = 0; i < 4; ++i) tmp[i] = s_[i] + T(0.5) * dt * k1[i];
        deriv_(tmp, u, k2);
        for (int i = 0; i < 4; ++i) tmp[i] = s_[i] + T(0.5) * dt * k2[i];
        deriv_(tmp, u, k3);
        for (int i = 0; i < 4; ++i) tmp[i] = s_[i] + dt * k3[i];
        deriv_(tmp, u, k4);
        for (int i = 0; i < 4; ++i)
            s_[i] += dt / T(6) * (k1[i] + T(2) * k2[i] + T(2) * k3[i] + k4[i]);

        heave_ = s_[3];
        vel_ = s_[1] - std::sqrt(T(2)) * wc_ * s_[3] - wc_ * wc_ * s_[2];
    }

    Config cfg_;
    std::array<T, 4> s_{};
    T wc_ = T(0.1), wc_target_ = T(0.1);
    T t_ = T(0), since_adapt_ = T(0);
    T bias_ = T(0);
    long bias_n_ = 0;
    T m_pp_ = T(0), m_vv_ = T(0), m_dd_ = T(0);
    T prev_u_ = T(0);
    bool have_prev_ = false;
    T heave_ = T(0), vel_ = T(0);
};

} // namespace heave_baselines
