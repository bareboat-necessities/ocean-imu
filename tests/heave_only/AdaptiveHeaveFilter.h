#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Shared engine of the adaptive double-integrating heave filters
  (GodhavnHeaveFilter.h, RichterHeaveFilter.h).  Equation numbers are those of

    M. Richter, K. Schneider, D. Walser, O. Sawodny, "Real-time heave motion
    estimation using adaptive filtering techniques", 19th IFAC World Congress,
    pp. 10119-10125, 2014.

  Both filters share the denominator D(s)^2 = (s^2 + 2 zeta w_c s + w_c^2)^2,
  zeta = 1/sqrt(2).  It is realised as two identical sections
  x'' + 2 zeta w_c x' + w_c^2 x = u, y = x', integrated with RK4, so the
  states carry across retuning.  The second section's position state is
  s/D^2 applied to the input and its rate is s^2/D^2; the Law combines them
  into the filter's output.

  Adaptation (Section 2.4, Fig. 4): every ident_period_s the FFT of the
  buffered acceleration gives the dominant heave frequency w_p (peak of the
  heave spectrum) and the buffer gives the mean wave height (Eq. 14)

      A_p = sqrt( 2/N sum_j a_j^2 / w_p^4 ).

  The Law turns (sigma_n^2, w_p, A_p) into the cutoff.  The determined cutoff
  is strongly low-pass filtered (here: relaxation of log w_c, param_tau_s).

  Two extensions, both off by default and not in the paper:
    bias_instability_term  adds the variance that an accelerometer bias random
                           walk of intensity q_b leaves after the filter
                           (Law::errorBound, exact for each structure) and
                           minimises the extended bound numerically;
    subtract_input_mean    removes a slow running mean of the input ahead of
                           the filter so a retune does not move its DC state.

  A Law provides:
    static T optimalCutoff(T noise_density, T wp, T Ap);      // published
    static T errorBound(T wc, T wp, T Ap, T noise_density, T bias_rw_density);
    static void output(T wc, T wp, const T s[4], T& heave, T& rate);
*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

#include "AccelSpectrum.h"

namespace heave_only {

template <typename T, typename Law>
class AdaptiveHeaveFilter {
    static_assert(std::is_floating_point<T>::value,
                  "AdaptiveHeaveFilter<T>: T must be a floating-point type.");

public:
    using LawType = Law;

    struct Config {
        bool adaptive = true;
        T wc_init = T(0.10);          // rad/s; the fixed cutoff when !adaptive
        T wc_min = T(0.005);          // rad/s
        T wc_max = T(1.0);            // rad/s
        // Accelerometer white-noise spectral density sigma_n^2, (m/s^2)^2/Hz
        // (= sample variance * sample period), from the datasheet (Eq. 8).
        T noise_density = T(1.1e-6);

        // Identification (Fig. 4).
        T fs_ident_hz = T(4);
        int fft_size = 512;           // 128 s at 4 Hz
        T ident_period_s = T(30);
        T f_min_hz = T(0.04);
        T f_max_hz = T(1.0);
        T param_tau_s = T(60);

        // Extensions (see the header comment).
        bool bias_instability_term = false;
        T bias_rw_density = T(2.5e-7);  // q_b, (m/s^2)^2/s
        bool subtract_input_mean = false;
        T input_mean_tau_s = T(300);
    };

    AdaptiveHeaveFilter() : AdaptiveHeaveFilter(Config{}) {}
    explicit AdaptiveHeaveFilter(const Config& cfg)
        : cfg_(cfg), spectrum_(cfg.fs_ident_hz, cfg.fft_size) { reset(); }

    void reset() {
        s_[0] = s_[1] = s_[2] = s_[3] = T(0);
        wc_ = wc_target_ = cfg_.wc_init;
        wp_ = Ap_ = T(NAN);
        since_ident_ = T(0);
        mean_ = T(0);
        mean_n_ = 0;
        heave_ = vel_ = T(0);
        spectrum_ = AccelSpectrum<T>(cfg_.fs_ident_hz, cfg_.fft_size);
    }

    // a_up: levelled vertical specific force minus gravity, m/s^2, z up.
    T update(T a_up, T dt) {
        if (!(dt > T(0)) || !std::isfinite(a_up)) return heave_;

        T u = a_up;
        if (cfg_.subtract_input_mean) {
            ++mean_n_;
            const T k = std::max(T(1) / T(mean_n_), std::min(T(1), dt / cfg_.input_mean_tau_s));
            mean_ += k * (a_up - mean_);
            u -= mean_;
        }
        step_(u, dt);

        if (cfg_.adaptive) {
            const bool pushed = spectrum_.add(a_up, dt);
            since_ident_ += dt;
            if (pushed && spectrum_.full() && since_ident_ >= cfg_.ident_period_s) {
                identify_();
                since_ident_ = T(0);
            }
            const T a = std::min(T(1), dt / cfg_.param_tau_s);
            wc_ = std::exp(std::log(wc_) + a * (std::log(wc_target_) - std::log(wc_)));
        }
        return heave_;
    }

    T heave() const { return heave_; }
    T heaveRate() const { return vel_; }
    T cutoff() const { return wc_; }
    T cutoffTarget() const { return wc_target_; }
    T dominantFrequency() const { return wp_; }
    T meanWaveHeight() const { return Ap_; }

private:
    void identify_() {
        const auto pk = spectrum_.peaks(cfg_.f_min_hz, cfg_.f_max_hz, T(0), 1);
        if (pk.empty()) return;
        const T wp = pk.front().w;

        // Eq. 14 on the buffer, mean removed.
        const auto buf = spectrum_.samples();
        T mean = T(0);
        for (T v : buf) mean += v;
        mean /= T(buf.size());
        T s2 = T(0);
        for (T v : buf) s2 += (v - mean) * (v - mean);
        const T Ap = std::sqrt(T(2) * s2 / T(buf.size())) / (wp * wp);
        if (!(Ap > T(0))) return;
        wp_ = wp;
        Ap_ = Ap;

        T target;
        if (!cfg_.bias_instability_term) {
            target = Law::optimalCutoff(cfg_.noise_density, wp, Ap);
        } else {
            constexpr int N = 160;
            const T lo = std::log(cfg_.wc_min), hi = std::log(cfg_.wc_max);
            T best = std::numeric_limits<T>::infinity();
            target = wc_target_;
            for (int i = 0; i < N; ++i) {
                const T w = std::exp(lo + (hi - lo) * T(i) / T(N - 1));
                const T j = Law::errorBound(w, wp, Ap, cfg_.noise_density, cfg_.bias_rw_density);
                if (j < best) { best = j; target = w; }
            }
        }
        wc_target_ = std::clamp(target, cfg_.wc_min, cfg_.wc_max);
    }

    // s_ = [x1, x1', x2, x2'].
    void deriv_(const T* s, T u, T* ds) const {
        const T a = std::sqrt(T(2)) * wc_;
        const T b = wc_ * wc_;
        ds[0] = s[1];
        ds[1] = u - a * s[1] - b * s[0];
        ds[2] = s[3];
        ds[3] = s[1] - a * s[3] - b * s[2];
    }

    void step_(T u, T dt) {
        T k1[4], k2[4], k3[4], k4[4], tmp[4];
        deriv_(s_, u, k1);
        for (int i = 0; i < 4; ++i) tmp[i] = s_[i] + T(0.5) * dt * k1[i];
        deriv_(tmp, u, k2);
        for (int i = 0; i < 4; ++i) tmp[i] = s_[i] + T(0.5) * dt * k2[i];
        deriv_(tmp, u, k3);
        for (int i = 0; i < 4; ++i) tmp[i] = s_[i] + dt * k3[i];
        deriv_(tmp, u, k4);
        for (int i = 0; i < 4; ++i)
            s_[i] += dt / T(6) * (k1[i] + T(2) * k2[i] + T(2) * k3[i] + k4[i]);
        Law::output(wc_, wp_, s_, heave_, vel_);
    }

    Config cfg_;
    AccelSpectrum<T> spectrum_;
    T s_[4] = {};
    T wc_ = T(0.1), wc_target_ = T(0.1);
    T wp_ = T(NAN), Ap_ = T(NAN);
    T since_ident_ = T(0);
    T mean_ = T(0);
    long mean_n_ = 0;
    T heave_ = T(0), vel_ = T(0);
};

} // namespace heave_only
