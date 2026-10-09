#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Literature baseline: the standard adaptive heave filter of

    J.-M. Godhavn, "Adaptive tuning of heave filter in motion sensor",
    OCEANS'98, vol. 1, pp. 174-178, 1998,

  implemented as specified by Richter, Schneider, Walser and Sawodny,
  "Real-time heave motion estimation using adaptive filtering techniques",
  19th IFAC World Congress, pp. 10119-10125, 2014, Sections 2.2-2.4 (the
  original paper is not openly available; equation numbers below are
  Richter's).

  Heave is the levelled vertical acceleration through (Eq. 3)

      H(s) = s^2 / (s^2 + 2 zeta w_c s + w_c^2)^2,   zeta = 1/sqrt(2),

  a double integrator in series with two second-order Butterworth high
  passes.  The double zero at s = 0 rejects a constant bias; the price is a
  phase lead in the wave band, |1 - s^2 H(j w_p)| -> 2 sqrt2 w_c / w_p
  (Eq. 7).  H is realised as two identical sections
  x'' + 2 zeta w_c x' + w_c^2 x = u, y = x', integrated with RK4; the states
  carry across retuning.

  Adaptation (Section 2.4, Fig. 4).  Every ident_period_s the FFT of the
  buffered acceleration gives the dominant heave frequency w_p (peak of the
  heave spectrum) and the buffer gives the mean wave height (Eq. 14)

      A_p = sqrt( 2/N sum_j a_j^2 / w_p^4 ).

  The cutoff minimises the total-error bound (Eq. 11)

      J = 8 A_p^2 (w_c / w_p)^2 + sigma_n^2 / (2^(7/2) w_c^3),

  whose minimiser is (Eq. 12)

      w_c,opt = 2^(-3/2) (3 sigma_n^2 w_p^2 / A_p^2)^(1/5),

  with sigma_n^2 the accelerometer white-noise spectral density from the
  datasheet (Eq. 8).  The determined parameters are strongly low-pass
  filtered (here: relaxation of log w_c with param_tau_s).

  Zero-displacement variant (zero_displacement, Richter et al. Section 3.2):
  one zero of H moves off the origin, Hzd(s) = s (s + a) / D(s)^2 (Eq. 23),
  which removes most of the phase lead; a follows Eq. 25,
  a = 2 sqrt2 w_c (1 - w_c^2 / w_p^2), and the cutoff minimises (Eq. 28)
  Jzd = 16 A_p^2 (w_c/w_p)^4 + (9 sqrt2/16) sigma_n^2 / w_c^3, i.e. (Eq. 29)
  w_c^7 = (27 sqrt2 / 1024) sigma_n^2 w_p^4 / A_p^2.  Since s/D^2 is already
  the second section's position state, p_zd = p + a x2.

  Two extensions, both off by default and not in the papers:
    bias_instability_term  adds the variance q_b / (2^(7/2) w_c^5) that an
                           accelerometer bias random walk of intensity q_b
                           leaves after H (exact for this H), and minimises
                           the extended J numerically (for Hzd the term is
                           q_b sqrt2/16 (w_c^-5 + 3 a^2 w_c^-7), also exact);
    subtract_input_mean    removes a slow running mean of the input ahead of
                           H so a retune does not move the DC state.
*/

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

#include "heave_baselines/AccelSpectrum.h"

namespace heave_baselines {

template <typename T = double>
class GodhavnHeaveFilter {
    static_assert(std::is_floating_point<T>::value,
                  "GodhavnHeaveFilter<T>: T must be a floating-point type.");

public:
    struct Config {
        bool adaptive = true;
        T wc_init = T(0.10);          // rad/s; the fixed cutoff when !adaptive
        T wc_min = T(0.005);          // rad/s
        T wc_max = T(1.0);            // rad/s
        // Accelerometer white-noise spectral density sigma_n^2, (m/s^2)^2/Hz
        // (= sample variance * sample period).
        T noise_density = T(1.1e-6);

        // Identification (Fig. 4).
        T fs_ident_hz = T(4);
        int fft_size = 512;           // 128 s at 4 Hz
        T ident_period_s = T(30);
        T f_min_hz = T(0.04);
        T f_max_hz = T(1.0);
        T param_tau_s = T(60);

        // Richter et al. 2014 zero-displacement filter (Section 3.2).
        bool zero_displacement = false;

        // Extensions (see the header comment).
        bool bias_instability_term = false;
        T bias_rw_density = T(2.5e-7);  // q_b, (m/s^2)^2/s
        bool subtract_input_mean = false;
        T input_mean_tau_s = T(300);
    };

    GodhavnHeaveFilter() : GodhavnHeaveFilter(Config{}) {}
    explicit GodhavnHeaveFilter(const Config& cfg)
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

    // Exact relative error |1 - s^2 H(s)| at s = j w (Eq. 7).
    static T relativeError(T wc, T w) {
        const T wc2 = wc * wc, w2 = w * w;
        return std::sqrt(wc2 * wc2 * wc2 * wc2 + T(8) * w2 * w2 * w2 * wc2) / (w2 * w2 + wc2 * wc2);
    }

    // Eq. 12.
    static T optimalCutoff(T noise_density, T wp, T Ap) {
        return std::pow(T(2), T(-1.5)) *
               std::pow(T(3) * noise_density * wp * wp / (Ap * Ap), T(0.2));
    }

    // Eq. 25.
    static T zeroDisplacementA(T wc, T wp) {
        const T r = std::isfinite(wp) && wp > T(0) ? wc / wp : T(0);
        return T(2) * std::sqrt(T(2)) * wc * std::max(T(0), T(1) - r * r);
    }

    // Eq. 29.
    static T optimalCutoffZd(T noise_density, T wp, T Ap) {
        return std::pow(T(27) * std::sqrt(T(2)) / T(1024) * noise_density *
                        wp * wp * wp * wp / (Ap * Ap), T(1) / T(7));
    }

    // Eq. 28, plus the optional bias-instability variance of Hzd.
    static T errorBoundZd(T wc, T wp, T Ap, T noise_density, T bias_rw_density) {
        const T r = wc / wp;
        const T a = zeroDisplacementA(wc, wp);
        const T wc3 = wc * wc * wc, wc5 = wc3 * wc * wc, wc7 = wc5 * wc * wc;
        return T(16) * Ap * Ap * r * r * r * r +
               T(9) * std::sqrt(T(2)) / T(16) * noise_density / wc3 +
               bias_rw_density * std::sqrt(T(2)) / T(16) * (T(1) / wc5 + T(3) * a * a / wc7);
    }

    // Eq. 11, plus the optional bias-instability variance.
    static T errorBound(T wc, T wp, T Ap, T noise_density, T bias_rw_density) {
        const T k = std::pow(T(2), T(3.5));
        const T r = wc / wp;
        return T(8) * Ap * Ap * r * r + noise_density / (k * wc * wc * wc) +
               bias_rw_density / (k * wc * wc * wc * wc * wc);
    }

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
            target = cfg_.zero_displacement ? optimalCutoffZd(cfg_.noise_density, wp, Ap)
                                            : optimalCutoff(cfg_.noise_density, wp, Ap);
        } else {
            constexpr int N = 160;
            const T lo = std::log(cfg_.wc_min), hi = std::log(cfg_.wc_max);
            T best = std::numeric_limits<T>::infinity();
            target = wc_target_;
            for (int i = 0; i < N; ++i) {
                const T w = std::exp(lo + (hi - lo) * T(i) / T(N - 1));
                const T j = cfg_.zero_displacement
                    ? errorBoundZd(w, wp, Ap, cfg_.noise_density, cfg_.bias_rw_density)
                    : errorBound(w, wp, Ap, cfg_.noise_density, cfg_.bias_rw_density);
                if (j < best) { best = j; target = w; }
            }
        }
        wc_target_ = std::clamp(target, cfg_.wc_min, cfg_.wc_max);
    }

    // s_ = [x1, x1', x2, x2'].  Heave = x2', heave rate = x2''.
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
        heave_ = s_[3];
        vel_ = s_[1] - std::sqrt(T(2)) * wc_ * s_[3] - wc_ * wc_ * s_[2];
        if (cfg_.zero_displacement) {
            const T a = zeroDisplacementA(wc_, wp_);
            vel_ += a * s_[3];
            heave_ += a * s_[2];
        }
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

} // namespace heave_baselines
