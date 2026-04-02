#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace WaveSpectrumShared {

// ------------------------------------------------------------
// Common small biquad
// ------------------------------------------------------------
struct Biquad {
    double b0 = 0.0, b1 = 0.0, b2 = 0.0;
    double a1 = 0.0, a2 = 0.0;
    double z1 = 0.0, z2 = 0.0;

    inline double process(double x) {
        const double y = b0 * x + z1;
        z1 = b1 * x - a1 * y + z2;
        z2 = b2 * x - a2 * y;
        return y;
    }

    inline void reset() {
        z1 = 0.0;
        z2 = 0.0;
    }
};

static constexpr double BIQUAD_Q = 0.707;

// ------------------------------------------------------------
// Basic helpers
// ------------------------------------------------------------
inline double safe_pos(double x, double floor = 1e-12) {
    return std::max(x, floor);
}

inline double ema_alpha_for_f(double f,
                              double fmin,
                              double fmax,
                              double alpha_low,
                              double alpha_high) {
    double t = (f - fmin) / std::max(1e-12, (fmax - fmin));
    t = std::clamp(t, 0.0, 1.0);
    return alpha_low + (alpha_high - alpha_low) * t;
}

// ------------------------------------------------------------
// Frequency grid
// ------------------------------------------------------------
template<int Nfreq>
inline void build_log_frequency_grid(std::array<double, Nfreq>& freqs,
                                     std::array<double, Nfreq + 1>& f_edges,
                                     std::array<double, Nfreq>& df,
                                     double f_min = 0.04,
                                     double f_max = 1.2) {
    for (int i = 0; i <= Nfreq; ++i) {
        const double t = double(i) / double(Nfreq);
        f_edges[i] = f_min * std::pow(f_max / f_min, t);
    }

    for (int i = 0; i < Nfreq; ++i) {
        freqs[i] = std::sqrt(f_edges[i] * f_edges[i + 1]);
        df[i] = f_edges[i + 1] - f_edges[i];
    }
}

// ------------------------------------------------------------
// Biquad design / response
// ------------------------------------------------------------
inline void designLowpassBiquad(Biquad& bq, double f_cut, double Fs) {
    const double Fc = f_cut / Fs;
    const double K = std::tan(M_PI * Fc);
    const double norm = 1.0 / (1.0 + K / BIQUAD_Q + K * K);

    bq.b0 = K * K * norm;
    bq.b1 = 2.0 * bq.b0;
    bq.b2 = bq.b0;
    bq.a1 = 2.0 * (K * K - 1.0) * norm;
    bq.a2 = (1.0 - K / BIQUAD_Q + K * K) * norm;
    bq.reset();
}

inline void designHighpassBiquad(Biquad& bq, double f_cut, double Fs) {
    if (f_cut <= 0.0) {
        bq.b0 = 1.0; bq.b1 = 0.0; bq.b2 = 0.0;
        bq.a1 = 0.0; bq.a2 = 0.0;
        bq.reset();
        return;
    }

    const double Fc = f_cut / Fs;
    const double K = std::tan(M_PI * Fc);
    const double norm = 1.0 / (1.0 + K / BIQUAD_Q + K * K);

    bq.b0 = 1.0 * norm;
    bq.b1 = -2.0 * norm;
    bq.b2 = 1.0 * norm;
    bq.a1 = 2.0 * (K * K - 1.0) * norm;
    bq.a2 = (1.0 - K / BIQUAD_Q + K * K) * norm;
    bq.reset();
}

inline double biquad_mag2_raw(const Biquad& bq, double omega_raw) {
    const double c1 = std::cos(omega_raw);
    const double s1 = std::sin(omega_raw);
    const double c2 = std::cos(2.0 * omega_raw);
    const double s2 = std::sin(2.0 * omega_raw);

    const double num_re = bq.b0 + bq.b1 * c1 + bq.b2 * c2;
    const double num_im = -(bq.b1 * s1 + bq.b2 * s2);

    const double den_re = 1.0 + bq.a1 * c1 + bq.a2 * c2;
    const double den_im = -(bq.a1 * s1 + bq.a2 * s2);

    const double num2 = num_re * num_re + num_im * num_im;
    const double den2 = den_re * den_re + den_im * den_im;

    return num2 / std::max(den2, 1e-16);
}

// ------------------------------------------------------------
// Simple 3-tap smoothing helper for std::array
// ------------------------------------------------------------
template<int N>
inline void smooth_3tap_array(const std::array<double, N>& in,
                              std::array<double, N>& out) {
    if (N <= 0) return;
    if (N == 1) {
        out[0] = in[0];
        return;
    }

    out[0] = 0.75 * in[0] + 0.25 * in[1];
    for (int i = 1; i < N - 1; ++i) {
        out[i] = 0.25 * in[i - 1] + 0.50 * in[i] + 0.25 * in[i + 1];
    }
    out[N - 1] = 0.25 * in[N - 2] + 0.75 * in[N - 1];
}

// ------------------------------------------------------------
// Light log-frequency smoother with configurable peak keep
// ------------------------------------------------------------
template<int Nfreq, typename SpectrumLike>
inline void smooth_logfreq_3tap_custom_inplace(SpectrumLike& spectrum,
                                               const std::array<double, Nfreq>& freqs,
                                               double k,
                                               double minC,
                                               double peak_ratio = 1.10,
                                               double peak_keep = 0.90) {
    if (Nfreq < 3) return;

    std::array<double, Nfreq> in{};
    std::array<double, Nfreq> out{};

    for (int i = 0; i < Nfreq; ++i) {
        in[i] = std::max(0.0, double(spectrum[i]));
        out[i] = in[i];
    }

    for (int i = 0; i < Nfreq; ++i) {
        const double xim1 = std::log(std::max((i > 0) ? freqs[i - 1] : freqs[i], 1e-12));
        const double xi   = std::log(std::max(freqs[i], 1e-12));
        const double xip1 = std::log(std::max((i < Nfreq - 1) ? freqs[i + 1] : freqs[i], 1e-12));

        const double dL = std::max(0.0, xi - xim1);
        const double dR = std::max(0.0, xip1 - xi);

        double wL = k * dL;
        double wR = k * dR;
        double wC = std::max(minC, 1.0 - (wL + wR));
        const double W = wL + wC + wR;

        wL /= W;
        wC /= W;
        wR /= W;

        const double Sm1 = (i > 0) ? in[i - 1] : in[i];
        const double Si  = in[i];
        const double Sp1 = (i < Nfreq - 1) ? in[i + 1] : in[i];

        out[i] = wL * Sm1 + wC * Si + wR * Sp1;
    }

    for (int i = 1; i < Nfreq - 1; ++i) {
        const double li = in[i - 1];
        const double ci = in[i];
        const double ri = in[i + 1];

        const bool strong_local_peak =
            (ci > peak_ratio * std::max(li, 1e-18)) &&
            (ci > peak_ratio * std::max(ri, 1e-18));

        if (strong_local_peak) {
            out[i] = std::max(out[i], peak_keep * ci);
        }
    }

    for (int i = 0; i < Nfreq; ++i) {
        spectrum[i] = std::max(0.0, out[i]);
    }
}

// ------------------------------------------------------------
// Shared acceleration-energy peak finder
// This is only a low-level helper; policy stays estimator-specific.
// ------------------------------------------------------------
template<int Nfreq>
inline int find_accel_peak_index(const std::array<double, Nfreq>& S_aa_true_arr,
                                 const std::array<double, Nfreq>& freqs,
                                 double lowfreq_cut_hz,
                                 double floor_mult = 1.03) {
    if (Nfreq < 3) return 0;

    std::array<double, Nfreq> Eaa{};
    std::array<double, Nfreq> Eaa_s{};

    for (int i = 0; i < Nfreq; ++i) {
        Eaa[i] = std::max(0.0, S_aa_true_arr[i]) * std::max(freqs[i], 1e-12);
    }

    smooth_3tap_array<Nfreq>(Eaa, Eaa_s);

    int k0 = 0;
    const double f0 = std::max(lowfreq_cut_hz, floor_mult * freqs[0]);
    while (k0 + 1 < Nfreq && freqs[k0 + 1] <= f0) ++k0;

    int k_peak = k0;
    double v_peak = -1.0;
    for (int i = k0; i < Nfreq; ++i) {
        if (Eaa_s[i] > v_peak) {
            v_peak = Eaa_s[i];
            k_peak = i;
        }
    }
    return k_peak;
}

// ------------------------------------------------------------
// Guarded peak-frequency estimate
// ------------------------------------------------------------
template<int Nfreq, typename SpectrumLike>
inline double estimate_fp_with_guard(const SpectrumLike& spectrum,
                                     const std::array<double, Nfreq>& freqs,
                                     double lowfreq_guard_hz) {
    if (Nfreq == 0) return 0.0;

    int k0 = 0;
    const double f_guard = std::max(lowfreq_guard_hz, freqs[0]);
    while (k0 + 1 < Nfreq && freqs[k0 + 1] < f_guard) ++k0;

    int k = k0;
    double vmax = -1.0;
    for (int i = k0; i < Nfreq; ++i) {
        const double s = std::max(0.0, double(spectrum[i]));
        if (s > vmax) {
            vmax = s;
            k = i;
        }
    }

    if (!(vmax > 0.0)) return freqs[k0];
    if (k <= 0 || k >= Nfreq - 1) return freqs[k];

    const double f0 = freqs[k - 1];
    const double f1 = freqs[k];
    const double f2 = freqs[k + 1];
    if (!(f0 > 0.0 && f1 > 0.0 && f2 > 0.0)) return f1;

    const double x0 = std::log(f0);
    const double x1 = std::log(f1);
    const double x2 = std::log(f2);

    const double y0 = std::log(std::max(double(spectrum[k - 1]), 1e-18));
    const double y1 = std::log(std::max(double(spectrum[k]),     1e-18));
    const double y2 = std::log(std::max(double(spectrum[k + 1]), 1e-18));

    const double L0a = 1.0 / ((x0 - x1) * (x0 - x2));
    const double L1a = 1.0 / ((x1 - x0) * (x1 - x2));
    const double L2a = 1.0 / ((x2 - x0) * (x2 - x1));

    const double a = y0 * L0a + y1 * L1a + y2 * L2a;
    const double b =
        y0 * (-(x1 + x2) * L0a) +
        y1 * (-(x0 + x2) * L1a) +
        y2 * (-(x0 + x1) * L2a);

    if (a >= 0.0) return f1;

    const double x_peak = -b / (2.0 * a);
    const double f_peak = std::exp(x_peak);

    if (f_peak <= std::min(f0, f2) || f_peak >= std::max(f0, f2)) return f1;
    return f_peak;
}

} // namespace WaveSpectrumShared
