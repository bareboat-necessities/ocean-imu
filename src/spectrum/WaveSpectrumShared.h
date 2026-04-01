#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace WaveSpectrumShared {

inline double safe_log(double v) {
    return std::log(std::max(v, 1e-18));
}

struct Biquad {
    double b0 = 0.0, b1 = 0.0, b2 = 0.0, a1 = 0.0, a2 = 0.0;
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

inline void designLowpassBiquad(Biquad& bq, double f_cut, double Fs, double Q = 0.707) {
    const double Fc = f_cut / Fs;
    const double K = std::tan(M_PI * Fc);
    const double norm = 1.0 / (1.0 + K / Q + K * K);

    bq.b0 = K * K * norm;
    bq.b1 = 2.0 * bq.b0;
    bq.b2 = bq.b0;
    bq.a1 = 2.0 * (K * K - 1.0) * norm;
    bq.a2 = (1.0 - K / Q + K * K) * norm;
}

inline void designHighpassBiquad(Biquad& bq, double f_cut, double Fs, double Q = 0.707) {
    if (f_cut <= 0.0) {
        bq.b0 = 1.0;
        bq.b1 = 0.0;
        bq.b2 = 0.0;
        bq.a1 = 0.0;
        bq.a2 = 0.0;
        bq.reset();
        return;
    }

    const double Fc = f_cut / Fs;
    const double K = std::tan(M_PI * Fc);
    const double norm = 1.0 / (1.0 + K / Q + K * K);

    bq.b0 = 1.0 * norm;
    bq.b1 = -2.0 * norm;
    bq.b2 = 1.0 * norm;
    bq.a1 = 2.0 * (K * K - 1.0) * norm;
    bq.a2 = (1.0 - K / Q + K * K) * norm;
    bq.reset();
}

inline double biquad_mag2_raw(const Biquad& bq, double Omega_raw) {
    const double c1 = std::cos(Omega_raw);
    const double s1 = std::sin(Omega_raw);
    const double c2 = std::cos(2.0 * Omega_raw);
    const double s2 = std::sin(2.0 * Omega_raw);

    const double num_re = bq.b0 + bq.b1 * c1 + bq.b2 * c2;
    const double num_im = -(bq.b1 * s1 + bq.b2 * s2);

    const double den_re = 1.0 + bq.a1 * c1 + bq.a2 * c2;
    const double den_im = -(bq.a1 * s1 + bq.a2 * s2);

    const double num2 = num_re * num_re + num_im * num_im;
    const double den2 = den_re * den_re + den_im * den_im;

    return num2 / std::max(den2, 1e-16);
}

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

inline double smoothstep01(double x) {
    x = std::clamp(x, 0.0, 1.0);
    return x * x * (3.0 - 2.0 * x);
}

inline double lowfreq_taper(double f, double f_cut) {
    if (!(f_cut > 0.0)) return 1.0;

    const double f0 = 0.72 * f_cut;
    const double f1 = 1.20 * f_cut;

    if (f <= f0) {
        const double r = std::clamp(f / std::max(f0, 1e-12), 0.0, 1.0);
        return std::pow(r, 6.0);
    }
    if (f >= f1) return 1.0;

    return smoothstep01((f - f0) / std::max(f1 - f0, 1e-12));
}

template<int Nfreq, typename SpectrumLike>
inline void smooth_logfreq_3tap_inplace(SpectrumLike& spectrum,
                                        const std::array<double, Nfreq>& freqs) {
    if (Nfreq < 3) return;

    std::array<double, Nfreq> S{};
    std::array<double, Nfreq> Sout{};

    for (int i = 0; i < Nfreq; ++i) {
        S[i] = static_cast<double>(spectrum[i]);
    }

    auto wpair = [&](int i) {
        const double eps = 1e-12;
        const double x_im1 = (i > 0)
            ? std::log(std::max(freqs[i - 1], eps))
            : std::log(std::max(freqs[i], eps));
        const double x_i = std::log(std::max(freqs[i], eps));
        const double x_ip1 = (i < Nfreq - 1)
            ? std::log(std::max(freqs[i + 1], eps))
            : std::log(std::max(freqs[i], eps));

        const double dL = std::max(0.0, x_i - x_im1);
        const double dR = std::max(0.0, x_ip1 - x_i);

        const double k = 0.35;
        const double minC = 0.40;

        const double wL = k * dL;
        const double wR = k * dR;
        const double wC = std::max(minC, 1.0 - (wL + wR));
        const double W = wL + wC + wR;

        return std::array<double, 3>{wL / W, wC / W, wR / W};
    };

    for (int i = 0; i < Nfreq; ++i) {
        const auto w = wpair(i);
        const double Sm1 = (i > 0) ? S[i - 1] : S[i];
        const double Sp1 = (i < Nfreq - 1) ? S[i + 1] : S[i];
        Sout[i] = w[0] * Sm1 + w[1] * S[i] + w[2] * Sp1;
    }

    for (int i = 0; i < Nfreq; ++i) {
        spectrum[i] = Sout[i];
    }
}

template<int Nfreq>
inline double estimate_lowfreq_cut_from_accel(const std::array<double, Nfreq>& S_aa_meas,
                                              const std::array<double, Nfreq>& freqs,
                                              double Tblk,
                                              double hp_f0_hz) {
    std::array<double, Nfreq> E{};
    std::array<double, Nfreq> Es{};

    for (int i = 0; i < Nfreq; ++i) {
        E[i] = std::max(0.0, S_aa_meas[i]) * std::max(freqs[i], 1e-12);
    }

    Es[0] = 0.75 * E[0] + 0.25 * E[1];
    for (int i = 1; i < Nfreq - 1; ++i) {
        Es[i] = 0.25 * E[i - 1] + 0.50 * E[i] + 0.25 * E[i + 1];
    }
    Es[Nfreq - 1] = 0.25 * E[Nfreq - 2] + 0.75 * E[Nfreq - 1];

    const double f_floor = std::max({
        1.35 * freqs[0],
        2.8 / std::max(Tblk, 1e-12),
        2.2 * hp_f0_hz
    });

    int i_floor = 0;
    while (i_floor + 1 < Nfreq && freqs[i_floor + 1] < f_floor) ++i_floor;

    int i_peak = std::max(1, i_floor);
    double e_peak = Es[i_peak];
    for (int i = i_peak + 1; i < Nfreq - 1; ++i) {
        if (Es[i] > e_peak) {
            e_peak = Es[i];
            i_peak = i;
        }
    }

    if (!(e_peak > 0.0) || i_peak <= i_floor) {
        return f_floor;
    }

    int i_valley = i_floor;
    double e_valley = Es[i_floor];
    for (int i = i_floor + 1; i < i_peak; ++i) {
        if (Es[i] < e_valley) {
            e_valley = Es[i];
            i_valley = i;
        }
    }

    const bool valley_is_good =
        (i_valley > i_floor) &&
        (e_valley < 0.82 * e_peak) &&
        (freqs[i_valley] < 0.88 * freqs[i_peak]);

    double f_cut = valley_is_good
        ? freqs[i_valley]
        : std::max(f_floor, 0.58 * freqs[i_peak]);

    f_cut = std::clamp(f_cut, f_floor, 0.90 * freqs[i_peak]);
    return f_cut;
}

template<int Nfreq, typename SpectrumLike>
inline void suppress_lowfreq_from_cut_inplace(SpectrumLike& spectrum,
                                              const std::array<double, Nfreq>& freqs,
                                              double f_cut_hz) {
    if (Nfreq < 3 || !(f_cut_hz > 0.0)) return;

    int i_cut = 0;
    while (i_cut + 1 < Nfreq && freqs[i_cut + 1] <= f_cut_hz) ++i_cut;

    std::array<double, Nfreq> E{};
    for (int i = 0; i < Nfreq; ++i) {
        E[i] = std::max(0.0, static_cast<double>(spectrum[i])) * std::max(freqs[i], 1e-12);
    }

    const int i_ref = std::min(i_cut + 1, Nfreq - 1);
    const double f_ref = std::max(freqs[i_ref], 1e-12);
    const double E_ref = std::max(E[i_ref], 1e-18);

    double prev = E_ref;
    constexpr double shape_pow = 4.8;

    for (int i = i_cut; i >= 0; --i) {
        const double r = std::clamp(freqs[i] / f_ref, 0.0, 1.0);
        const double cap = E_ref * std::pow(r, shape_pow);

        double v = std::max(0.0, E[i]);
        v = std::min(v, cap);
        v = std::min(v, prev);

        E[i] = v;
        prev = v;
    }

    for (int i = 0; i < Nfreq; ++i) {
        spectrum[i] = E[i] / std::max(freqs[i], 1e-12);
    }
}

template<int Nfreq, typename SpectrumLike>
inline double estimate_fp_with_guard(const SpectrumLike& spectrum,
                                     const std::array<double, Nfreq>& freqs,
                                     double f_guard_hz) {
    int k0 = 0;
    const double f_guard = std::max(f_guard_hz, 1.05 * freqs[0]);
    while (k0 + 1 < Nfreq && freqs[k0 + 1] <= f_guard) ++k0;

    int k = k0;
    double vmax = -1.0;
    for (int i = k0; i < Nfreq; ++i) {
        const double s = std::max(0.0, static_cast<double>(spectrum[i]));
        if (s > vmax) {
            vmax = s;
            k = i;
        }
    }

    if (k <= 0 || k >= Nfreq - 1) return freqs[k];

    const double f0 = freqs[k - 1];
    const double f1 = freqs[k];
    const double f2 = freqs[k + 1];
    if (f0 <= 0.0 || f1 <= 0.0 || f2 <= 0.0) return f1;

    const double x0 = std::log(f0);
    const double x1 = std::log(f1);
    const double x2 = std::log(f2);

    const double y0 = safe_log(std::max(0.0, static_cast<double>(spectrum[k - 1])));
    const double y1 = safe_log(std::max(0.0, static_cast<double>(spectrum[k])));
    const double y2 = safe_log(std::max(0.0, static_cast<double>(spectrum[k + 1])));

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
