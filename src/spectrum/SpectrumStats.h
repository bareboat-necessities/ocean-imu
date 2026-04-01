#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace SpectrumStats {

template <typename T>
inline double safe_log(T v) {
    return std::log(std::max(1e-18, static_cast<double>(v)));
}

template <int Nfreq, typename SpectrumLike>
inline double compute_hs(const SpectrumLike& spectrum, const std::array<double, Nfreq>& df) {
    double m0 = 0.0;
    for (int i = 0; i < Nfreq; ++i) {
        m0 += static_cast<double>(spectrum[i]) * df[i];
    }
    return 4.0 * std::sqrt(std::max(m0, 0.0));
}

// Find an adaptive low-frequency guard bin using E(f)=fS(f).
template <int Nfreq, typename SpectrumLike>
inline int lowfreq_guard_bin(const SpectrumLike& spectrum,
                             const std::array<double, Nfreq>& freqs) {
    if (Nfreq < 4) return 0;

    std::array<double, Nfreq> E{};
    std::array<double, Nfreq> Es{};

    for (int i = 0; i < Nfreq; ++i) {
        E[i] = std::max(0.0, static_cast<double>(spectrum[i])) * std::max(freqs[i], 1e-12);
    }

    Es[0] = 0.75 * E[0] + 0.25 * E[1];
    for (int i = 1; i < Nfreq - 1; ++i) {
        Es[i] = 0.25 * E[i - 1] + 0.50 * E[i] + 0.25 * E[i + 1];
    }
    Es[Nfreq - 1] = 0.25 * E[Nfreq - 2] + 0.75 * E[Nfreq - 1];

    int k_peak = 1;
    double v_peak = Es[1];
    for (int i = 2; i < Nfreq - 1; ++i) {
        if (Es[i] > v_peak) {
            v_peak = Es[i];
            k_peak = i;
        }
    }
    if (!(v_peak > 0.0) || k_peak <= 1) return 0;

    int k_valley = 0;
    double v_min = Es[0];
    for (int i = 1; i < k_peak; ++i) {
        if (Es[i] < v_min) {
            v_min = Es[i];
            k_valley = i;
        }
    }

    const bool left_edge_raised =
        (E[0] > 1.25 * std::max(v_min, 1e-18)) ||
        (Nfreq > 2 && E[1] > 1.10 * std::max(E[2], 1e-18));

    const bool valley_is_good =
        (k_valley > 0) &&
        (v_min < 0.72 * v_peak) &&
        (freqs[k_valley] < 0.85 * freqs[k_peak]);

    if (left_edge_raised && valley_is_good) {
        return k_valley;
    }

    // Fallback relative guard
    const double f_guard = 0.42 * freqs[k_peak];
    int k_guard = 0;
    for (int i = 0; i < k_peak; ++i) {
        if (freqs[i] <= f_guard) k_guard = i;
    }
    return k_guard;
}

template <int Nfreq, typename SpectrumLike>
inline double estimate_fp(const SpectrumLike& spectrum,
                          const std::array<double, Nfreq>& freqs) {
    int k0 = lowfreq_guard_bin<Nfreq>(spectrum, freqs);

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
    if (k <= k0) return freqs[k];

    const double f0 = freqs[k - 1], f1 = freqs[k], f2 = freqs[k + 1];
    if (f0 <= 0.0 || f1 <= 0.0 || f2 <= 0.0) return f1;

    const double x0 = std::log(f0), x1 = std::log(f1), x2 = std::log(f2);
    const double y0 = safe_log(spectrum[k - 1]);
    const double y1 = safe_log(spectrum[k]);
    const double y2 = safe_log(spectrum[k + 1]);

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

template <int Nfreq, typename SpectrumLike>
inline double estimate_tp(const SpectrumLike& spectrum,
                          const std::array<double, Nfreq>& freqs) {
    const double fp = estimate_fp<Nfreq>(spectrum, freqs);
    return (fp > 0.0) ? (1.0 / fp) : 0.0;
}

} // namespace SpectrumStats
