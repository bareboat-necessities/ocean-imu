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

template <int Nfreq, typename SpectrumLike>
inline double estimate_fp(const SpectrumLike& spectrum, const std::array<double, Nfreq>& freqs) {
    int k = 0;
    double vmax = 0.0;
    for (int i = 0; i < Nfreq; ++i) {
        const double s = static_cast<double>(spectrum[i]);
        if (s > vmax) {
            vmax = s;
            k = i;
        }
    }

    if (k == 0 || k == Nfreq - 1) return freqs[k];

    const double f0 = freqs[k - 1], f1 = freqs[k], f2 = freqs[k + 1];
    if (f0 <= 0.0 || f1 <= 0.0 || f2 <= 0.0) return f1;

    const double x0 = std::log(f0), x1 = std::log(f1), x2 = std::log(f2);
    const double y0 = safe_log(spectrum[k - 1]);
    const double y1 = safe_log(spectrum[k]);
    const double y2 = safe_log(spectrum[k + 1]);

    const double dx01 = x0 - x1, dx02 = x0 - x2, dx12 = x1 - x2;
    const double denom = (dx01 * dx02 * dx12);
    if (std::abs(denom) < 1e-18) return f1;

    const double L0a = 1.0 / (dx01 * dx02);
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
inline double estimate_tp(const SpectrumLike& spectrum, const std::array<double, Nfreq>& freqs) {
    const double fp = estimate_fp<Nfreq>(spectrum, freqs);
    return (fp > 0.0) ? (1.0 / fp) : 0.0;
}

} // namespace SpectrumStats
