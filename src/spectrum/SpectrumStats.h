#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace SpectrumStats {

// Keep this tiny helper only if other code still wants it.
// Safe to leave here; it is generic and harmless.
template <typename T>
inline double safe_log(T v) {
    return std::log(std::max(1e-18, static_cast<double>(v)));
}

// Zeroth spectral moment:
//   m0 = integral S(f) df
// with per-bin widths df[i].
template <int Nfreq, typename SpectrumLike>
inline double compute_m0(const SpectrumLike& spectrum,
                         const std::array<double, Nfreq>& df) {
    double m0 = 0.0;
    for (int i = 0; i < Nfreq; ++i) {
        const double s = std::max(0.0, static_cast<double>(spectrum[i]));
        m0 += s * df[i];
    }
    return m0;
}

// Significant wave height:
//   Hs = 4 * sqrt(m0)
template <int Nfreq, typename SpectrumLike>
inline double compute_hs(const SpectrumLike& spectrum,
                         const std::array<double, Nfreq>& df) {
    const double m0 = compute_m0<Nfreq>(spectrum, df);
    return 4.0 * std::sqrt(std::max(m0, 0.0));
}

} // namespace SpectrumStats
