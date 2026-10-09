#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Heave-only filter: the zero-displacement heave filter of

    M. Richter, K. Schneider, D. Walser, O. Sawodny, "Real-time heave motion
    estimation using adaptive filtering techniques", 19th IFAC World Congress,
    pp. 10119-10125, 2014, Section 3.2 (equation numbers are the paper's).

  The standard filter's main real-time error is its phase lead (Eq. 7).
  Moving one of its zeros off the origin (Eq. 23)

      Hzd(s) = s (s + a) / (s^2 + 2 zeta w_c s + w_c^2)^2,   zeta = 1/sqrt(2),

  with (Eq. 25)

      a = 2 sqrt2 w_c (1 - w_c^2 / w_p^2)  ->  2 sqrt2 w_c   for w_p >> w_c,

  cuts the error on a tone to |e| -> 4 (w_c/w_p)^2 (Eq. 26) at the cost of
  about nine times the white-noise variance, (9 sqrt2/16) sigma_n^2 / w_c^3
  (Eq. 27).  The cutoff minimises (Eq. 28)

      Jzd = 16 A_p^2 (w_c/w_p)^4 + (9 sqrt2/16) sigma_n^2 / w_c^3,

  whose minimiser is (Eq. 29)

      w_c^7 = (27 sqrt2 / 1024) sigma_n^2 w_p^4 / A_p^2.

  With a single zero at the origin Hzd passes far more bias drift than H.
  The bias-instability extension adds its exact variance for a bias random
  walk of intensity q_b,

      q_b sqrt2/16 (w_c^-5 + 3 a^2 w_c^-7).

  s/D^2 is the second section's position state of the shared engine, so the
  output is the standard heave plus a times that state.  Identification,
  retuning and the extensions are in AdaptiveHeaveFilter.h.
*/

#include <algorithm>
#include <cmath>

#include "AdaptiveHeaveFilter.h"

namespace heave_only {

template <typename T>
struct RichterZeroDisplacementLaw {
    // Eq. 25; before w_p is identified the asymptotic value is used.
    static T zeroOffset(T wc, T wp) {
        const T r = std::isfinite(wp) && wp > T(0) ? wc / wp : T(0);
        return T(2) * std::sqrt(T(2)) * wc * std::max(T(0), T(1) - r * r);
    }

    // Eq. 29.
    static T optimalCutoff(T noise_density, T wp, T Ap) {
        return std::pow(T(27) * std::sqrt(T(2)) / T(1024) * noise_density *
                        wp * wp * wp * wp / (Ap * Ap), T(1) / T(7));
    }

    // Eq. 28, plus the bias-instability variance of Hzd.
    static T errorBound(T wc, T wp, T Ap, T noise_density, T bias_rw_density) {
        const T r = wc / wp;
        const T a = zeroOffset(wc, wp);
        const T wc3 = wc * wc * wc, wc5 = wc3 * wc * wc, wc7 = wc5 * wc * wc;
        return T(16) * Ap * Ap * r * r * r * r +
               T(9) * std::sqrt(T(2)) / T(16) * noise_density / wc3 +
               bias_rw_density * std::sqrt(T(2)) / T(16) * (T(1) / wc5 + T(3) * a * a / wc7);
    }

    // Heave = (s^2 + a s)/D^2 u = x2' + a x2, heave rate = x2'' + a x2'.
    static void output(T wc, T wp, const T s[4], T& heave, T& rate) {
        const T a = zeroOffset(wc, wp);
        heave = s[3] + a * s[2];
        rate = s[1] - std::sqrt(T(2)) * wc * s[3] - wc * wc * s[2] + a * s[3];
    }
};

template <typename T = double>
using RichterHeaveFilter = AdaptiveHeaveFilter<T, RichterZeroDisplacementLaw<T>>;

} // namespace heave_only
