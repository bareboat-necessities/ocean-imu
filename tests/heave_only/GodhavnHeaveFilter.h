#pragma once

/*
  Copyright 2026, Mikhail Grushinskiy

  Heave-only filter: the standard adaptive heave filter of

    J.-M. Godhavn, "Adaptive tuning of heave filter in motion sensor",
    OCEANS'98, vol. 1, pp. 174-178, 1998,

  implemented as specified by Richter, Schneider, Walser and Sawodny,
  "Real-time heave motion estimation using adaptive filtering techniques",
  19th IFAC World Congress, pp. 10119-10125, 2014, Sections 2.2-2.4 (the
  original paper is not openly available; equation numbers are Richter's).

  Heave is the levelled vertical acceleration through (Eq. 3)

      H(s) = s^2 / (s^2 + 2 zeta w_c s + w_c^2)^2,   zeta = 1/sqrt(2),

  a double integrator in series with two second-order Butterworth high
  passes.  The double zero at s = 0 rejects a constant bias; the price is a
  phase lead in the wave band, |1 - s^2 H(j w_p)| -> 2 sqrt2 w_c / w_p
  (Eq. 7).

  The cutoff minimises the total-error bound (Eq. 11)

      J = 8 A_p^2 (w_c / w_p)^2 + sigma_n^2 / (2^(7/2) w_c^3),

  whose minimiser is (Eq. 12)

      w_c,opt = 2^(-3/2) (3 sigma_n^2 w_p^2 / A_p^2)^(1/5).

  The bias-instability extension adds q_b / (2^(7/2) w_c^5), the variance an
  accelerometer bias random walk of intensity q_b leaves after H (exact).
  Identification, retuning and the extensions are in AdaptiveHeaveFilter.h.
*/

#include <cmath>

#include "AdaptiveHeaveFilter.h"

namespace heave_only {

template <typename T>
struct GodhavnLaw {
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

    // Eq. 11, plus the bias-instability variance.
    static T errorBound(T wc, T wp, T Ap, T noise_density, T bias_rw_density) {
        const T k = std::pow(T(2), T(3.5));
        const T r = wc / wp;
        return T(8) * Ap * Ap * r * r + noise_density / (k * wc * wc * wc) +
               bias_rw_density / (k * wc * wc * wc * wc * wc);
    }

    // Heave = s^2/D^2 u = x2', heave rate = x2''.
    static void output(T wc, T /*wp*/, const T s[4], T& heave, T& rate) {
        heave = s[3];
        rate = s[1] - std::sqrt(T(2)) * wc * s[3] - wc * wc * s[2];
    }
};

template <typename T = double>
using GodhavnHeaveFilter = AdaptiveHeaveFilter<T, GodhavnLaw<T>>;

} // namespace heave_only
