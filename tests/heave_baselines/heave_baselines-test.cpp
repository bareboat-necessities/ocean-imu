#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <random>

/*
  Copyright 2026, Mikhail Grushinskiy
*/

#define EIGEN_NON_ARDUINO

#include "heave_baselines/GodhavnHeaveFilter.h"
#include "heave_baselines/KuchlerHeaveEKF.h"

using namespace heave_baselines;

namespace {

int failures = 0;

void check(bool ok, const char* what, double value, double limit) {
    std::cout << (ok ? "PASS " : "FAIL ") << what << " = " << value
              << " (limit " << limit << ")\n";
    if (!ok) ++failures;
}

// Two-tone heave with a constant accelerometer bias and white noise at the
// simulator's 200 Hz sensor level.  Returns the steady-state heave RMS error
// over the last half of a 20-minute record, relative to the heave RMS.
template <typename Filter>
double relative_rms(Filter& f, double bias) {
    const double dt = 0.005;
    const double A1 = 1.0, w1 = 2.0 * std::numbers::pi / 8.0;
    const double A2 = 0.3, w2 = 2.0 * std::numbers::pi / 5.0;
    std::mt19937 rng(7);
    std::normal_distribution<double> n(0.0, 0.0148);
    const int steps = static_cast<int>(1200.0 / dt);
    double se = 0.0, sr = 0.0;
    for (int k = 0; k < steps; ++k) {
        const double t = k * dt;
        const double z = A1 * std::cos(w1 * t + 0.3) + A2 * std::cos(w2 * t + 1.1);
        const double a = -A1 * w1 * w1 * std::cos(w1 * t + 0.3)
                         - A2 * w2 * w2 * std::cos(w2 * t + 1.1);
        const double zh = f.update(a + bias + n(rng), dt);
        if (t > 600.0) { se += (zh - z) * (zh - z); sr += z * z; }
    }
    return std::sqrt(se / sr);
}

} // namespace

int main() {
    {
        // The error budget's sinusoid term tends to (2 sqrt2 w_c / w_p)^2 A^2/2
        // for w_c << w_p (Richter et al. 2014, Eq. 9).
        const double wc = 0.01, wp = 1.0, A = 1.0;
        const double J = GodhavnHeaveFilter<double>::errorBudget(wc, A, wp, 0.0, 0.0);
        const double ref = 0.5 * A * A * 8.0 * (wc / wp) * (wc / wp);
        check(std::abs(J / ref - 1.0) < 0.02, "godhavn sinusoid error term ratio", J / ref, 0.02);
    }
    {
        // The bias is rejected and what remains is the high-pass phase lead,
        // |1 - s^2 H| ~ 2 sqrt2 w_c / w per tone, at the adapted cutoff.
        GodhavnHeaveFilter<double> g;
        const double r = relative_rms(g, 0.05);
        const double wc = g.cutoff();
        check(wc > 0.04 && wc < 0.12, "godhavn adapted cutoff rad/s", wc, 0.12);
        const double e1 = 1.0 * 2.0 * std::sqrt(2.0) * wc / (2.0 * std::numbers::pi / 8.0);
        const double e2 = 0.3 * 2.0 * std::sqrt(2.0) * wc / (2.0 * std::numbers::pi / 5.0);
        const double theory = std::sqrt((e1 * e1 + e2 * e2) / (1.0 + 0.09));
        check(std::abs(r / theory - 1.0) < 0.10, "godhavn heave rms / phase-lead theory",
              r / theory, 1.10);
    }
    {
        // Two undamped tones are exactly the EKF's model: with small process
        // noise it must lock both frequencies and the offset.
        KuchlerHeaveEKF<double>::Config cfg;
        cfg.zeta_q = 0.0005;
        cfg.omega_rw = 0.0005;
        KuchlerHeaveEKF<double> k(cfg);
        const double r = relative_rms(k, 0.05);
        check(r < 0.03, "kuchler relative heave rms", r, 0.03);
        check(std::abs(k.offset() - 0.05) < 0.005, "kuchler offset error m/s^2",
              std::abs(k.offset() - 0.05), 0.005);
        check(k.modeCount() >= 1, "kuchler active modes", k.modeCount(), 1);
    }
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
