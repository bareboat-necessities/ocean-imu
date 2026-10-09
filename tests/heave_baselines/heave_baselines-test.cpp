#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <random>
#include <stdexcept>

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

constexpr double kDt = 0.005;
constexpr double kW1 = 2.0 * std::numbers::pi / 8.0;
constexpr double kW2 = 2.0 * std::numbers::pi / 5.0;

// Two-tone heave (1 m at 8 s, 0.3 m at 5 s) with a constant accelerometer
// bias and white noise at the simulator's 200 Hz sensor level.  Returns the
// heave RMS error over the last half of a 20-minute record, relative to the
// heave RMS.
template <typename Filter>
double relative_rms(Filter& f, double bias) {
    std::mt19937 rng(7);
    std::normal_distribution<double> n(0.0, 0.0148);
    const int steps = static_cast<int>(1200.0 / kDt);
    double se = 0.0, sr = 0.0;
    for (int k = 0; k < steps; ++k) {
        const double t = k * kDt;
        const double z = std::cos(kW1 * t + 0.3) + 0.3 * std::cos(kW2 * t + 1.1);
        const double a = -kW1 * kW1 * std::cos(kW1 * t + 0.3)
                         - 0.3 * kW2 * kW2 * std::cos(kW2 * t + 1.1);
        const double zh = f.update(a + bias + n(rng), kDt);
        if (t > 600.0) { se += (zh - z) * (zh - z); sr += z * z; }
    }
    return std::sqrt(se / sr);
}

} // namespace

int main() {
    using G = GodhavnHeaveFilter<double>;
    {
        // Richter et al. 2014, Eq. 7: e -> 2 sqrt2 w_c / w_p for w_p >> w_c.
        const double r = G::relativeError(0.01, 1.0) / (2.0 * std::sqrt(2.0) * 0.01);
        check(std::abs(r - 1.0) < 0.01, "godhavn Eq. 7 asymptote ratio", r, 0.01);
    }
    {
        // Eq. 12 is the minimiser of Eq. 11.
        const double q = 1.1e-6, wp = 0.8, Ap = 1.2;
        const double wc = G::optimalCutoff(q, wp, Ap);
        const double J0 = G::errorBound(wc, wp, Ap, q, 0.0);
        const bool ok = J0 < G::errorBound(wc * 1.02, wp, Ap, q, 0.0) &&
                        J0 < G::errorBound(wc / 1.02, wp, Ap, q, 0.0);
        check(ok, "godhavn Eq. 12 minimises Eq. 11, w_c", wc, 0.0);
    }
    {
        // Richter et al. 2014, Eq. 29 minimises Eq. 28.
        const double q = 1.1e-6, wp = 0.8, Ap = 1.2;
        const double wc = G::optimalCutoffZd(q, wp, Ap);
        const double J0 = G::errorBoundZd(wc, wp, Ap, q, 0.0);
        const bool ok = J0 < G::errorBoundZd(wc * 1.02, wp, Ap, q, 0.0) &&
                        J0 < G::errorBoundZd(wc / 1.02, wp, Ap, q, 0.0);
        check(ok, "zero-displacement Eq. 29 minimises Eq. 28, w_c", wc, 0.0);
    }
    {
        // Eq. 26: with a from Eq. 25 the zero-displacement filter's error on a
        // tone falls to |e| -> 4 (w_c/w_p)^2, far below the standard filter.
        // Noise-free: Hzd amplifies white noise ~9x more (Eq. 27), which would
        // swamp the structural error this checks.
        G::Config cfg;
        cfg.adaptive = false;
        cfg.wc_init = 0.07;
        cfg.zero_displacement = true;
        G g(cfg);
        double se = 0.0, sr = 0.0;
        for (int k = 0; k < static_cast<int>(1200.0 / kDt); ++k) {
            const double t = k * kDt;
            const double z = std::cos(kW1 * t);
            const double zh = g.update(-kW1 * kW1 * z + 0.05, kDt);
            if (t > 600.0) { se += (zh - z) * (zh - z); sr += z * z; }
        }
        // a uses the unidentified w_p (asymptotic Eq. 25), so allow 2x of Eq. 26.
        const double r = std::sqrt(se / sr);
        const double theory = 4.0 * (0.07 / kW1) * (0.07 / kW1);
        check(r < 2.0 * theory && r < 0.5 * G::relativeError(0.07, kW1),
              "zero-displacement tone error vs Eq. 26", r, 2.0 * theory);
    }
    {
        // Fixed cutoff: the bias is rejected and what remains is the phase
        // lead of Eq. 7 on each tone.
        G::Config cfg;
        cfg.adaptive = false;
        cfg.wc_init = 0.07;
        G g(cfg);
        const double r = relative_rms(g, 0.05);
        const double e1 = G::relativeError(0.07, kW1);
        const double e2 = 0.3 * G::relativeError(0.07, kW2);
        const double theory = std::sqrt((e1 * e1 + e2 * e2) / 1.09);
        check(std::abs(r / theory - 1.0) < 0.05, "godhavn heave rms / Eq. 7 theory",
              r / theory, 1.05);
    }
    {
        // Adaptation identifies the dominant heave frequency and the Eq. 14
        // wave height of the 8 s tone.
        G g;
        relative_rms(g, 0.05);
        check(std::abs(g.dominantFrequency() / kW1 - 1.0) < 0.03,
              "godhavn identified w_p / true", g.dominantFrequency() / kW1, 1.03);
        // Eq. 14 attributes all acceleration energy to w_p.
        const double a_rms2 = 0.5 * (std::pow(kW1, 4) + std::pow(0.3 * kW2 * kW2, 2));
        const double Ap = std::sqrt(2.0 * a_rms2) / (kW1 * kW1);
        check(std::abs(g.meanWaveHeight() / Ap - 1.0) < 0.05,
              "godhavn Eq. 14 A_p / expected", g.meanWaveHeight() / Ap, 1.05);
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
    {
        // The radix-2 FFT needs a power-of-two length.
        bool threw = false;
        try {
            KuchlerHeaveEKF<double>::Config cfg;
            cfg.fft_size = 500;
            KuchlerHeaveEKF<double> k(cfg);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "non-power-of-two fft_size rejected", threw, 1);
    }
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
