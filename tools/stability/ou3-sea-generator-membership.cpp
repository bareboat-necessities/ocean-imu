// Execute the pinned physical generator through its public API. This observes
// its spectrum; it is neither an estimator nor a replacement BRMM provider.
// This audit utility requires the external pinned generator headers. CI builds
// it explicitly with warnings as errors, separately from the filter host TUs.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include "PiersonMoskowitzStokes3D_Waves.h"

int main()
{
    // Exact argument types/values of data-sim/waves_sim.cpp's medium scenario.
    constexpr float hs = 1.5f;
    constexpr float tp = 5.7f;
    constexpr float direction = -30.0f;
    constexpr float gravity = 9.80665f;
    auto spreading = std::make_shared<Cosine2sRandomizedDistribution>(
        direction * std::numbers::pi / 180.0, 10.0, 42u);
    PMStokesN3dWaves<128, 3> model(hs, tp, spreading, 0.02, 0.8, gravity, 42u);
    const double f = model.frequencies()(127);
    const double a = model.amplitudes()(127); // after Stokes renormalization
    const double omega = 2.0 * std::numbers::pi * f;
    const double k = omega * omega / gravity;
    const double c3 = (3.0 / 8.0) * a * std::pow(k * a, 2);
    if (!(a > 0 && c3 > 0 && std::isfinite(c3))) return 1;
    std::cout << std::setprecision(17)
              << "{\"frequency_count\":128,\"order\":3,"
              << "\"highest_fundamental_hz\":" << f << ','
              << "\"highest_third_harmonic_hz\":" << 3.0 * f << ','
              << "\"highest_fundamental_amplitude_m\":" << a << ','
              << "\"highest_third_harmonic_amplitude_m\":" << c3 << ','
              << "\"highest_third_harmonic_acceleration_amplitude_mps2\":"
              << 9.0 * omega * omega * c3 << ','
              << "\"gravity_mps2\":" << gravity << ','
              << "\"atoms\":[";
    for (int i = 0; i < 128; ++i) {
        const double fi = model.frequencies()(i);
        const double ai = model.amplitudes()(i);
        const double wi = 2.0 * std::numbers::pi * fi;
        const double ki = wi * wi / gravity;
        if (i != 0) std::cout << ',';
        std::cout << "{\"frequency_hz\":" << fi << ",\"amplitude_m\":" << ai
                  << ",\"second_coefficient_m\":" << 0.5 * ki * ai * ai
                  << ",\"third_coefficient_m\":" << (3.0 / 8.0) * ki * ki * ai * ai * ai << '}';
    }
    std::cout << "]}\n";
}
