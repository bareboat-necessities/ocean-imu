#pragma once
// Replay truth only: the declared track of oceanography-waves-lib v1.3.0.
// Never supply this reference to an estimator or remove its acceleration
// from the sensor input. World axes are north, west, up.
#include <array>
#include <cmath>

struct UnderwayReference {
    std::array<double, 3> displacement{}, velocity{}, acceleration{};

    static UnderwayReference at(double time) {
        constexpr double pi = 3.14159265358979323846;
        constexpr double speed = 5.0 * 1852.0 / 3600.0;
        constexpr double amplitude = 3.0 * pi / 180.0;
        constexpr double frequency = 2.0 * pi / 120.0;
        // Fourier-Bessel primitive of (cos(A sin(wt)), sin(A sin(wt))).
        // Same analytic track as the upstream generator, without quadrature
        // drift or subtracting a fitted trend from a scored filter output.
        static const auto bessel = [] {
            std::array<double, 13> values{};
            for (int n = 0; n <= 12; ++n) {
                double term = 1.0;
                for (int j = 1; j <= n; ++j) term *= amplitude / (2 * j);
                values[static_cast<std::size_t>(n)] = term;
                for (int j = 1; j <= 16; ++j) {
                    term *= -amplitude * amplitude / (4 * j * (j + n));
                    values[static_cast<std::size_t>(n)] += term;
                }
            }
            return values;
        }();
        const double phase = frequency * time;
        double x = bessel[0] * time, y = 0.0;
        for (int n = 1; n <= 12; ++n) {
            const double c = 2 * bessel[static_cast<std::size_t>(n)] / (n * frequency);
            if (n % 2 == 0) x += c * std::sin(n * phase);
            else y += c * (1 - std::cos(n * phase));
        }
        const double heading = amplitude * std::sin(phase);
        const double rate = amplitude * frequency * std::cos(phase);
        return {{speed * x, speed * y, 0.0},
                {speed * std::cos(heading), speed * std::sin(heading), 0.0},
                {-speed * rate * std::sin(heading), speed * rate * std::cos(heading), 0.0}};
    }
};
