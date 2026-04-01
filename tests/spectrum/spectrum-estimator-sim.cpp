#include <iostream>

#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>

#include "spectrum/WaveSpectrumEstimator.h"

#include "SpectrumEstimatorSimShared.h"

const float g_std = 9.80665f;

int main() {
    constexpr float dt = 1.0f / 200.0f;
    constexpr int Nfreq = 32;
    constexpr int Nblock = 256;
    using Estimator = WaveSpectrumEstimator<Nfreq, Nblock>;

    std::cout << "Spectrum estimator simulation (IMU -> estimator vs ref spectrum files)\n";

    for (const auto& fname : SpectrumEstimatorSimShared::discover_wave_data_files()) {
        SpectrumEstimatorSimShared::process_wave_file<Estimator, Nfreq>(
            fname, dt, 1234u, "spectrum_estimator_");
    }

    return 0;
}
