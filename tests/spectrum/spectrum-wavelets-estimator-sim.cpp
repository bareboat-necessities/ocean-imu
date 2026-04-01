#include <iostream>

#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>

#include "spectrum/WaveSpectrumEstimator_Wavelets.h"

#include "SpectrumEstimatorSimShared.h"

const float g_std = 9.80665f;

int main() {
    constexpr float dt = 1.0f / 200.0f;
    constexpr int Nfreq = 32;
    constexpr int Nblock = 384;
    using Estimator = WaveSpectrumEstimator<Nfreq, Nblock>;

    std::cout << "Wavelet spectrum estimator simulation (IMU -> estimator vs ref spectrum files)\n";

    bool all_quality_ok = true;
    for (const auto& fname : SpectrumEstimatorSimShared::discover_wave_data_files()) {
        all_quality_ok = SpectrumEstimatorSimShared::process_wave_file<Estimator, Nfreq>(
                             fname, dt, 4321u, "spectrum_wavelets_")
                         && all_quality_ok;
    }

    return all_quality_ok ? 0 : 1;
}
