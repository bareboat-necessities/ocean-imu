// Observe the pinned release's incident PM components and exact vessel RAO.
// Point diagnostics do not certify continuum BRMM/Normal-Live membership.
#define EIGEN_NON_ARDUINO
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include "PiersonMoskowitzStokes3D_Waves.h"
#include "VesselRao.h"
int main()
{
    auto spread = std::make_shared<Cosine2sRandomizedDistribution>(
        -30.0f * std::numbers::pi / 180.0, 10.0, 42u);
    PMStokesN3dWaves<128, 3> incident(1.5f, 5.7f, spread, .02, .8, 9.80665f, 42u);
    const auto waves = incident.incidentHarmonics();
    const VesselRao vessel(waves);
    std::cout << std::setprecision(17)
              << "{\"response_model\":\"VESSEL_RAO_28FT\",\"frequency_count\":"
              << waves.size() << ",\"particle_harmonics_applied\":false,\"atoms\":[";
    for (std::size_t i=0; i<waves.size(); ++i) {
        if (i) std::cout << ',';
        const auto& w=waves[i];
        const auto h=vessel.transfer(w.omega,w.wavenumber,w.direction);
        std::cout << "{\"frequency_hz\":" << w.omega/(2*std::numbers::pi)
                  << ",\"amplitude_m\":" << w.amplitude
                  << ",\"translation_gain_norm\":" << std::sqrt(std::norm(h[0])+std::norm(h[1])+std::norm(h[2]))
                  << ",\"yaw_gain_norm\":" << std::abs(h[5]) << '}';
    }
    std::cout << "]}\n";
}
