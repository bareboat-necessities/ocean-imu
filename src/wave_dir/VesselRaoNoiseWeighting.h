#pragma once

#include <algorithm>
#include <cmath>
#include "VesselRaoEqualizer.h"

namespace wave_direction {

// Optional measurement weighting for a known vessel response. The wave-band
// scale and frequency must come from the previous measurement-only schedule.
// This is a regularization weight, not a new estimate of sensor noise.
struct VesselRaoNoiseWeighting {
    VesselRaoEqualizer::Config response{};
    float max_std_scale = 1.0f; // disabled; a deployment must opt in
    float transition_snr = 2.0f;

    Eigen::Vector3f scales(float vertical_std, float frequency_hz,
                           const Eigen::Vector3f& nominal_noise_std) const {
        Eigen::Vector3f result = Eigen::Vector3f::Ones();
        if (!(std::isfinite(max_std_scale) && max_std_scale > 1.0f &&
              std::isfinite(transition_snr) && transition_snr > 0.0f &&
              std::isfinite(vertical_std) && vertical_std >= 0.0f &&
              std::isfinite(frequency_hz) && frequency_hz > 0.0f)) return result;
        const double fT = frequency_hz * response.heave_period_s;
        const double inverse_heave = std::hypot(1.0 - fT * fT,
            2.0 * response.heave_damping * fT);
        const double omega = 6.283185307179586 * frequency_hz;
        const double taus[2] = {response.horizontal_x_tau_s, response.horizontal_y_tau_s};
        const float maximum = std::min(max_std_scale, 4.0f);
        for (int i = 0; i < 2; ++i) {
            if (!(std::isfinite(taus[i]) && taus[i] >= 0.0 &&
                  std::isfinite(inverse_heave) && nominal_noise_std[i] > 0.0f)) continue;
            // Unit horizontal projection is conservative when direction is
            // unobservable. Do not use simulated angle or elevation truth.
            const double wt = omega * taus[i];
            const double snr = vertical_std * inverse_heave /
                ((1.0 + wt * wt) * nominal_noise_std[i]);
            const double u = std::clamp(snr / transition_snr - 1.0, 0.0, 1.0);
            const double weight = 1.0 - u * u * (3.0 - 2.0 * u);
            result[i] = static_cast<float>(std::sqrt(1.0 +
                (maximum * maximum - 1.0) * weight));
        }
        return result;
    }
};

} // namespace wave_direction
