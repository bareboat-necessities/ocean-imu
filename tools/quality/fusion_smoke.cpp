// Bounded host replay for the ASan/UBSan gate. Compile one production wrapper per job.
#define EIGEN_NON_ARDUINO
#include <Eigen/Dense>
#include <cmath>
#include <cstdio>

#if defined(FILTER_OU2)
#include "kalman_ou_ii/SeaStateFusionFilter_OU_II.h"
using Fusion = SeaStateFusion_OU_II<TrackerType::KALMANF>;
#elif defined(FILTER_OU3)
#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
#elif defined(FILTER_TFG)
#include "kalman_tfg/SeaStateFusionFilter_TFG.h"
using Fusion = ocean_imu::tfg::SeaStateFusionFilter_TFG<>;
#elif defined(FILTER_CARMA)
#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"
using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false>;
#elif defined(FILTER_TFG_CARMA)
#include "kalman_carma/SeaStateFusionFilter_TFG_CARMA3.h"
using Fusion = ocean_imu::carma::SeaStateFusion_TFG_CARMA3<TrackerType::KALMANF, 4>;
#elif defined(FILTER_SHARED)
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"
using Fusion = ocean_imu::carma_shared_modes::SeaStateFusion_CarmaSharedModes;
#elif defined(FILTER_TFG_SHARED)
#include "kalman_carma_shared_modes/SeaStateFusionFilter_TfgCarmaSharedModes.h"
using Fusion = ocean_imu::carma_shared_modes::SeaStateFusion_TfgCarmaSharedModes;
#elif defined(FILTER_TANGENT)
#include "kalman_carma_shared_modes/KalmanTangentCarmaSharedModes.h"
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"
using Fusion = ocean_imu::carma_shared_modes::SharedModesFusion<
    ocean_imu::carma_shared_modes::KalmanTangentCarmaSharedModes<float>>;
#else
#error "Select a fusion filter for the sanitizer replay"
#endif

extern const float g_std = 9.80665f;

template<class F> auto position(F& f) {
    if constexpr (requires { f.raw(); }) return f.raw().mekf().get_position();
    else return f.mekf().get_position();
}

int main() {
    using V = Eigen::Vector3f;
    constexpr float dt = 0.005f;
    Fusion filter;
    filter.begin(Fusion::Config{});
    int live_samples = 0;
    for (int k = 0; k < 12000; ++k) {
        const float t = k * dt;
        const float roll = 0.1f * std::sin(t);
        const Eigen::AngleAxisf to_body(-roll, V::UnitX());
        const V gyro(0.1f * std::cos(t), 0, 0);
        const V accel = to_body * V(0, 0, 0.4f * std::sin(t) - g_std);
        filter.update(dt, gyro, accel, 35.0f);
        if (k % 2 == 0) filter.updateMag(to_body * V(18, 1.5f, 43));
        if (filter.isLive()) {
            if (!position(filter).allFinite()) return 1;
            ++live_samples;
        }
    }
    std::printf("Live IMU samples: %d / 12000\n", live_samples);
    return live_samples >= 2000 ? 0 : 1;
}
