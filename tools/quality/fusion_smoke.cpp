// Bounded host replay for the ASan/UBSan gate. Compile one production wrapper per job.
//
// 600 s of one long-crested sea travelling toward 30 degrees, with exact
// sensor kinematics: horizontal and vertical wave accelerations in quadrature,
// roll and pitch at the wave frequency, constant heading, 200 Hz IMU and
// 100 Hz magnetometer. The duration outlasts the shared-mode coupling warm-up
// (120 s) and monitor horizon (60 s), so for the shared-mode cores the coupling
// confidence, band regression and excitation-cap paths also run under the
// sanitizers, and the replay fails unless that coupling engages.
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

template<class F> auto& core(F& f) {
    if constexpr (requires { f.raw(); }) return f.raw().mekf();
    else return f.mekf();
}

// Shared-mode cores must have engaged their estimated coupling by the end.
template<class F> bool coupling_engaged(F& f) {
    const auto& m = core(f);
    if constexpr (requires { m.couplingConfidence(); }) {
        const auto st = m.coupling();
        const float applied = st.applied_coherence_xz + st.applied_coherence_yz;
        std::printf("Coupling confidence %.3f, applied coherence %.3f\n",
                    static_cast<double>(m.couplingConfidence()), static_cast<double>(applied));
        return m.couplingConfidence() > 0.5f && applied > 0.1f;
    } else {
        return true;
    }
}

int main() {
    using V = Eigen::Vector3d;
    using Eigen::AngleAxisd;
    constexpr int kSamples = 600 * 200;
    constexpr double dt = 0.005;
    constexpr double pi = 3.14159265358979323846;
    constexpr double w = 2.0 * pi / 7.0;      // 7 s wave period
    constexpr double heave = 0.5;             // m amplitude
    constexpr double to_rad = 30.0 * pi / 180.0;
    constexpr double yaw = 20.0 * pi / 180.0;
    const V travel(std::cos(to_rad), std::sin(to_rad), 0.0);
    const V mag_world(18.0, 1.5, 43.0);       // uT, NED
    Fusion filter;
    filter.begin(Fusion::Config{});
    int live_samples = 0;
    for (int k = 0; k < kSamples; ++k) {
        const double t = k * dt, s = std::sin(w * t), c = std::cos(w * t);
        // Linear wave, NED: up-displacement heave*cos(wt), so the downward
        // acceleration is heave w^2 cos(wt); the horizontal displacement
        // heave*sin(wt) along the travel direction leads it by a quarter period.
        const V a_world = -heave * w * w * s * travel + V(0.0, 0.0, heave * w * w * c);
        const double roll = 0.06 * std::sin(w * t + 0.4), roll_rate = 0.06 * w * std::cos(w * t + 0.4);
        const double pitch = 0.04 * std::sin(w * t - 0.7), pitch_rate = 0.04 * w * std::cos(w * t - 0.7);
        const Eigen::Matrix3d R_bw = (AngleAxisd(yaw, V::UnitZ()) * AngleAxisd(pitch, V::UnitY()) *
                                      AngleAxisd(roll, V::UnitX())).toRotationMatrix();
        // Body rates of R_bw = Rz(yaw) Ry(pitch) Rx(roll) at constant yaw.
        const V gyro(roll_rate, pitch_rate * std::cos(roll), -pitch_rate * std::sin(roll));
        const V accel = R_bw.transpose() * (a_world - V(0.0, 0.0, g_std));
        filter.update(static_cast<float>(dt), gyro.cast<float>(), accel.cast<float>(), 35.0f);
        if (k % 2 == 0) filter.updateMag((R_bw.transpose() * mag_world).cast<float>());
        if (filter.isLive()) {
            if (!core(filter).get_position().allFinite()) return 1;
            ++live_samples;
        }
    }
    std::printf("Live IMU samples: %d / %d\n", live_samples, kSamples);
    const bool coupled = coupling_engaged(filter);
    return live_samples >= kSamples / 2 && coupled ? 0 : 1;
}
