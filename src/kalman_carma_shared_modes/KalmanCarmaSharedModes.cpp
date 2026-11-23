/*
    Copyright (c) 2026 Mikhail Grushinskiy
    Released under the MIT License
*/

#if defined(ARDUINO_ARCH_ESP32) && defined(__XTENSA__)
// Compile the 200 Hz estimator at speed optimization in its own object.
// Calibration, graphics and the rest of the Arduino application keep -Os.
#if defined(__GNUC__) && !defined(__clang__) && defined(__OPTIMIZE_SIZE__)
#pragma GCC optimize ("O2", "unroll-loops")
#endif
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"
#include "kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h"

namespace ocean_imu::kalman {
template Kalman3D_Wave_CARMA3<
    float, 4, false, 4, false, true, true,
    carma_shared_modes::SharedModesMotion<float>>::Kalman3D_Wave_CARMA3(
        const Eigen::Matrix<float,3,1>&, const Eigen::Matrix<float,3,1>&,
        const Eigen::Matrix<float,3,1>&, float, float, float, float, float,
        const carma::CarmaShape&, const carma::CarmaShape&);
template Kalman3D_Wave_TFG_CARMA3<float,6,
    carma_shared_modes::SharedModesMotion<float>>::Kalman3D_Wave_TFG_CARMA3(
        const Eigen::Matrix<float,3,1>&, const Eigen::Matrix<float,3,1>&,
        const Eigen::Matrix<float,3,1>&, float,float,float,float,float,const carma::CarmaShape&);
template void carma_shared_modes::KalmanTfgCarmaSharedModes<float>::Base::time_update(
    const Eigen::Vector3f&, float);
template void carma_shared_modes::KalmanTfgCarmaSharedModes<float>::Base::measurement_update_acc_only(
    const Eigen::Vector3f&, float);
template void carma_shared_modes::KalmanTfgCarmaSharedModes<float>::Base::measurement_update_mag_only(
    const Eigen::Vector3f&);
}  // namespace ocean_imu::kalman

namespace ocean_imu::carma_shared_modes {
template void KalmanTfgCarmaSharedModes<float>::time_update(const Eigen::Vector3f&, float);
template void KalmanTfgCarmaSharedModes<float>::measurement_update_acc_only(const Eigen::Vector3f&, float);
}
#endif
