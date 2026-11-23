#pragma once
// Copyright (c) 2026 Mikhail Grushinskiy. MIT License.
// TwoFrameGroup(3,6,2), 27 invariant errors. The common TFG core
// owns prediction, bias/noise correlations, sensor updates and group reset;
// the existing six-state motion policy owns the minimal CARMA(4,1) process.
#include "kalman_carma/Kalman3D_Wave_TFG_CARMA3.h"
#include "SharedModesEstimator.h"
#if defined(ARDUINO_ARCH_ESP32) && defined(__XTENSA__)
namespace ocean_imu::kalman {
// Keep the constructor's literal pool in a library object on Xtensa.
extern template Kalman3D_Wave_TFG_CARMA3<float,6,
    carma_shared_modes::SharedModesMotion<float>>::Kalman3D_Wave_TFG_CARMA3(
        const Eigen::Matrix<float,3,1>&, const Eigen::Matrix<float,3,1>&,
        const Eigen::Matrix<float,3,1>&, float,float,float,float,float,const carma::CarmaShape&);
}
#endif
namespace ocean_imu::carma_shared_modes {
template<typename T = float>
using KalmanTfgCarmaSharedModes = SharedModesEstimator<T,
    kalman::Kalman3D_Wave_TFG_CARMA3<T, 6, SharedModesMotion<T>>>;
#if defined(ARDUINO_ARCH_ESP32) && defined(__XTENSA__)
extern template void KalmanTfgCarmaSharedModes<float>::time_update(const Eigen::Vector3f&, float);
extern template void KalmanTfgCarmaSharedModes<float>::measurement_update_acc_only(const Eigen::Vector3f&, float);
#endif
}  // namespace ocean_imu::carma_shared_modes

#if defined(ARDUINO_ARCH_ESP32) && defined(__XTENSA__)
namespace ocean_imu::kalman {
extern template void carma_shared_modes::KalmanTfgCarmaSharedModes<float>::Base::time_update(
    const Eigen::Vector3f&, float);
extern template void carma_shared_modes::KalmanTfgCarmaSharedModes<float>::Base::measurement_update_acc_only(
    const Eigen::Vector3f&, float);
extern template void carma_shared_modes::KalmanTfgCarmaSharedModes<float>::Base::measurement_update_mag_only(
    const Eigen::Vector3f&);
}
#endif

namespace ocean_imu::carma {
struct TransportedMagCarmaConfig;
template<class Core> struct CarmaFusionConfigFor;
template<class T> struct CarmaFusionConfigFor<carma_shared_modes::KalmanTfgCarmaSharedModes<T>> {
    using type = TransportedMagCarmaConfig;
};
}  // namespace ocean_imu::carma
