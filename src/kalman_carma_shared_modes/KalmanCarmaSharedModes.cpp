/*
    Copyright (c) 2026 Mikhail Grushinskiy
    Released under the MIT License
*/

#if defined(ARDUINO_ARCH_ESP32) && defined(__XTENSA__)
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
}  // namespace ocean_imu::kalman
#endif
