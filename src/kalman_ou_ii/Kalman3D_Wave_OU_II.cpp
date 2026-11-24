// Keep the deployed OU-II prediction and its Xtensa literal pool in their own
// translation unit. The sketch declares the matching extern instantiation;
// the equations remain defined once in the public template header.
#if defined(ARDUINO) && !defined(EIGEN_NON_ARDUINO)
#ifndef EIGEN_STACK_ALLOCATION_LIMIT
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#endif
#include "AtomS3R/AtomS3R_FilterOptimization.h"
ATOMS3R_BEGIN_FILTER_MATH
#include "Kalman3D_Wave_OU_II.h"

template void ocean_imu::kalman::Kalman3D_Wave_OU_II<float, true, true>::time_update(
    const Eigen::Vector3f&, float);
ATOMS3R_END_FILTER_MATH
#endif
