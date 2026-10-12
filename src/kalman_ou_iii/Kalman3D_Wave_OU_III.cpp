#if defined(ARDUINO_ARCH_ESP32) && defined(__XTENSA__)
#include "Kalman3D_Wave_OU_III.h"

// Arduino compiles library sources separately from the sketch, keeping the
// prediction's Xtensa literal pool close to the instructions that use it.
template void ocean_imu::kalman::Kalman3D_Wave_OU_III<float, true, true>::time_update(
    const Eigen::Matrix<float, 3, 1>&, float);
#endif
