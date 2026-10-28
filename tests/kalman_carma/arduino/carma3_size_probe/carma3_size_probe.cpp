/*
  carma3_size_probe -- ESP32-S3 code-size and stack probe of the EXPERIMENTAL
  three-axis CARMA filters against the shipping OU-III (docs/carma-three-axis.md).

  Not a sketch and not a link: one translation unit per variant
  (CARMA3_PROBE_VARIANT), compiled with the Espressif Xtensa GCC for the
  ESP32-S3 at -Os, instantiating the complete fusion wrapper's sample path
  (update + magnetometer update) so every template the firmware would need is
  emitted.  build.sh reports the object's code and data sizes, the wrapper
  and filter object sizes (sizeof on the target) and the largest stack frame
  on the update path.  No hardware is involved: these are build-time numbers,
  not on-target timing.

    0  SeaStateFusion_OU_III                       (shipping baseline)
    1  SeaStateFusion_CARMA3<KALMANF, 4,F, 4,F>    CARMA(4,1) x/y/z
    2  SeaStateFusion_CARMA3<KALMANF, 5,F, 5,F>    CARMA(5,2) x/y/z
    3  SeaStateFusion_CARMA3<KALMANF, 4,F, 5,F>    CARMA(4,1) x/y + CARMA(5,2) z
    4  SeaStateFusion_CARMA3<KALMANF, 5,T, 5,F>    CARMA(5,2) x/y/z + S_x, S_y
    5  SeaStateFusion_CARMA3<KALMANF, 5,F, 4,F>    CARMA(5,2) x/y + CARMA(4,1) z
    6  SeaStateFusion_CARMA3<KALMANF, 4,T, 4,F>    CARMA(4,1) x/y/z + S_x, S_y
    7  SeaStateFusion_TFG_CARMA3<KALMANF, 4>       CARMA(4,1) x/y/z on the two-frame
                                                   Lie group (docs/tfg-carma3-research.md)
*/

// Built only by the Xtensa target compiler or by `build.sh host`; host
// compilers otherwise (the repository's warnings and clang-tidy gates compile
// every tests/*.cpp) see an empty translation unit.
#if defined(__XTENSA__) || defined(CARMA3_PROBE_HOST)

#define EIGEN_NON_ARDUINO
#define EIGEN_NO_DEBUG

#include <Eigen/Dense>

#include "kalman_ou_iii/SeaStateFusionFilter_OU_III.h"
#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"
#include "kalman_carma/SeaStateFusionFilter_TFG_CARMA3.h"

#if CARMA3_PROBE_VARIANT == 0
using Fusion = SeaStateFusion_OU_III<TrackerType::KALMANF>;
#elif CARMA3_PROBE_VARIANT == 1
using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false>;
#elif CARMA3_PROBE_VARIANT == 2
using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 5, false, 5, false>;
#elif CARMA3_PROBE_VARIANT == 3
using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 5, false>;
#elif CARMA3_PROBE_VARIANT == 4
using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 5, true, 5, false>;
#elif CARMA3_PROBE_VARIANT == 5
using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 5, false, 4, false>;
#elif CARMA3_PROBE_VARIANT == 6
using Fusion = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, true, 4, false>;
#else
using Fusion = ocean_imu::carma::SeaStateFusion_TFG_CARMA3<TrackerType::KALMANF, 4>;
#endif

const float g_std = 9.80665f;

// Read back by build.sh from the object (extern: external linkage, so each
// constant gets its own .rodata.<name> section).
extern "C" {
extern const unsigned carma3_probe_sizeof_fusion;
extern const unsigned carma3_probe_sizeof_mekf;
const unsigned carma3_probe_sizeof_fusion = sizeof(Fusion);
const unsigned carma3_probe_sizeof_mekf = sizeof(typename std::remove_reference<
    decltype(std::declval<Fusion&>().raw().mekf())>::type);
}

extern "C" Fusion* carma3_probe_make() {
    auto* f = new Fusion();
    typename Fusion::Config cfg;
    cfg.sigma_a = Eigen::Vector3f::Constant(0.0294f);
    cfg.sigma_g = Eigen::Vector3f::Constant(3.1e-5f);
    cfg.sigma_m = Eigen::Vector3f::Constant(7.7f);
    f->begin(cfg);
    return f;
}

extern "C" float carma3_probe_step(Fusion* f, float dt, const float* gyr, const float* acc, const float* mag) {
    f->update(dt, Eigen::Vector3f(gyr[0], gyr[1], gyr[2]), Eigen::Vector3f(acc[0], acc[1], acc[2]));
    if (mag) f->updateMag(Eigen::Vector3f(mag[0], mag[1], mag[2]));
    return f->raw().mekf().get_position().z();
}

#endif  // __XTENSA__ || CARMA3_PROBE_HOST
