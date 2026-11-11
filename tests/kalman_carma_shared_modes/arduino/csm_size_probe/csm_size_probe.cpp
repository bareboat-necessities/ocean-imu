/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    ESP32-S3 (AtomS3R) size and stack probe of the EXPERIMENTAL standalone
    shared-modes QMEKF (src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h,
    27 states) and, for reference, the CARMA3 baseline wrapper A41.  Compiled
    by build.sh with Espressif's Xtensa GCC and arduino-esp32's C++ flags; no
    sketch, no hardware.  Variants: 0 A41, 1 CSM.
*/

#if defined(__XTENSA__) || defined(CSM_PROBE_HOST)

#define EIGEN_NON_ARDUINO
#include <type_traits>
#include <utility>

#include "kalman_carma/SeaStateFusionFilter_CARMA3.h"
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"
#include "kalman_carma_shared_modes/SeaStateFusionFilter_TfgCarmaSharedModes.h"
#include "kalman_carma/SeaStateFusionFilter_TFG_CARMA3.h"

const float g_std = 9.80665f;

#if CSM_PROBE_VARIANT == 0
using Top = ocean_imu::carma::SeaStateFusion_CARMA3<TrackerType::KALMANF, 4, false, 4, false>;
using Wrapper = Top;
#elif CSM_PROBE_VARIANT == 1
using Top = ocean_imu::carma_shared_modes::SeaStateFusion_CarmaSharedModes;
using Wrapper = Top;
#elif CSM_PROBE_VARIANT == 2
using Top = ocean_imu::carma::SeaStateFusion_TFG_CARMA3<TrackerType::KALMANF,4>;
using Wrapper = Top;
#elif CSM_PROBE_VARIANT == 3
using Top = ocean_imu::carma_shared_modes::SeaStateFusion_TfgCarmaSharedModes;
using Wrapper = Top;
#endif
using Mekf = typename std::remove_reference<decltype(std::declval<Wrapper&>().raw().mekf())>::type;

extern "C" {
extern const unsigned csm_probe_sizeof_top;
extern const unsigned csm_probe_sizeof_wrapper;
extern const unsigned csm_probe_sizeof_mekf;
const unsigned csm_probe_sizeof_top = sizeof(Top);
const unsigned csm_probe_sizeof_wrapper = sizeof(Wrapper);
const unsigned csm_probe_sizeof_mekf = sizeof(Mekf);
}

extern "C" Top* csm_probe_make() {
    auto* f = new Top();
    typename Top::Config cfg;
    cfg.sigma_a = Eigen::Vector3f::Constant(0.0294f);
    cfg.sigma_g = Eigen::Vector3f::Constant(3.1e-5f);
    cfg.sigma_m = Eigen::Vector3f::Constant(7.7f);
    f->begin(cfg);
    return f;
}

extern "C" float csm_probe_step(Top* f, float dt, const float* gyr, const float* acc, const float* mag) {
    f->update(dt, Eigen::Vector3f(gyr[0], gyr[1], gyr[2]), Eigen::Vector3f(acc[0], acc[1], acc[2]));
    if (mag) f->updateMag(Eigen::Vector3f(mag[0], mag[1], mag[2]));
#if CSM_PROBE_VARIANT == 0 || CSM_PROBE_VARIANT == 2
    return f->raw().mekf().get_position().z();
#else
    return f->raw().mekf().get_position().z() + f->raw().mekf().couplingConfidence();
#endif
}

#endif
