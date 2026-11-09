#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy
    Released under the MIT License

    KalmanCarmaSharedModes<T> -- EXPERIMENTAL QMEKF marine-motion filter with
    shared directional stochastic modes (Option C of the coupled-prior study;
    docs/kalman-carma-shared-modes.md, doc/kalman_carma_shared_modes/).  Not
    a deployed estimator and not a default: SeaStateFusionFilter_OU_III stays
    the shipping filter and the independent-axis CARMA3 QMEKF
    (kalman_carma/Kalman3D_Wave_CARMA3.h) is unchanged.

    The class IS the CARMA3 QMEKF core Kalman3D_Wave_CARMA3<T, 4, false, 4,
    false, true, true, SharedModesMotion<T>> (attitude error, gyroscope bias,
    accelerometer bias, Joseph updates, multiplicative attitude reset, the
    shared-mode motion block: one complexified CARMA(4,1) mode per axis in its
    minimal realization, six real states each (CarmaSharedModesMotion.h),
    axes sharing its dynamics and driven by circular complex noise with a
    Hermitian excitation estimated from the filter's own world-frame
    acceleration) -- 27 error states -- extended with the coupling
    confidence (CarmaCouplingConfidence.h):

      time_update(gyr, Ts)                 base prediction; remembers the rate
      measurement_update_acc_only(acc, T)  base update, then one confidence
                                           step; at 0.1 s the applied coupling
                                           gain becomes c and, when c falls,
                                           the unexplained horizontal-vertical
                                           cross-covariance is shrunk by
                                           c_new / c_old (PSD-preserving)

    so it drops into any schedule that drives a Kalman3D_Wave_CARMA3, the
    CARMA3 wrapper in particular (SeaStateFusionFilter_CarmaSharedModes.h).
    At c = 0 the filter is the independent-axis CARMA3 baseline exactly, in
    law and through every model commit.

    Scalar T: float on the target (the default), double for analysis.  Fixed
    size throughout; nothing is allocated after construction.  Ground truth,
    oracle coupling and RAO knowledge never enter the class;
    setExternalCoupling() exists for research priors and the robustness
    tests, and the confidence checks it against the IMU like the estimated one.
*/

#include <algorithm>
#include <cmath>
#include <complex>

#include "kalman_carma/Kalman3D_Wave_CARMA3.h"
#include "kalman_carma_shared_modes/CarmaSharedModesMotion.h"
#include "kalman_carma_shared_modes/SharedModesEstimator.h"

#if defined(ARDUINO_ARCH_ESP32) && defined(__XTENSA__)
namespace ocean_imu::kalman {
// Emit the target's constructor in its own library object. In the large
// Arduino sketch object its literal pool can land beyond Xtensa L32R reach.
// Other scalar types and motion configurations remain header-only.
extern template Kalman3D_Wave_CARMA3<
    float, 4, false, 4, false, true, true,
    carma_shared_modes::SharedModesMotion<float>>::Kalman3D_Wave_CARMA3(
        const Eigen::Matrix<float,3,1>&, const Eigen::Matrix<float,3,1>&,
        const Eigen::Matrix<float,3,1>&, float, float, float, float, float,
        const carma::CarmaShape&, const carma::CarmaShape&);
}  // namespace ocean_imu::kalman
#endif

namespace ocean_imu::carma_shared_modes {
template<typename T = float>
using KalmanCarmaSharedModes = SharedModesEstimator<T,
    kalman::Kalman3D_Wave_CARMA3<T, 4, false, 4, false, true, true, SharedModesMotion<T>>>;
}  // namespace ocean_imu::carma_shared_modes
