#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy
    Released under the MIT License

    KalmanCarmaSharedModesMinimal<T, R> -- EXPERIMENTAL 27-state realization
    of the shared-modes QMEKF (docs/kalman-carma-shared-modes-minimal.md,
    doc/kalman_carma_shared_modes/kalman_carma_shared_modes-minimal.tex).

    It is the shared-modes QMEKF core KalmanCarmaSharedModes<T, Motion>
    (KalmanCarmaSharedModes.h: the CARMA3 QMEKF with attitude error,
    gyroscope and accelerometer bias, Joseph updates, multiplicative reset,
    magnetometer, the coupling estimate, confidence and PSD-preserving
    fallback) with the six-state-per-axis motion block
    CarmaMinimalSharedModeMotion<T, R> (CarmaSharedModesMinimalMotion.h)
    instead of the eight-state complexified one:

        [dtheta(3), b_g(3), x(6), y(6), z(6), b_a(3)]   27 error states.

    Every per-sample call, configuration, diagnostic and the fallback are the
    core's, unchanged; only the motion block differs.  At fixed kappa the
    filter's estimates equal the 33-state filter's (the removed directions
    are unobservable); R selects what a model commit holds fixed:
      Cascade (default)  the CARMA3 baseline's [p, v, a, l] plus the
                         normalized lossless colouring (coupling off: the
                         CARMA3 baseline A41 in law, adaptation included);
      Quotient           the 33-state block's observable quotient with a
                         regression lift at kappa commits (amplitude and
                         coupling changes: the 33-state filter exactly);
      QuotientHold       the quotient without the lift (ablation).

    The 33-state filter (KalmanCarmaSharedModes<T>) is unchanged and remains
    the shared-modes reference; this is a research alternative, not a
    default.
*/

#include "kalman_carma_shared_modes/CarmaSharedModesMinimalMotion.h"
#include "kalman_carma_shared_modes/CarmaSharedModesMotion.h"
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"

namespace ocean_imu::carma_shared_modes {

template<typename T = float, MinimalRealization R = MinimalRealization::Cascade>
using MinimalSharedModesMotion = SharedModesMotion<T, CarmaMinimalSharedModeMotion<T, R>>;

template<typename T = float, MinimalRealization R = MinimalRealization::Cascade>
using KalmanCarmaSharedModesMinimal = KalmanCarmaSharedModes<T, MinimalSharedModesMotion<T, R>>;

}  // namespace ocean_imu::carma_shared_modes
