#pragma once

/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy
    Released under the MIT License

    SharedModesMotion -- the motion block of the EXPERIMENTAL shared-modes
    QMEKF (KalmanCarmaSharedModes.h).

    It IS the validated Option C block of the coupled-prior study
    (kalman_carma/CarmaCoupledMotion.h, CarmaSharedModeMotion =
    CarmaKroneckerMotion<T, true>: one complexified CARMA(4,1) mode per axis,
    eight real states each, axes sharing its dynamics and driven by circular
    complex noise with the Hermitian excitation of
    kalman_carma/CarmaCrossSpectrumEstimator.h), unchanged: same dynamics,
    discretization, excitation mapping and readouts.  It adds only a second,
    short-horizon band cross-spectrum estimate (the MONITOR) fed with the same
    world-frame acceleration the coupling estimator sees, which the coupling
    confidence (CarmaCouplingConfidence.h) compares with the prior.  The
    monitor never enters the prior, so a filter whose confidence is held at
    one reproduces the research filter bit for bit.
*/

#include "kalman_carma/CarmaCoupledMotion.h"

namespace ocean_imu::carma_shared_modes {

template<typename T>
class SharedModesMotion : public kalman::carma::CarmaKroneckerMotion<T, true> {
public:
    using Base = kalman::carma::CarmaKroneckerMotion<T, true>;
    using Vector3 = typename Base::Vector3;
    using Estimator = kalman::carma::CarmaCrossSpectrumEstimator<T>;
    using Base::Base;

    // Same call as the base (the QMEKF's accelerometer update), plus the
    // monitor.
    void observe(const Vector3& a_world) {
        monitor_.update(a_world, this->step(), this->operating_frequency_hz());
        Base::observe(a_world);
    }

    Estimator& monitor() noexcept { return monitor_; }
    const Estimator& monitor() const noexcept { return monitor_; }

private:
    Estimator monitor_;
};

}  // namespace ocean_imu::carma_shared_modes
