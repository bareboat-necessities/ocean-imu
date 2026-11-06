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

    // Amplitude commit deadband (relative; 0 = every commit, the default and
    // the research filter's behaviour).  Screened on seeds 1001-1020 and
    // rejected (docs/kalman-carma-shared-modes.md section 2); kept so that
    // screening is reproducible.  A commit whose amplitude is within tol of
    // the one in force keeps it; kappa always applies.
    void set_amplitude_tolerance(T tol) { if (std::isfinite(tol) && tol >= T(0)) atol_ = tol; }
    bool set_vertical(T kappa, T accel_std) {
        return Base::set_vertical(kappa, hold_(sz_, accel_std, atol_));
    }
    bool set_horizontal(T kappa, T sx, T sy, T r) {
        return Base::set_horizontal(kappa, hold_(sx_, sx, atol_), hold_(sy_, sy, atol_), r);
    }

private:
    static T hold_(T& cur, T k, T tol) {
        if (tol > T(0) && cur > T(0) && std::isfinite(k) && std::fabs(k / cur - T(1)) < tol) return cur;
        cur = k;
        return k;
    }
    Estimator monitor_;
    T atol_ = T(0), sz_ = T(0), sx_ = T(0), sy_ = T(0);
};

}  // namespace ocean_imu::carma_shared_modes
