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

    // Frequency-scale deadband (relative; 0 = every commit, as the research
    // filter).  The complexified mode's latent and imaginary states change
    // meaning with kappa; the CARMA3 schedule recommits kappa every 0.1 s, and
    // on small seas the tuner's frequency chatters.  A commit whose kappa is
    // within tol of the one in force keeps it; amplitudes always apply.
    void set_kappa_tolerance(T tol) { if (std::isfinite(tol) && tol >= T(0)) tol_ = tol; }
    T kappa_tolerance() const noexcept { return tol_; }
    // Amplitude deadband (relative; 0 = every commit).
    void set_amplitude_tolerance(T tol) { if (std::isfinite(tol) && tol >= T(0)) atol_ = tol; }
    bool set_vertical(T kappa, T accel_std) {
        return Base::set_vertical(hold_(kv_, kappa, tol_), hold_(sz_, accel_std, atol_));
    }
    bool set_horizontal(T kappa, T sx, T sy, T r) {
        return Base::set_horizontal(hold_(kh_, kappa, tol_), hold_(sx_, sx, atol_), hold_(sy_, sy, atol_), r);
    }

private:
    static T hold_(T& cur, T k, T tol) {
        if (tol > T(0) && cur > T(0) && std::isfinite(k) && std::fabs(k / cur - T(1)) < tol) return cur;
        cur = k;
        return k;
    }
    Estimator monitor_;
    T tol_ = T(0), atol_ = T(0), kv_ = T(0), kh_ = T(0), sz_ = T(0), sx_ = T(0), sy_ = T(0);
};

}  // namespace ocean_imu::carma_shared_modes
