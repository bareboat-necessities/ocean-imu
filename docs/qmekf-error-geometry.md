# Quaternion MEKF error geometry

The right attitude error is `q_true = q_ref * Exp(delta_theta)`. For constant
bias-corrected body rate `w`, its transition is `Exp(-[w]x dt)`, not the
three-by-three block of the four-dimensional half-angle quaternion propagator.
The gyro-bias cross block is `-dt * J_r(w dt)`. After injecting `delta_theta`,
the covariance reset uses `J_r(delta_theta)` and transports the attitude/bias
cross blocks too. These corrections remain in `src/ahrs/KalmanQMEKF.h`.

`tests/ahrs/qmekf_geometry-test.cpp` checks the transition against finite
differences and quadrature, checks the injection reset derivative, and checks
axial covariance retention. This is an algebraic/numerical implementation
check, not a global nonlinear capture or dynamic-acceleration qualification.
No motion detector, covariance erasure or automatic bias repair follows from
these identities.
