# TFG mathematical hardening

This note records the corrections made after review of the two-frame Lie-group filter. It is deliberately narrower than `tfg-design.md`: it states what is now guaranteed by code/tests and what is **not** claimed.

## Estimator classification

The implementation is a **right-invariant Lie-group error-state EKF with invariant-style measurements** using the convention

- `R = R_bw` (body to world),
- right-invariant error `eta = Xhat o X^-1`,
- left correction `Xhat+ = Exp_G(-K r) o Xhat-`.

The complete marine estimator is **not claimed to be an autonomous/log-linear InEKF**. The process linearization contains state-dependent world/bias couplings (for example terms built from `-[x]x Bmap`), so the strongest group-affine Barrau--Bonnabel log-linear property is not assumed by the implementation or by its tests.

## 1. Finite retraction covariance reset

The measurement update first forms the ordinary Joseph posterior covariance in the pre-injection tangent coordinates. For a finite estimated error

`a = K r`,

the state is then corrected with `Exp_G(-a)`. The covariance is transported by the full retraction Jacobian

`J_reset = d/d(delta) Log_G(Exp_G(-a) Exp_G(a + delta)) |_(delta=0)`

and the stored covariance becomes

`P+ = J_reset P_Joseph J_reset^T`.

For Level 2 this is the right Jacobian of the complete two-frame group, including attitude/world, gyro-bias and accelerometer-bias blocks. Level 1 keeps its additive body-bias blocks Euclidean while using the same group reset for attitude/world coordinates.

`covariance_transport-test.cpp` independently finite-differences the defining map above and compares both the analytic Jacobian and the final covariance against that oracle. The previous assertion that Joseph covariance must remain unchanged after injection has been removed.

## 2. Accelerometer-bias model parity with OU-III

The residual accelerometer bias now uses the same slow first-order Gauss--Markov/OU model as current OU-III:

`beta_a(h) = exp(-h/tau_ba) beta_a(0)`

with exact discrete driving covariance

`Q_ba,d = (tau_ba/2) (1 - exp(-2h/tau_ba)) Q_ba,c`.

The default correlation time is 5000 s. `set_Q_bacc_rw()` retains its historical name but, as in OU-III, now accepts a continuous driving-noise **standard deviation per sqrt(s)** and squares it internally. `set_acc_bias_ou_stationary_std()` is also available.

For Level 2 the bias-error transition is

`Phi_beta_a,beta_a = exp(-h/tau_ba) Exp(omega_world h)`

and the discrete process covariance is rotated into the end-of-step world frame. When startup policy freezes accelerometer-bias updates, both the physical body-bias mean and its process noise are frozen; Level-2 tangent coordinates still transform consistently with the changing attitude.

## 3. Magnetic first lock is a gauge transform

The first magnetic lock is treated as choosing the horizontal world-frame gauge, not as a measurement correction. `apply_world_yaw_gauge()` rotates coherently:

- `R`,
- all world columns `v, p, S, a_w`,
- world-referred stationary/tuning covariances,
- an existing magnetic world reference, and
- the corresponding tangent/covariance blocks (including Level-2 bias coordinates).

Body-frame bias states themselves are not rotated. `SeaStateFusionFilter_TFG::updateMag()` now uses this operation instead of assigning only `state().R`.

## 4. Gyro-bias stochastic discretization during rotation

The old `h^2/2` and `h^3/3` gyro-bias terms are exact only when rotation is effectively constant in the world frame. With finite angular rate, a body-frame gyro-bias noise impulse at time `s` contributes to attitude through

`C(s) = integral_s^h R(t) dt`.

The implementation now evaluates

- `Q_phi,phi^(b) = integral C(s) Q_b C(s)^T ds`,
- `Q_phi,beta_g = -integral C(s) Q_b B_final^T ds`,

with five-point Gauss--Legendre quadrature. The final gyro-bias block is still analytic. A separate high-rate midpoint integration in `covariance_transport-test.cpp` is the stochastic oracle, so the production quadrature is not tested against itself.

The gyro-white-noise leakage into world-vector blocks remains the documented first-order approximation inherited from the surrounding estimator; this hardening does not claim that the entire `Q_d` is a closed-form exact discretization.

## Validation scope

The branch adds `.github/workflows/tfg-validation.yml`, which builds and runs the existing group, convention, OU-chain, propagation, Jacobian, covariance and orchestrator tests on TFG changes. The covariance test additionally checks:

1. full-group reset Jacobian against finite differences;
2. post-update covariance against the independently transported Joseph covariance;
3. OU accelerometer-bias mean/covariance against the analytic solution;
4. rotating gyro-bias `Q_d` against a high-rate numerical oracle;
5. complete world-yaw gauge congruence.

These changes remove the four mathematical/implementation gaps identified in review. They do **not** by themselves establish that TFG supersedes OU-III empirically. The paired wave simulations should be rerun after the corrected dynamics before changing performance or replacement claims.
