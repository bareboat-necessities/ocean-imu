# ALT finite runtime primitive status

## Scope

This note records the finite real-arithmetic runtime descriptors now composed on
PR #523. They are proof relations for the current shipping implementation, not
COMPLETE-BRMM/BIAS admission, contraction, startup capture or deployment
finite-precision certificates.

## Shared OU and accelerometer-bias roots

For one shipping step `h`, post-clamp time constant `tau`, and scalar decay
`alpha=exp(-h/tau)`, `finite_ou_runtime_primitives.py` generates the same mean
coefficients used by shipping:

- `phi_va=tau(1-alpha)`;
- the literal `|h/tau|<0.01` polynomial branch for `phi_pa,phi_Sa`;
- otherwise the closed form using `expm1(-h/tau)=alpha-1`.

That tuple feeds both the finite mean recurrence and `F_LL`; a caller at the
runtime entry cannot provide a detached second set.

For active BA, one `(tau_b,phi_b,Q_bacc)` root supplies both the mean factor and

`Q_BB = Q_bacc * tau_b/2 * (1-phi_b^2)`,

which is exactly the shipping `-tau_b/2*expm1(-2h/tau_b)` relation when
`phi_b=exp(-h/tau_b)`. Both BA cross-covariance blocks use the same `phi_b`.
H18 remains the exact held branch `phi_b=1,Q_BB=0`.

## Attitude/gyro-bias covariance runtime

`finite_attitude_runtime.py` materializes the shipping constant-rate attitude
covariance primitives for `with_gyro_bias=true`:

- the small-rate and trigonometric `R(w,t),B(w,t)` branches;
- `F_AA=[[R,B],[0,I]]`;
- structured `Q_AA` with isotropic-Qg fast integral or anisotropic Simpson
  `R Q R'`, Simpson `B Qbg B'`, and the closed `integral_B` cross term;
- the optional fast `Qbase*h` branch;
- the finite-valued 6x6 PSD-hygiene control branches.

Trig values and Eigen branch outcomes are explicit witnesses. They are not
silently asserted to be the deployed binary values.

## Integrated-OU Qaxis runtime

`finite_qaxis_runtime.py` materializes the literal
`IntegratedOUChain<T,3>::process_covariance` relation used by
`QdAxis4x1_analytic`:

- the small-`x` polynomial formulas;
- the general alpha-dependent formulas;
- the nested 3x3 v/p/a marginal regularization;
- S cross/marginal formulas;
- the final 4x4 `regularize_psd_if_needed` branch.

`finite_prediction_runtime.py` composes this with the attitude and OU/BA
runtime descriptors. Its highest prediction entry accepts no precomputed
shipping `F_AA,Q_AA,F_LL,Q_LL,Q_BB` or free decay factor. It also requires the
attitude covariance angular rate to equal the SAME bias-corrected gyro used by
the nominal quaternion prediction.

## Post-prediction and measurement control

`finite_post_prediction.py` materializes the literal post-prediction order:

1. optional pending `a_w` covariance synchronization by adding the positive
   spectral part of `target-P_aw` only to the `a_w` block;
2. covariance symmetry hygiene;
3. periodic S scheduler elapsed-credit update;
4. if due, S=0 service through the shipping safe-LDLT runtime branch.

The floor covers no-pending, eigensolver-failure and successful positive-part
branches. The scheduler covers due/not-due and overshoot remainder. A not-due S
edge is an identity suffix.

`finite_measurement_runtime.py` wraps the existing finite measurement algebra in
shipping's `safe_ldlt3_` control:

- first LDLT success: no innovation shift;
- first failure: exactly one bump
  `max(machine_epsilon,1e-6*(noise_scale+1))` and retry;
- second failure: reject without state/covariance update.

On an accepted retry the same shifted innovation is used by the inverse-free
gain and Joseph covariance relation.

## Remaining source/runtime obligations

The finite formulas above remove detached transition/process matrices but do not
yet supply the universal runtime history. Remaining work includes:

- bind exp/trig evaluations, matrix factorization/eigensolver outcomes and
  machine epsilon to rigorous deployment arithmetic;
- bind `tau,Sigma_aw,tau_b,Q_bacc`, a_w floor target, scheduler period,
  `R_acc,R_mag,R_S` and all relevant guards to the SAME frontend/tuner state;
- materialize WPE, adaptive bandpass, sigma statistic, candidate/active tuner
  state and staged commits;
- materialize asynchronous magnetometer continuation, all measurement/nonfinite
  rejection branches and every H18/A21 hybrid edge/prefix;
- attach corrected COMPLETE-BRMM and BIAS0/1/2 admission to that graph.

Until that complete source-uniform finite word is closed, storage feasibility,
rho, retained basin, startup capture and all ALT final gates remain fail-closed.
