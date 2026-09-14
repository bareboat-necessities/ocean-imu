# ALT finite runtime primitive status

## Scope

This note records the finite real-arithmetic runtime descriptors currently composed. They are proof relations for the current shipping implementation, not
COMPLETE-BRMM/BIAS admission, contraction, startup capture or deployment
finite-precision certificates.

## Shared OU and accelerometer-bias roots

For one shipping step `h`, post-clamp time constant `tau`, and scalar decay
`alpha=exp(-h/tau)`, `finite_ou_runtime_primitives.py` generates the same mean
coefficients used by shipping: `phi_va=tau(1-alpha)`, the literal
`|h/tau|<0.01` polynomial branch for `phi_pa,phi_Sa`, or the closed form using
`expm1(-h/tau)=alpha-1`. The same tuple feeds the finite mean recurrence and
`F_LL`; the runtime entry cannot provide a detached second set.

For active BA, one `(tau_b,phi_b,Q_bacc)` root supplies both the mean factor and
`Q_BB=Q_bacc*tau_b/2*(1-phi_b^2)`, exactly the shipping relation when
`phi_b=exp(-h/tau_b)`. Both BA cross-covariance blocks use the same `phi_b`.
H18 remains the exact held branch `phi_b=1,Q_BB=0`.

The held mean/error recurrence is important: holding `b_a_hat` does NOT mean
`e_ba` receives the physical driver alone. Since `e_ba=beta_true-b_a_hat`,

`e_ba+ = e_ba + (phi_true-1) beta_true + u_b`.

A regression in this PR initially asserted `e_ba+=e_ba+u_b`; focused CI caught
that test error. The finite predictor was already using the correct recurrence.
The corrected regression now checks invariance of `b_a_hat=beta_true-e_ba`.

## Attitude/gyro-bias covariance runtime

`finite_attitude_runtime.py` materializes the shipping constant-rate attitude
covariance primitives for `with_gyro_bias=true`: the small-rate and trigonometric
`R(w,t),B(w,t)` branches; `F_AA=[[R,B],[0,I]]`; structured `Q_AA` with isotropic
or Simpson gyro-noise integration, Simpson bias-noise integration and the closed
`integral_B` cross term; the optional fast `Qbase*h` branch; and finite-valued
6x6 PSD-hygiene control. Trig values and Eigen branch outcomes remain explicit
runtime witnesses, not deployment certificates.

## Integrated-OU Qaxis runtime

`finite_qaxis_runtime.py` materializes the literal
`IntegratedOUChain<T,3>::process_covariance` relation used by
`QdAxis4x1_analytic`: small-`x` polynomial and general alpha-dependent formulas,
the nested 3x3 v/p/a marginal regularization, S cross/marginal terms, and final
4x4 `regularize_psd_if_needed` branch.

`finite_prediction_runtime.py` composes this with attitude and OU/BA runtime
descriptors. Its highest prediction entry accepts no precomputed shipping
`F_AA,Q_AA,F_LL,Q_LL,Q_BB` or free decay factor. It additionally enforces that
the attitude covariance angular rate equals the SAME bias-corrected gyro used by
the nominal quaternion prediction.

## Post-prediction and measurement control

`finite_post_prediction.py` materializes optional pending `a_w` covariance
positive-part synchronization, covariance symmetry hygiene, periodic S scheduler
elapsed credit, and due S=0 service. No-pending, eigensolver-failure and success
branches are explicit; not-due S is an identity suffix.

`finite_measurement_runtime.py` wraps the finite measurement algebra in shipping
`safe_ldlt3_`: first success, one bump
`max(machine_epsilon,1e-6*(noise_scale+1))` and retry, or second failure with no
state/covariance update. On an accepted retry the same shifted innovation is used
by gain and Joseph arithmetic.

## Asynchronous magnetic product state

`finite_live_interleave.py` joins the exact startup bridge to successive Live
IMU, magnetic and external-hold events. `finite_live_magnetic_word.py` binds a
fixed qualified physical field/hard-iron model, raw continuous accumulation,
refinement, reference/yaw writes, hard-iron/reference application and the same
corrected packet's measurement/count suffix. Full covariance and all frontend /
calibration/scheduler memory persist between event types.

`finite_continuous_mag_runtime.py` derives fit and reference from one moment
state; failed due solves replace the estimate, and unsuccessful applications
retain the shipping anchor/clock mutations. The finite-real bound and precise
arithmetic/source limitations are in `ou3-alt-live-magnetic-word.md`.

## Remaining source/runtime obligations

The finite formulas above remove detached transition/process matrices but do not
yet supply the universal runtime history. Remaining work includes binding
exp/trig evaluations and numerical factorization/eigensolver outcomes; binding
`tau,Sigma_aw,tau_b,Q_bacc`, a_w floor target, scheduler period,
`R_acc,R_mag,R_S` and guards to one frontend/tuner predecessor; qualifying the composed
WPE/bandpass/sigma/tuner candidate/active staged commits and asynchronous
magnetometer continuation, remaining reject/nonfinite branches and H18/A21
source conditions; and attaching corrected COMPLETE-BRMM plus BIAS0/1/2 admission.

Until that complete source-uniform finite word is closed, storage feasibility,
rho, retained basin, startup capture and all ALT final gates remain fail-closed.
