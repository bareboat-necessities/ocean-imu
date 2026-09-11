# OU-III parallel ALT proof handover

## Resume point and independent scope

Read `AGENTS.md`, `docs/ou3-alt-proof-plan.md`, this handover,
`docs/ou3-alt-finite-measurement-proof.md`, `docs/ou3-alt-finite-core-composition.md`,
`docs/ou3-alt-runtime-primitives.md`, the ALT section of
`docs/ou3-proof-research-state.md`, and `docs/ou3-brmm-main-handover.md`.
Continue PR #523 on its branch. The original P2/P3/P4/P5 route remains
independently continuable and unchanged; `P3=1e-18` remains frozen.

## Immutable physical contracts

Retain joint24 `z=(c,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta)`, full 21-state covariance,
all motion/bias cross terms, one persistent Live S origin and the corrected
COMPLETE-BRMM physical envelope including acceleration <=8.8 m/s^2 and all-time
centered primitive `D_S<=1100 m*s`. Preserve same-signal WPE/bandpass/sigma/tau/
T_S ancestry, staged tuner state and separate BIAS0/1/2 histories. Unknown state
errors may not be renamed as bounded disturbances.

## Finite graph now represented

The supplying finite identities retain actual continuous physical attitude
increments and angular/model defect, correlated physical translation moments,
gyro-bias drift, one shared physical accelerometer-bias driver, finite
accelerometer/magnetometer/S residual secants, inverse-free `K Sigma=N`, both
quaternion injection branches, same-beta radial projection and full Joseph/reset
covariance. H18 retains latent BA covariance; no 18-state marginal substitute is
used.

`finite_core.py` recursively carries finite joint24 mean/error, nominal attitude,
full covariance and physical predecessor through represented prefixes.

The prediction runtime is composed through:

- `finite_prediction_covariance.py`: full AA/linear/BA covariance congruence;
- `finite_ou_runtime_primitives.py`: one OU decay root for mean/F_LL and one BA
  root for mean, BA cross factors and Q_BB; H18 held branch exact;
- `finite_attitude_runtime.py`: constant-rate R/B, structured/fast attitude Q,
  isotropic/Simpson branches, integral-B cross term and 6x6 PSD hygiene;
- `finite_qaxis_runtime.py`: literal small/general IntegratedOUChain<3> Qaxis
  formulas, nested 3x3 marginal hygiene and final 4x4 hygiene;
- `finite_prediction_runtime.py`: highest prediction entry, accepting no
  precomputed shipping transition/process matrices and requiring the SAME
  bias-corrected gyro for nominal and covariance attitude propagation.

Post-prediction and measurement control includes:

- `finite_post_prediction.py`: pending a_w positive-part covariance sync,
  symmetry hygiene, scheduler credit/due recurrence and due S service;
- `finite_measurement_runtime.py`: safe-LDLT first success, one deterministic
  bump/retry and double-failure rejection; retry uses the same shifted innovation
  in gain/Joseph, rejection preserves state/covariance.

`finite_tuner_commit.py` now materializes the next-sample staged commit boundary.
One pending `TuneState` supplies the OU tau, stationary Sigma_aw, realized
pseudo-update period, and—when Live—the anisotropic applied R_S. The period uses
the same applied tau. The cubic cadence normalization is represented by an
explicit same-period square-root witness. Periodic online adaptation uses
`apply_ou_tune_(false)`, so it changes the stationary process covariance but does
not queue posterior a_w synchronization; discrete sync remains a separate branch.

A focused-CI regression initially expected held H18 `e_ba += u_b`; that was a
test error, not a predictor error. The correct held estimator relation is
`e_ba+ = e_ba + (phi_true-1) beta + u_b`, equivalently constant
`b_a_hat=beta-e_ba`. The regression was corrected to assert that invariant.

## What remains open

**The complete source-uniform finite 600-step runtime word is still NOT
materialized.** Current limiter remains C/E same-history source/runtime
attachment, not demonstrated instability or failed common storage.

Remaining major attachments:

1. derive the `TuneState` candidate/smoothed state from the actual same-history
   WPE, adaptive bandpass, sigma statistic and commit cadence;
2. bind exp/trig and matrix factorization/eigensolver witnesses to deployed
   arithmetic and later rigorous finite-precision residuals;
3. bind applied R_acc/R_mag, band-noise floor, scheduler tolerance and all guards
   to the same frontend/runtime predecessor;
4. materialize asynchronous magnetometer continuation, all nonfinite/rejection
   branches, actual H18->A21 release and every literal runtime prefix;
5. attach corrected COMPLETE-BRMM and BIAS0/1/2 admission to the finite graph.

Only after that master passes `assert_finite_storage_master` may the required
high-precision feasibility diagnostic and first common joint24 storage search
run. No rho/metric search is authorized earlier.

## Gates

`ALT_LIVE_PASS=false`

`ALT_STARTUP_PASS=false`

`ALT_END_TO_END_PASS=false`

Original P4/P5 gates remain untouched. There is no certified rho, ultimate
bound, retained basin or capture time yet.
