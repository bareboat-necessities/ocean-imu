# OU-III parallel ALT proof handover

## Resume point and independent scope

Read `AGENTS.md`, `docs/ou3-alt-proof-plan.md`, this handover,
`docs/ou3-alt-finite-measurement-proof.md`, `docs/ou3-alt-finite-core-composition.md`,
`docs/ou3-alt-runtime-primitives.md`, `docs/ou3-alt-contraction.md`, the ALT section
of `docs/ou3-proof-research-state.md`, and `docs/ou3-brmm-main-handover.md`.
Continue PR #523 on its branch. The original P2/P3/P4/P5 route remains
independently continuable and unchanged; `P3=1e-18` remains frozen.

## Target and immutable physical contracts

Retain `z=(c,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta)` in R^24, all motion/bias cross
terms, one persistent Live S origin and coercive storage eventually satisfying
`V_next <= rho V + supply`, `0<rho<1`. Unknown motion errors may not be renamed
as bounded inputs. Corrected COMPLETE-BRMM retains acceleration <=8.8 m/s^2,
all-time centered primitive `D_S<=1100 m*s`, zero lever arm, same-signal
WPE/bandpass/sigma/tau/T_S ancestry, staged tuner state and separate BIAS0/1/2
histories.

## Finite graph now proved/represented

The supplying finite identities retain the actual continuous physical quaternion
increment and angular/model defect, correlated physical translation moments,
gyro-bias drift, one shared physical accelerometer-bias driver, finite
accelerometer/magnetometer/S residual secants, inverse-free `K Sigma=N`, both
quaternion injection branches, same-beta radial projection and full 21-state
Joseph/reset covariance. H18 retains latent BA uncertainty; no 18-state marginal
replacement is used.

`finite_core.py` recursively carries finite joint24 mean/error, nominal attitude,
full covariance and physical predecessor. Every represented accepted event uses
the immediately preceding pair; there is no Jacobian cocycle or fresh covariance.

The prediction path is now composed through:

- `finite_prediction_covariance.py`: full AA/linear/BA block congruence retaining
  every cross covariance;
- `finite_ou_runtime_primitives.py`: one `(h,tau,alpha)` root supplies the mean
  OU coefficients and `F_LL`; one `(tau_b,phi_b,Q_bacc)` root supplies the BA
  mean factor, BA cross factors and `Q_BB`; H18 is exactly held;
- `finite_attitude_runtime.py`: constant-rate `R/B`, structured/fast attitude-Q,
  isotropic versus Simpson branches, `integral_B`, and 6x6 PSD-hygiene control;
- `finite_qaxis_runtime.py`: literal small/general `IntegratedOUChain<3>` Qaxis
  formulas, nested 3x3 marginal hygiene and final 4x4 hygiene;
- `finite_prediction_runtime.py`: highest prediction entry. It accepts no
  precomputed shipping `F_AA/Q_AA/F_LL/Q_LL/Q_BB` or phi factor and enforces the
  SAME bias-corrected gyro for nominal quaternion and attitude covariance.

Post-prediction/runtime control is represented by:

- `finite_post_prediction.py`: pending a_w covariance synchronization as the
  positive spectral part of `target-P_aw`, including no-pending, eigensolver
  failure and success branches; symmetry hygiene; scheduler elapsed-credit and
  due/not-due recurrence; due S service routed through the measurement runtime;
- `finite_measurement_runtime.py`: literal shipping `safe_ldlt3_` first success,
  one deterministic diagonal-bump retry, and double-failure rejection. Accepted
  retry uses the SAME shifted innovation in gain and Joseph; rejection preserves
  state/covariance.

These are exact real-arithmetic/control-flow descriptors. Trig/exp values,
Eigen branch outcomes and machine epsilon are explicit witnesses, not silently
asserted source facts or finite-precision certificates.

## Implementation correspondence, not source qualification

`shipping_finite_identity.py` still passively compares instrumented and unchanged
shipping wrapper executions from startup and checks selected H18, A21 and
H18->A21 600-step windows. It forces no estimator root. This remains a regression
for implementation correspondence only, never a universal source proof.

The focused CI job runs all finite modules plus correspondence and keeps
`ALT_LIVE_PASS`, `ALT_STARTUP_PASS`, `ALT_END_TO_END_PASS`, source-uniform word
qualification, storage feasibility and deployment-roundoff closure fail-closed.
The inherited broad ALT suite is separate and may still fail the unchanged
Mahony prerequisite; do not weaken it.

## Current controlling gap

**The complete source-uniform finite 600-step runtime word is still NOT
materialized.** The remaining limiter is C/E (same-history runtime/source
attachment), not demonstrated filter instability or failed common storage.

The largest remaining attachments are:

1. bind the scalar exp/trig and PSD/LDLT/eigensolver witnesses above to actual
   deployed arithmetic and, later, rigorous finite-precision residuals;
2. bind `tau`, `Sigma_aw`, `tau_b`, `Q_bacc`, a_w floor target, applied
   `R_acc/R_mag/R_S`, and scheduler period/tolerance to the SAME frontend/tuner
   predecessor rather than independent operands;
3. materialize WPE, adaptive bandpass, sigma statistic, candidate/active tuner
   state and staged commits, preserving same-signal ancestry;
4. materialize asynchronous magnetometer state, all measurement reject/nonfinite
   branches, H18->A21 release/guard and every literal runtime prefix;
5. attach COMPLETE-BRMM and BIAS0/1/2 admission to that graph.

Only after the finite master closes may `assert_finite_storage_master` admit it
to the required high-precision complete-word feasibility diagnostic and then a
common joint24 dissipativity-storage search. No rho/metric search is authorized
before that point.

## Gates

`ALT_LIVE_PASS=false`

`ALT_STARTUP_PASS=false`

`ALT_END_TO_END_PASS=false`

Original `P4_PASS` / `P5_MAY_START` remain untouched. There is no certified rho,
ultimate bound, retained basin or capture time yet.
