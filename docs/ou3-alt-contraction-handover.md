# OU-III parallel ALT proof handover

## Resume point and independent scope

Continue PR #523 on its current branch, preserving concurrent commits. Read
`AGENTS.md`, `docs/ou3-alt-proof-plan.md`, this handover, the ALT research ledger,
`docs/ou3-alt-finite-measurement-proof.md`,
`docs/ou3-alt-finite-core-composition.md`, `docs/ou3-alt-runtime-primitives.md`,
`docs/ou3-alt-mahony-binary32.md`, and `docs/ou3-brmm-main-handover.md`.
The original P2/P3/P4/P5 route remains independently continuable; its premises
and gates must not be weakened. `P3=1e-18` remains frozen.

## Immutable physical contracts

Retain joint24 `z=(c,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta)`, full 21-state covariance,
all motion/bias cross terms, one persistent Live S origin and the corrected
COMPLETE-BRMM envelope including acceleration <=8.8 m/s^2 and all-time centered
primitive `D_S<=1100 m*s`. Preserve same-signal WPE/bandpass/sigma/tau/T_S
ancestry, staged tuner state and separate BIAS0/1/2 histories. Unknown state
errors may not be renamed as bounded disturbances.

## Available finite graph components

The local identities retain actual continuous physical attitude increments and
angular/model defect, correlated physical translation moments, gyro-bias drift,
one shared accelerometer-bias driver, finite accelerometer/magnetometer/S
residual secants, inverse-free `K Sigma=N`, quaternion injection, same-beta
projection and full Joseph/reset covariance. H18 retains latent BA covariance;
there is no 18-state A21 storage substitute.

`finite_core.py` carries finite mean/error, nominal attitude, full covariance and
physical predecessor. Prediction modules generate AA/linear/BA covariance
blocks from mean/runtime coefficients and retain OU/BA roots, structured/fast
attitude Q, IntegratedOUChain Qaxis, masks and numerical-hygiene branches.
`finite_post_prediction.py` retains pending a_w sync, scheduler credit/due
recurrence and S service. `finite_measurement_runtime.py` retains first LDLT
success, deterministic bump/retry and double-failure rejection. Numerical
factorization witnesses are conditional, not freely admitted source variables.

Source-side modules represent raw packet/de-heel/temperature-model relations,
private vertical observation, WPE, adaptive band/noise gain, debiased variance,
stillness and tuner candidate/commit boundaries. Temporal composers and active-
parameter interfaces retain these outputs across represented events. A type,
provenance token, local equality or completed trace does not qualify the entire
physical source history or establish complete runtime-prefix coverage.

The staged commit must keep one candidate TuneState as the ancestor of applied
OU tau, stationary Sigma_aw, pseudo-update cadence and anisotropic R_S. Online
`apply_ou_tune_(false)` does not queue posterior a_w synchronization; discrete
sync remains separate. H18 holds the estimated bias, so its error satisfies
`e_ba+ = e_ba + (phi_true-1) beta + u_b`, not merely `e_ba += u_b`.

## Initialized private-Mahony numerical binding

`finite_binary32_mahony.py` supplies an exact finite initialized observer graph
under `binary32-rne-gradual-no-fma-eigen3-scalar`. It computes, rather than
accepts, both normalization reciprocals from the same rounded norm sums using
the literal integer seed and Newton program. An independent midpoint-cell
checker verifies basic-operation rounding; each exact rounding defect retains
its own operands. Zero/subnormal norm sums are represented.

The existing vertical API accepts `arithmetic_profile` for this path, returns
the same `V.Result`, and rejects supplied seed/reciprocal witnesses. The raw
packet bridge forwards that option; subsequent WPE/band consumers can share the
one output. The default real helper remains conditional. Initial seeding,
nonfinite/overflow branches and actual target compiler/FMA/reduction-profile
qualification remain OPEN. Do not replace the norm-defective float quaternion
by a unit quaternion, or promote host correspondence to deployment proof.

Run the new exact/native checks with:

```sh
PYTHONPATH="$PWD:$PWD/tools/stability:$PWD/tests/ou3_alt_contraction" \
  python3 -m unittest test_finite_binary32_mahony \
  test_finite_vertical_complementary_runtime test_core test_bias_families -v
```

Native correspondence uses g++ and Eigen, observes real startup/update calls,
and never forces a runtime root. Focused CI requires the native dependency;
local runs can specify `EIGEN_INCLUDE_DIR`. Tests are regressions, not physical
source admission or contraction evidence.

## Decisive remaining object

**The complete source-uniform finite 600-step runtime word is still open.**
The limiter remains C/E source/runtime attachment, not demonstrated instability
or failed common storage. Complete the remaining same-predecessor operands and
branches behind the available components, including applied R_acc/R_mag,
transcendental and solver arithmetic, asynchronous magnetic continuation,
all H18/A21 transitions and every literal prefix. Attach the corrected physical
and bias admission to that graph. A finite branch map alone proves neither that
startup reaches it nor that subsequent execution remains there.

Only after the complete master passes `assert_finite_storage_master` may the
high-precision feasibility diagnostic and common joint24 storage search run.
Then close useful supply/ultimate bound, every-prefix retention, hybrid landing,
fresh entry, both startup paths and deployment finite precision. No rho search
is authorized by the numerical binding's successful tests.

## Gate state

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.
Original P4/P5 are not promoted by ALT. No certified rho, retained basin,
ultimate bound or finite capture time is available yet.
