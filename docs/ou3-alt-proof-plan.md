# OU-III ALT proof plan

## Independent scope

This is a parallel proof architecture for the current shipping OU-III filter.
It does not replace, weaken, delete or become a prerequisite of the existing
P2/P3/P4/P5 proof track. Keep that route independently continuable.

The target remains an end-to-end theorem for every admitted corrected
COMPLETE-BRMM history, every admitted BIAS0/BIAS1/BIAS2 history and every
declared disturbance/finite-precision history inside the declared ALT deployment
scope. ALT preserves the actual Mahony/proxy/startup path, H18/A21 hybrid logic,
frontend/tuner memory, full 21-state covariance and joint24
motion/error/true-bias state.

## Deployment-scope exclusion: no wind heel

The optional shipping wind-heel retarget feature is excluded from ALT.
The certified deployment language requires `wind_heel_rad_ == 0` from
construction onward and zero calls to `update_wind_heel()`.

This is a theorem-scope restriction, not a shipping change: shipping initializes
`wind_heel_rad_` to zero. On this scope B'=B, `deheel_vector_` is the identity,
and there is no wind-heel/body-frame retarget event in the hybrid language.
Executions with nonzero wind heel or any dynamic wind-heel update are outside
ALT until a future proof explicitly widens the theorem.

`deployment_scope.py` encodes this fail-closed assumption, and the storage guards
now require `zero_wind_heel_scope_enforced=true`; a later storage search cannot
silently reintroduce the excluded branch.

## Intended architecture

Construct the actual finite same-history runtime word first, then prove a
coercive dissipativity/storage inequality of the form

`V_next <= rho V + w'Gamma w + c'Beta c`, `0<rho<1`,

with bounded neutral/source supply. Use exact finite descriptors and inverse-free
innovation relations where possible. Unknown state errors may not be relabeled
as disturbances. Physical/source correlations must survive composition.

Only after the finite runtime relation is complete may a high-precision
feasibility diagnostic or common joint24 storage search start.

## Mandatory no-dead-end guards

The following are not theorem work and may not be substituted for a missing
runtime/source relation:

- more random seeds, captured traces or replay-fitted roots;
- frozen gain/covariance observers;
- pointwise Jacobian products advertised as a finite physical map;
- independently boxed coefficients that destroy same-history ancestry;
- wordwise S reset, position reanchor or a convenience entry set;
- Gaussian/high-probability replacements for deterministic source admission;
- measurement covariance or sensor statistics substituted for a deterministic
  theorem disturbance envelope;
- storage/rho searches before the complete finite master passes the
  representation guard;
- changing the shipping filter solely to make the proof easier;
- reintroducing the excluded wind-heel branch without an explicit theorem-scope
  widening and proof of its retarget event.

If a proposed experiment cannot falsify or close a stated theorem obligation,
do not spend proof effort on it. The two-strike architecture-review rule applies
to genuine completed proof attempts, not incomplete source graphs or test bugs.

## Finite representation guard

Before any common-M, rho, endpoint refinement or high-precision storage attempt,
`proof_plan.assert_finite_storage_master` must validate the representation. A
Jacobian cocycle, source token, successful execution or collection of metadata
flags cannot replace an exact finite runtime graph or a proved complete anchored
mean-value relation.

The finite graph must retain at least:

- actual continuous physical attitude increments and angular/model defects;
- correlated physical v/p/S/a moments and one persistent Live S origin;
- joint24 `(c,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta)` with all cross information;
- full 21-state covariance in H18 and A21;
- one same-history BIAS driver/root per word;
- literal prediction coefficients/process covariance and all hygiene branches;
- applied measurement R, safe-LDLT accept/retry/reject, Joseph/reset/projection;
- frontend/WPE/bandpass/sigma/tuner memory and staged commits;
- asynchronous magnetic state, scheduler credit and every H18/A21 edge;
- every completed literal runtime prefix;
- deployment arithmetic residuals before final theorem promotion;
- explicit zero-wind-heel deployment-scope ancestry.

## Current finite-runtime advancement

PR #523 has exact/local descriptors for physical prediction, accepted finite
measurements, full covariance composition, runtime OU/BA roots, attitude F/Q,
integrated-OU Qaxis, pending a_w covariance synchronization, scheduler/S service
and safe-LDLT branches. The highest prediction entry accepts no precomputed
shipping transition/process matrices and enforces one bias-corrected gyro across
nominal and covariance attitude propagation.

The startup magnetic path is structurally attached through the persistent
Mahony proxy, world-frame gravity-alignment gate, literal asynchronous admission
clock, default MagAutoTuner accept/reject/ready recurrence, yaw-stripped tilt
frame, raw true-field/hard-iron/residual source identity, and the same main
`PhysicalKinematics` endpoint used by the finite physical word.

The startup north-ready path retains the generation-zero `setMagWorldRef_`
write, pending absolute-yaw gauge, proxy-to-MEKF handoff seed, shipping 0.035 rad
tilt covariance, 0.087 rad gauged yaw covariance versus 1.5708 rad free-yaw
covariance, and the real-arithmetic core of `initialize_from_attitude`. That
reset zeroes the attitude error state, replaces the attitude 3x3 covariance by
the tilt/yaw projector split, zeroes only attitude<->gyro-bias covariance, and
preserves the remaining state/covariance entries.

Under the zero-wind-heel scope, the remaining B' correspondence disappears:
for a unit accepted boat quaternion q=(w,x,y,z), the covariance yaw axis is now
attached directly as

`u_down_body = R(q)^T e_z = (2(xz-wy), 2(yz+wx), 1-2(x^2+y^2))`.

No nontrivial body-prime rotation is admitted in ALT.

## MAG-BMM150-DET-v1 deterministic magnetic admission

ALT declares the following engineering source class for a commissioned BMM150
installation. These are theorem admission limits, not Bosch guarantees for
arbitrary mounting environments:

- `20 uT <= ||B_W||_2 <= 75 uT`;
- `||(B_W.x, B_W.y)||_2 >= 15 uT`;
- `||b_HI_body||_2 <= 5 uT`;
- `||n_mag_body||_2 <= 2 uT` for every accepted theorem sample.

The horizontal lower bound is essential for deterministic north/yaw
observability. The hard-iron and residual limits are commissioned source
requirements, not replacements for `R_mag`.

## Deterministic startup north and real-arithmetic attitude capture closed

The canonical operating domain declares startup world-averaged gravity direction
error `<= 0.02 rad`. Using `sin(x) <= x`,

`2 * 75 uT * sin(0.02/2) <= 1.5 uT`.

Therefore the complete deterministic perturbation of the startup magnetic mean
is

`E <= 5 + 2 + 1.5 = 8.5 uT`.

No `1/sqrt(N)` factor is used. Since the true horizontal field is at least
15 uT, the learned horizontal mean cannot vanish and

`|sin(delta_yaw)| <= 8.5/15 = 17/30`.

For `x=0.61`, the alternating Taylor lower bound gives

`sin(x) >= x - x^3/6 > 17/30`,

so `|delta_yaw| < 0.61 rad`. SO(3) geodesic triangle inequality with the 0.02 rad
tilt bound gives total attitude error `< 0.63 rad < pi/4`. Thus the
source-qualified real-arithmetic startup attitude enters the declared 45-degree
radius. Deployment `atan2`/AngleAxis/quaternion normalization and binary32
correspondence remain open and prevent promotion.

## Fresh real-arithmetic H18 Live entry now composed

The zero-heel gauged handoff is now composed through shipping `goLive` and
`enterLive_` at the exact real-arithmetic state/covariance level.

The composition retains shipping order:

1. `initialize_from_attitude` installs the accepted handoff quaternion and
   attitude covariance;
2. `enterLive_` commits the qualified tuner operating point;
3. `reset_aw_covariance_to_stationary()` replaces `P_aw,aw` by the SAME committed
   `Sigma_aw` and clears every a_w cross-covariance;
4. because wrapper `goLive` passes `allow_acc_bias=false` and the startup bias
   lock is still engaged, accelerometer-bias learning remains disabled;
5. Live `R_S` is required from the same committed parameter state;
6. startup stage becomes Live with stage timer zero.

Therefore the fresh scoped real-arithmetic entry is explicitly an H18 entry,
not an assumed covariance-consistency set. This still does not prove the first
binary32/deployed Live prefix or the eventual H18->A21 release.

## Current blockers

The immediate blockers are now:

- close deployment/binary32 correspondence for startup yaw extraction, `atan2`,
  AngleAxis, quaternion normalization, norm gates and the handoff reset;
- attach the wrapper `stage_=Live`, `live_time_sec_=t_` and inner-stage transition
  to one same-history finite event including binary32 clock arithmetic;
- source-qualify the complete asynchronous magnetometer call schedule through
  tuner-ready, provisional lock, Live refinement and subsequent updates;
- connect the fresh H18 state above to the already materialized first literal
  Live IMU/pseudo/mag prefix without a gap in frontend/tuner/scheduler ancestry;
- close the eventual H18->A21 release edge with the required accepted magnetic
  count/time guard and same-history source continuation;
- complete deployment arithmetic residuals and remaining solver/hygiene branches
  before any storage search.

None of these may be replaced by trace replay, statistical concentration, or an
assumed fresh-entry covariance-consistency condition.

## Promotion state

`ALT_LIVE_PASS=false`

`ALT_STARTUP_PASS=false`

`ALT_END_TO_END_PASS=false`

Existing `P4_PASS` / `P5_MAY_START` are not controlled by ALT and remain
untouched by this plan.
