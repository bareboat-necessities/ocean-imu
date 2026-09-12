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
require `zero_wind_heel_scope_enforced=true`; a later storage search cannot
silently reintroduce the excluded branch.

## Intended architecture and anti-dead-end guard

Construct the actual finite same-history runtime word first, then prove a
coercive dissipativity/storage inequality

`V_next <= rho V + w'Gamma w + c'Beta c`, `0<rho<1`.

Do not substitute random seeds, captured/replayed traces, frozen gains,
pointwise Jacobian products, independent coefficient boxes, covariance
consistency, statistical source events, wordwise resets or storage/rho searches
for a missing runtime/source relation. Unknown state errors may not be relabeled
as bounded disturbances. The complete finite master must pass
`proof_plan.assert_finite_storage_master` before any common-M search begins.

The finite graph must retain actual physical attitude increments and defects,
correlated v/p/S/a moments, one Live S origin, joint24 state, full 21 covariance,
one BIAS root, literal prediction/measurement/hygiene branches, frontend/WPE/
band/tuner memory, async magnetic state, scheduler credit, every H18/A21 edge,
deployment arithmetic residuals and explicit zero-heel scope ancestry.

## Current finite-runtime advancement

PR #523 has finite descriptors for physical prediction, accepted/rejected
measurements, full covariance, runtime OU/BA roots, attitude F/Q, integrated-OU
Qaxis, pending a_w synchronization, S scheduling/service, SafeLDLT branches,
accelerometer vibration guard/Racc, held-sample forcing, private Mahony,
WPE/band/stillness and staged tuner commits.

The startup magnetic path is structurally attached through persistent Mahony,
world-frame gravity admission, async wrapper clocks, default MagAutoTuner,
yaw-stripped tilt frame, physical magnetic source and the same
`PhysicalKinematics` ancestry as the finite physical word.

Under zero heel, the accepted boat quaternion directly supplies the covariance
yaw axis

`u_down_body = R(q)^T e_z = (2(xz-wy), 2(yz+wx), 1-2(x^2+y^2))`.

## MAG-BMM150-DET-v1 and startup capture

ALT admits commissioned installations satisfying:

- `20 <= ||B_W|| <= 75 uT`;
- horizontal field `>=15 uT`;
- body hard iron `<=5 uT`;
- deterministic residual `<=2 uT` per theorem sample.

With startup gravity-direction error <=0.02 rad, tilt contributes <=1.5 uT,
so total deterministic horizontal perturbation is <=8.5 uT and
`|sin(delta_yaw)| <= 17/30`. Since
`sin(0.61) >= 0.61-0.61^3/6 > 17/30`, yaw error is <0.61 rad. The SO(3)
triangle inequality gives total startup attitude error <0.63 rad <pi/4. No
`1/sqrt(N)` statistical reduction is used. Deployment atan2/AngleAxis/
normalization correspondence remains open.

## Fresh H18 and first Live sample now composed

The zero-heel gauged handoff is composed through shipping
`goLive -> initialize_from_attitude -> enterLive_`. It installs the attitude
covariance, seats `P_aw,aw` on the SAME committed `Sigma_aw`, clears every a_w
cross covariance, requires the same committed Live `R_S` and keeps BA learning
disabled.

`finite_fresh_joint24_entry.py` derives the actual H18 joint24 coordinates from
the SAME physical `Reference` and estimator state:

`c=Cayley(q_true_WB*conjugate(q_hat_WB))`,
`e_bg=b_g-b_g_hat`, `e_v=v-v_hat`, `e_p=p-p_hat`,
`e_S=S_centered-S_hat`, `e_aw=a-a_w_hat`,
`e_ba=beta-b_a_hat`, final coordinates `beta=beta_true`.

The one-time Live origin and fresh centered physical S=0 are enforced. No
fresh-entry covariance-consistency assumption or independent error box is used.

`finite_startup_live_runtime_bridge.py` preserves TunerReady Mahony/WPE/band/
stats/stillness/vibration-guard memory across goLive, derives active parameters
from the SAME carried TuneState, retargets the persistent S scheduler and
preserves any online pending-tune bit into the first Live boundary.

`finite_startup_first_live_step.py` now substitutes that exact bridge into the
existing Live IMU prefix. The first represented prediction, S-service decision,
held accelerometer event and tuner/WPE suffix therefore start from startup
ancestry rather than a synthetic Live root. This is still a conditional
real-arithmetic prefix, not the complete source-uniform word.

## H18/A21 hybrid language

`MAG-CALL-SCHEDULE-v1` requires first post-Live mag call <=40 ms and later gaps
<=40 ms. Shipping counts attempted post-delay `updateMag()` calls independent of
innovation acceptance, so its 250-count internal lock clears within 10 s and
the strict >1 s guard is automatically met.

Do not assume eventual A21 under arbitrary external hold. The graph retains:
no hold -> exact H18->A21 floor edge; held -> H18 may persist indefinitely;
first later release while Live -> exact H18->A21 edge; asserting hold in A21 ->
H18 with BA cross-covariances zeroed.

## Current blockers

The immediate blockers are now:

- deployment/binary32 correspondence for startup yaw extraction, atan2,
  AngleAxis, quaternion normalization, handoff setters and clocks;
- compose asynchronous magnetometer events and the firing Live tilt-reset edge
  into startup-rooted successive Live prefixes, not only as separate local maps;
- source-qualify the complete provisional/refinement/continuous-hard-iron
  magnetic schedule required by the declared scope;
- carry COMPLETE-BRMM/BIAS ancestry and every solver/hygiene/finite-precision
  branch through an arbitrary 600-step word, including both H18 and A21/hold
  continuations;
- prove the exact fresh/source-produced states land in a retained storage basin;
  only after the complete finite master passes its guard may common joint24
  storage/rho feasibility be attempted.

None may be replaced by trace replay, statistical concentration, or a
convenience entry set.

## Promotion state

`ALT_LIVE_PASS=false`

`ALT_STARTUP_PASS=false`

`ALT_END_TO_END_PASS=false`

Existing `P4_PASS` / `P5_MAY_START` remain independent and untouched.
