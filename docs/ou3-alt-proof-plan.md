# OU-III ALT proof plan

## Independent scope

This is a parallel proof architecture for the current shipping OU-III filter.
It does not replace, weaken, delete or become a prerequisite of the existing
P2/P3/P4/P5 proof track. Keep that route independently continuable.

The target remains an end-to-end theorem for every admitted corrected
COMPLETE-BRMM history, every admitted BIAS0/BIAS1/BIAS2 history and every
declared disturbance/finite-precision history. ALT must preserve the actual
Mahony/proxy/startup path, H18/A21 hybrid logic, frontend/tuner memory, full
21-state covariance and joint24 motion/error/true-bias state.

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
- changing the shipping filter solely to make the proof easier.

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
- deployment arithmetic residuals before final theorem promotion.

## Current finite-runtime advancement

PR #523 now has exact/local descriptors for physical prediction, accepted finite
measurements, full covariance composition, runtime OU/BA roots, attitude F/Q,
integrated-OU Qaxis, pending a_w covariance synchronization, scheduler/S service
and safe-LDLT branches. The highest prediction entry accepts no precomputed
shipping transition/process matrices and enforces one bias-corrected gyro across
nominal and covariance attitude propagation.

The startup magnetic path is now structurally attached through the persistent
Mahony proxy, world-frame gravity-alignment gate, literal asynchronous admission
clock, default MagAutoTuner accept/reject/ready recurrence, yaw-stripped tilt
frame, raw true-field/hard-iron/residual source identity, and the same main
`PhysicalKinematics` endpoint used by the finite physical word. The learned
magnetic reference cannot be freely reselected per event.

The startup north-ready path also retains the literal generation-zero
`setMagWorldRef_` write, pending absolute-yaw gauge, proxy-to-MEKF handoff seed,
shipping 0.035 rad tilt covariance, 0.087 rad gauged yaw covariance versus
1.5708 rad free-yaw covariance, and the real-arithmetic core of
`initialize_from_attitude`. The latter zeroes the attitude error state, replaces
the attitude 3x3 covariance by the tilt/yaw projector split, zeroes only the
attitude<->gyro-bias covariance, and preserves the remaining state/covariance
entries.

This materially advances the finite-map stage but does NOT satisfy the universal
source quantifier. Exp/trig and numerical factorization witnesses, deployment
roundoff, remaining hybrid/runtime branches and source admission still require
same-history closure.

## MAG-BMM150-DET-v1 deterministic magnetic admission

ALT declares the following engineering source class for a commissioned BMM150
installation. These are theorem admission limits, not Bosch guarantees for
arbitrary mounting environments:

- `20 uT <= ||B_W||_2 <= 75 uT`;
- `||(B_W.x, B_W.y)||_2 >= 15 uT`;
- `||b_HI_body||_2 <= 5 uT`;
- `||n_mag_body||_2 <= 2 uT` for every accepted theorem sample.

The total-field range encloses ordinary terrestrial geomagnetic magnitudes with
margin while remaining far inside the BMM150 electrical range. The horizontal
lower bound is separate and essential: without a nonzero horizontal field no
uniform deterministic north/yaw capture theorem is possible. The 5 uT hard-iron
limit is a commissioned-placement requirement. The 2 uT deterministic residual
limit remains materially wider than the BMM150 RMS output-noise values of the
normal presets, but unlike RMS noise it is a hard theorem admission limit. A
history that violates either condition is outside the theorem rather than being
silently absorbed into `R_mag`.

`finite_mag_source_qualification` checks the named assumption with exact norm
witnesses, and the theorem-facing startup magnetic edge refuses an unqualified
sample before it can mutate the tuner or magnetic wrapper clock. This closes the
prior missing-magnetic-envelope E blocker, but it does NOT yet qualify the whole
asynchronous magnetic schedule.

## Deterministic startup north capture now obtained

The canonical operating domain already declares startup world-averaged gravity
direction error `<= 0.02 rad`. Using `sin(x) <= x`, the tilt-frame chord term is
therefore bounded without a floating transcendental witness by

`2 * 75 uT * sin(0.02/2) <= 1.5 uT`.

Under MAG-BMM150-DET-v1 the complete deterministic perturbation of the startup
magnetic mean is consequently

`E <= 5 + 2 + 1.5 = 8.5 uT`.

No `1/sqrt(N)` factor is used: deterministic hard iron/residual may remain
coherent across all tuner samples. Since the true horizontal field is at least
15 uT, the learned horizontal mean cannot vanish and the exact source-uniform
direction relation gives

`|sin(delta_yaw)| <= 8.5/15 = 17/30`.

This is now a proved finite capture supply relation. It is not yet the complete
startup theorem: the corresponding `atan2`/binary32 angle enclosure and its
composition with the 0.02 rad tilt error must still be certified against the
declared 45-degree full-SO(3) fresh-entry radius.

## Current startup blockers

The immediate startup blockers are now:

- attach the `initialize_from_attitude` world-down axis to the accepted boat
  quaternion through the actual wind-heel/body-prime conversion and Eigen
  normalization;
- close binary32 correspondence for startup yaw extraction, `atan2`, AngleAxis,
  quaternion normalization, norm gates and the handoff reset;
- source-qualify the asynchronous magnetometer call schedule through tuner ready
  and the subsequent refinement/Live path;
- turn `|sin(delta_yaw)| <= 17/30` plus the 0.02 rad tilt bound into a rigorous
  binary32/full-SO(3) `<45 deg` handoff-entry certificate;
- compose `goLive`/`enterLive_`, fresh H18 entry and the existing Live word on the
  same physical/source history.

None of these may be replaced by trace replay, statistical concentration, or an
assumed fresh-entry covariance-consistency condition.

## Promotion state

`ALT_LIVE_PASS=false`

`ALT_STARTUP_PASS=false`

`ALT_END_TO_END_PASS=false`

Existing `P4_PASS` / `P5_MAY_START` are not controlled by ALT and remain
untouched by this plan.
