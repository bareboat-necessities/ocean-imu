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

This materially advances the finite-map stage but does NOT satisfy the universal
source quantifier. Exp/trig and numerical factorization witnesses, deployment
roundoff, remaining hybrid/runtime branches and source admission still require
same-history closure.

## Explicit magnetic source-specification blocker

The current canonical COMPLETE-BRMM/BIAS assumptions bound vessel motion and the
bias families, but the ALT source audit has not found a declared deterministic
envelope for all three quantities entering the startup magnetic identity:

`m_raw_B = R_true B_world + b_HI_body + n_mag_body`.

In particular, no current canonical theorem assumption has been identified that
gives a finite bound for the world magnetic-field magnitude, body-fixed hard
iron, and deterministic magnetometer residual. `R_mag`, bench statistics,
simulation values and sensor datasheet typical noise are not interchangeable
with such deterministic source assumptions.

Therefore `finite_mag_source_qualification` is intentionally fail-closed:
local finite identities may continue, but magnetic source qualification,
complete-word qualification and storage search must remain false until named
source/theorem assumptions provide those envelopes. Do not invent numerical
values merely to unblock the proof. This is an E-type source-specification
obligation unless an already-authoritative assumption is located and attached.

## Promotion state

`ALT_LIVE_PASS=false`

`ALT_STARTUP_PASS=false`

`ALT_END_TO_END_PASS=false`

Existing `P4_PASS` / `P5_MAY_START` are not controlled by ALT and remain
untouched by this plan.
