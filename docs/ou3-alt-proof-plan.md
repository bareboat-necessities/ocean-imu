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

## Pre-rho attitude representation

`ou3-alt-attitude-atlas.md` supplies the four-chart cover and exact transport
used by the current finite runtime. It covers ungauged timeout entry, including
true yaw pi, while preserving joint24 and the full 21-state shipping covariance.
`ou3-alt-startup-pre-rho.md` retains the exact single-Cayley obstruction as a
regression. Unknown attitude terms remain same-state rational graph functions;
they may not be relabeled as independent bounded disturbances.

The interleaver retains ungauged acquisition across Live and uses MEKF tilt for
initial north acquisition after Live. Continuous calibration/refinement retain
their own private-observer frame. Later north starts the magnetic service clock
without resetting the physical Live/S origin or duplicating packet statistics.

The guard still requires universal source/control reachability and source-uniform
arithmetic. Finite-word magnetic counter safety is closed by saturation for
all finite call counts; see `ou3-alt-deployment-prerequisite.md`. Measurement,
statistics and release continue at the cap, including after threshold changes.
The default RN32 timeout crossing is
sample 30,002; shared finite budgets include 600 more updates, conditional on
actual source-produced alignment. Coercive storage must respect the fact that
coordinate zero in charts 1..3 is a 180-degree error, not zero error.

## Current finite-runtime advancement

The graph has finite descriptors for physical prediction, accepted/rejected
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

The conditional magnetic calculation needs a full accumulation-to-handoff
frame bound, including physical heading variation. Gravity-direction error
<=0.02 rad alone does not supply it. With an independently established
half-frame sine bound <=0.01 and handoff tilt <=0.02, the 8.5 uT perturbation,
0.61-rad yaw and 0.63-rad full-angle bounds follow. Source qualification of those
premises and deployment atan2/AngleAxis/normalization remain open.

The finite atlas word does not require this small-angle accuracy certificate.
It retains the full frame, nonzero quaternion and actual magnetic discrepancy.
The universal proper-rotation chord bound 2 gives a finite mean-perturbation
image bound of 157 uT; no north-nonvanishing or contraction claim follows.
Small magnetic-frame accuracy is therefore a later basin/usefulness question,
not an independent prerequisite for the current finite-word representation.

## Fresh H18 and first Live sample now composed

The zero-heel gauged handoff is composed through shipping
`goLive -> initialize_from_attitude -> enterLive_`. It installs the attitude
covariance, seats `P_aw,aw` on the SAME committed `Sigma_aw`, clears every a_w
cross covariance, requires the same committed Live `R_S` and keeps BA learning
disabled.

`finite_fresh_joint24_entry.py` derives the actual H18 joint24 coordinates from
the SAME physical `Reference` and estimator state:

`(chart,c)=Atlas(q_true_WB*conjugate(q_hat_WB))`,
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

`finite_live_interleave.py` now extends that startup bridge to successive IMU,
magnetic and external-hold events. The magnetic graph includes the default
continuous estimator carried from startup, literal refinement/reset clocks,
reference/yaw writes, coupled hard-iron/reference application and same-packet
measurement/count logic. `ou3-alt-live-magnetic-word.md` supplies the finite
composition theorem and uniform finite-real calibration bounds. These remove
detached calibration/reference operands; they do not qualify the complete word
or permit a storage search.

The interleaved IMU edge now uses `step_from_shipping_operands` in
`finite_live_tilt_prefix.py`. `finite_tilt_reset_runtime.py` derives the watchdog
predicate from the actual post-accelerometer nominal quaternion and derives a
firing `initialize_from_acc_preserve_yaw` quaternion/down axis from that same
pre-reset attitude and the exact guarded accelerometer. The ideal-real
preserve-yaw graph is algebraic through normalized gravity and half-angle
identities; the final reset quaternion is no longer a free theorem operand.
The >70-degree test uses a rigorous rational enclosure of the shipping
`acos(cos_tilt)*57.295779513f` threshold. The tiny enclosure boundary and actual
binary32/libm branch correspondence remain explicitly fail-closed.

## Joined machine startup attachment

`finite_startup_joined_machine_history.py` roots the guard/private-observer/
LPF/stillness/band/TuneState product at reset. Ordinary first-sample seeding is
computed from the same accelerometer by `finite_binary32_mahony_startup.py`;
near-antiparallel Eigen SVD remains unqualified. Startup samples bind the same
machine source to the consumed band and sigma events. Pending boundary and
conditional goLive preserve upstream memory, and `admitted_live` substitutes
that same result through the strong scheduler/prediction/measurement/Racc chain
without new Live snapshots or a second goLive execution.

Racc consumes pre-candidate TuneState sigma every sample, independently of the
pending active-parameter commit, and uses the raw preupdate WPE getter/prior.
The joined Racc/accelerometer relation consumes the same machine-conditioned
operand; its source-uniform arithmetic supply remains open.
Conditional attachment does not close universal startup capture, ungauged
handoff, deployment arithmetic or the source-uniform 600-transition word.
All finite-master/storage/theorem gates remain fail-closed.

Startup raw-IMU observability must be qualified separately from the post-Live
arbitrary-bounded ISS quantifier. The latter permits no inference about the
first-sample norm. `ou3-alt-startup-disturbance-contract.md` proves that an
arbitrary bounded residual can keep that norm below the initialization
threshold while preserving the correct gravity direction. A longer timeout,
more Mahony bounds or magnetic service cannot initialize that proxy. Retain
the obstruction and obtain the intended startup sensor contract; do not
silently choose a residual cap or reinterpret covariance as one.

## Physical source-prefix graph

`finite_source_bound_live_word.py` carries physical source and sensor ancestry
with the runtime state. Source qualification cannot follow from matching tokens.
The finite constructor now checks existing physical vector caps, the joint
nine-component acceleration-moment IQC, necessary rotation/rate constraints and
BIAS component/norm envelopes. Consecutive segments retain one actual phi.
`finite_brmm_moment_prefix.py` derives cumulative moments and the prefix budget
from the same segments using the exact Gramian concatenation identity; the
proof is in `ou3-alt-source-continuation.md`.

These executable checks are necessary outer constraints, not a runtime
membership oracle. `finite_complete_brmm_restriction.py` closes the universal
physical-source implication separately: the primary theorem source is the same
bounded physical history with bounded centered primitive, and restricting any
admitted history to the 5 ms grid yields the finite caps, primitive recurrence,
coupled moment IQC and rotation bound. A common spectral/shaping generator is
explicitly not required because those constructions are only sufficient
certificate methods. Matching labels still cannot admit an arbitrary runtime
trace. The frequency/frontend relation, actual BIAS generating functions and
deployment arithmetic remain open. None of this authorizes rho search, shrinks
the filter-error domain, or replaces native arithmetic closure.

The source-owning Live word also has a dedicated fresh-origin endpoint. It checks
the actual sample-zero `Reference` at `time == live_origin`, requires centered
S=0, retains the same history/BIAS root and physical vector/bias caps, and permits
a magnetic call before transition 1 without inventing a 5 ms segment. That call
does not advance source ordinal or time. This closes only the checked outer
endpoint topology; complete source membership at sample zero remains false.

## H18/A21 hybrid language

`MAG-CALL-SCHEDULE-v1` requires first post-Live mag call <=40 ms and later gaps
<=40 ms. Shipping counts attempted post-delay `updateMag()` calls independent of
innovation acceptance, so its internal lock clears within 10 s. The 250th call can occur before the
strict >1 s guard: the corrected proof uses
`first_gap + max(249*gap, 1+gap)` and continued, locally finite call coverage.
Finite-prefix checks do not prove the infinite schedule. The source-bound
saturation invariant preserves every configurable count threshold for all finite
prefixes. Other arithmetic operations and floating clocks require their own
qualification before total shipping execution can be claimed.

Do not assume eventual A21 under arbitrary external hold. The graph retains:
no hold -> exact H18->A21 floor edge; held -> H18 may persist indefinitely;
first later release while Live -> exact H18->A21 edge; asserting hold in A21 ->
H18 with BA cross-covariances zeroed.

## Current blockers

The attitude representation, conditional ungauged continuation, certified empty
magnetic startup prefix and signed magnetic-counter safety are attached.
The unresolved qualifications are:

- deployment/binary32 correspondence for startup yaw extraction, atan2,
  AngleAxis, quaternion normalization, handoff setters and clocks;
- close deployment correspondence for the now-bound Live tilt reset, including
  the threshold-boundary sliver, normalization cutoff, sqrt/atan2/asin/
  AngleAxis/libm execution and nonfinite branches; do not reintroduce a free
  final reset quaternion;
- source-qualify the complete provisional/refinement/continuous-hard-iron
  history and its remaining arithmetic/solver branches; bounded finite-real
  calibration is not proof of calibration accuracy or a useful ISS margin;
- carry COMPLETE-BRMM/BIAS ancestry and every solver/hygiene/finite-precision
  branch through an arbitrary 600-step word, including both H18 and A21/hold
  continuations; retain the proved saturated-counter projection and qualify
  floating-clock behavior without substituting unbounded rational clocks;
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
