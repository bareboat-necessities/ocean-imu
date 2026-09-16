# Persistent machine CORE continuation

The complete source-uniform 600-step word requires each compiler event to
consume the preceding compiler event's full successor.  This includes joint24,
all 21 covariance rows in H18 and A21, nominal attitude and its atlas chart,
the same physical/bias reference, and the control state used by subsequent
branches.  Closing scalar arithmetic bounds alone does not prove this link.

The strongest WPE/clock interleaver now carries one complete CORE and event
control runtime for each compiler mode. It supplies those predecessors to
`finite_machine_active_prediction_roots.build` and
`finite_machine_prediction_displacement.compare`. The exact predictor continues
to consume the exact-shadow predecessor; each compiler predictor consumes its
own previous successor. Angular roots read the carried gyro-bias error, and the
accel-bias branch reads the carried H18/A21 mode. No attitude, covariance or
bias-error equality between the compiler and exact predecessors is required.
The same `finite_live_input_contract.PacketAdmission` binds RN32 of the source
gyro and accelerometer to both modes. Prediction retains the exact additive
gyro projection defect, `D * (RN32(raw_gyro) - raw_gyro)`, in the physical
rate/bias/residual identity. The finite machine predictor consumes that rounded
gyro, and the ungauged gravity gate consumes the same rounded gyro/accelerometer.
The exact shadow continues to consume its original source packet.

The regression retains the defect in the older local-only entry. A
pending stationary-variance difference produces 48 nonzero full covariance
discrepancies, with maximum approximately 0.0146197582243.  The following local
prediction starts from the shadow covariance and drops that predecessor
discrepancy. The repaired distinct-predecessor entry instead agrees with direct
prediction from the preceding compiler covariance. A separate nonzero
gyro-bias-error regression proves that its angular roots and attitude follow
the carried precursor too. These are conditional component fixtures, not
admitted-source counterexamples or native arithmetic measurements.

After the per-mode covariance floor, S service, applied R_S, Racc and
accelerometer correction, `finite_machine_core_continuation.advance_imu`
executes each mode's own watchdog predicate and timer, preserve-yaw reset when
it fires, and ungauged gravity gate. Reset witnesses bind to that mode's
current nominal attitude and the same executed machine-conditioned
accelerometer. The runtime retains its magnetic memory, control, call schedule
and clock. `advance_mag` executes the existing dual-clock magnetic algorithm
from that runtime and the same physical packet; `advance_hold` executes the
actual H18/A21 covariance transitions. Neither copies the shadow successor or
requires the modes to take the same branch.

The source-rooted startup constructor also executes the literal per-compiler
`reset_aw_covariance_to_stationary` effect from the already-produced goLive
commits: all nine a_w covariance entries receive that mode's Sigma_aw, and the
108 cross entries in its a_w rows/columns become zero. The startup goLive
producer remains attached to the carried CORE history. The actual startup
machine private observer is retained for sample-zero magnetic calls as well as
later IMU events.

`finite_machine_startup_core.py` now supplies a sealed default construction
root, retained by the actual startup `State` through every IMU and pending
boundary. The selected application leaves Pq0=5e-4, Pb0=1e-6 and the initial
linear/accel-bias settings at their header defaults. It overrides sensor noise,
which this retained-state argument does not identify with the wrapper defaults.
Before Live, MEKF prediction/measurement is disabled; constructor zero means,
gyro-bias variance RN32(1e-6), linear variances 1/400/2500 and the rounded square
of the 0.004 accel-bias standard deviation persist. The final a_w block comes
from the actual goLive commit, so no unused initial noise value is substituted.

For a construction-rooted ungauged startup, the canonical CORE constructor
consumes the carried machine private observer directly. It executes both
quaternion normalizations in the zero-heel setter, derives and normalizes the
world-down covariance axis, executes the covariance fallback if needed, and
builds all joint24 coordinates from the same physical reference and literal
zero estimator means. No replacement quaternion, reciprocal or covariance is
accepted. The resulting initialization operation ledger is retained across
the Live events. Its scalar Eigen reduction/contraction schedule is explicit;
  target compilation correspondence remains required.

The outer product separately carries the already implemented compiler-specific
WPE, frontend, TuneState, active parameters, scheduler, a_w and Racc histories.
The reused CORE/control runtime receives the actual machine private observer;
other tuner fields in its exact-event template are not consumed by MAG/HOLD.
Source ordinals, physical endpoints and the persistent bias-family history stay
owned by the source interleaver. Completion requires 600 actual attached
compiler events. The older diagnostic-only `observe_*` entries record missing
successors and cannot satisfy completion.

This repairs successor substitution for the conditional finite event graphs.
It does not prove that all source inputs have valid witnesses or that the
finite-real CORE algorithms equal rounded Eigen execution. The following
obligations remain explicit:

- The gauged startup path still lacks machine magnetic pending-yaw ancestry.
  It retains the conditional exact q_hat/joint24/non-a_w covariance premise
  and the separate machine a_w write. Equality of the existing exact `entry`
  and `fresh` objects does not establish machine initialization. The ungauged
  producer does not resolve this different branch or prove startup capture.
  Component startup fixtures without the sealed construction root also remain
  conditional. The default constructor proof does not certify arbitrary public
  constructor settings or state/covariance setters.
- The initialized observer's existing upper norm bound does not alone prove
  the handoff's strict quaternion-norm lower check. The new producer executes
  that check from the actual quaternion; a universal noncancellation/lower-norm
  argument remains required, as does source/control reachability of goLive.
- Prediction, matrix products, solvers, accepted measurement updates, reset
  transcendental operations and gravity/magnetic arithmetic need their actual
  target rounding and totality bounds.
- API vector rounding is now attached in prediction, private Mahony, guarded
  accel and ungauged gravity. Subsequent de-heel, gyro-bias subtraction and
  gravity rotation/filter arithmetic are still conditional finite-real graphs;
  attaching their inputs does not qualify target rounding of their operations.
- The user-authorized raw MEMS cap does not itself bound the MEKF's accumulated
  gyro-bias estimate or prove all-event finite arithmetic over 600 updates.

This composition supplies no storage, rho estimate, startup theorem or new source
restriction.  Both proof tracks and all final ALT gates retain their scope.
