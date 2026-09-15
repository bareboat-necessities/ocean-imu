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

- The initial compiler CORE is conditionally attached to the startup runtime;
  a target-qualified startup CORE/covariance correspondence is still required.
- Prediction, matrix products, solvers, accepted measurement updates, reset
  transcendental operations and gravity/magnetic arithmetic need their actual
  target rounding and totality bounds.
- The predictor and ungauged gravity path still read the source's exact raw
  packet representation. Their binary32 API projection must be accounted for
  in target correspondence; the private Mahony and guarded accel paths already
  consume their explicit machine operands.
- The user-authorized raw MEMS cap does not itself bound the MEKF's accumulated
  gyro-bias estimate or prove all-event finite arithmetic over 600 updates.

This composition supplies no storage, rho estimate, startup theorem or new source
restriction.  Both proof tracks and all final ALT gates retain their scope.
