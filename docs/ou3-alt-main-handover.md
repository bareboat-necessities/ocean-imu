# OU-III ALT continuation entry point

Read `AGENTS.md`, the normative [`ou3-alt-proof-plan.md`](ou3-alt-proof-plan.md),
[`ou3-alt-contraction-handover.md`](ou3-alt-contraction-handover.md), the ALT
section of `ou3-proof-research-state.md`, and
[`ou3-alt-live-magnetic-word.md`](ou3-alt-live-magnetic-word.md).
Continue an explicitly named active PR; after merge, start a new PR from latest
`main`. The original P2/P3/P4/P5 route remains independently continuable through
`ou3-brmm-main-handover.md`; `P3=1e-18` is frozen.

## Current checkpoint

The strong startup/goLive/Live product carries a persistent machine adaptive
band and DebiasedEMA statistics state for each tuner compiler track. Cold
samples advance those states although tau/sigma/R_S candidates remain held;
MAG/HOLD and goLive preserve them by identity. Live entry consumes the startup
histories, not a newly synthesized band/statistics ledger.

At sample entry, a pending commit reads each carried band's p11 and the same
compiled bench sigma. After the band step, the new band output drives the
statistics update, and the new variance/noise readouts drive sigma. The two
readout times must not be exchanged. The scalar transaction still uses
`finite_tuner_machine_boundary_commit` and retains rounded tau, cadence,
stationary covariance squares and R_S. Ordinary boundaries clear pending;
goLive commits unconditionally but preserves the original pending bit.

The frequency graph retains both distinct clamps:

`WPE/prior -> statistics clamp/store -> outer tuning clamp -> tau target`.

Band corners instead read the previous statistics frequency, falling back to
the sample-entry WPE/prior when unavailable. Exact bounds remain distinct from
their compiled binary32 values. The statistics variance readout retains both
separate multiply/subtract and contracted FMA possibilities. Ready p11=0 has
sqrt=0; a zero bench sigma does not waive sqrt ancestry checks.

Prediction roots derive the same pending transaction that the Live prefix
executes before prediction. They may not use stale pre-boundary tau/Sigma.
Executed regressions carry a changed pending covariance through both machine
prediction-root families, the per-track scheduler and the full joint24/21x21
prediction-displacement relation. The higher wrappers count the actual lower
IMU ordinal. No second physical segment or Live-origin reset is introduced.

The 938-test finite-map selection passes locally. Native build validation is
blocked by missing Eigen headers in this environment; this is not a full CI or
native correspondence PASS. The represented arithmetic relations still take
conditional libm/contraction witnesses. Machine vertical-input and stillness
production, compiler/nonfinite correspondence, source-uniform supply bounds,
and post-prediction S/accelerometer propagation remain open. No storage search
or stability promotion is authorized.

ALT excludes wind heel: `wind_heel_rad_==0` from construction onward and no
`update_wind_heel()` events. Shipping is unchanged. The product composer also
rejects nonidentity IMU de-heel maps at every represented boundary.

The finite graph joins the exact gauged H18 startup bridge to successive IMU,
asynchronous magnetic and hold events. The magnetic path includes raw physical
source qualification, the default continuous hard-iron estimator,
refinement/reset clocks, same-mean reference and yaw writes, coupled offset /
reference application, full measurement/covariance successor, and count/hold
logic. Continuous magnetic memory starts before startup admission and is not
restarted at Live. All full-21 covariance blocks and joint24 coordinates persist.

The finite-real calibration lemma gives raw norm <=82 uT, accepted/applied
continuous offset <=28.7 uT and active reference/corrected observation <=110.7
uT. These are boundedness results, not contraction or calibration accuracy.
The exact same-history magnetic discrepancy is retained for later storage work.

The corrected timing theorem proves internal unlock by ten seconds on the
existing default call schedule. Call 250 need not satisfy the strict one-second
guard; a later call can do so. External hold may keep H18 indefinitely.
Finite-prefix deadline checks never certify infinite coverage.

The theorem-facing firing tilt-watchdog edge no longer accepts a free tilt angle
or final reset quaternion. Its >70-degree predicate is derived from the same
post-accelerometer nominal attitude, and preserve-yaw reconstruction is tied to
the same predecessor and guarded accelerometer. The covariance reseed follows
shipping order exactly: the anisotropic yaw axis is taken from the accel-only
intermediate qref before yaw is restored; the later `set_quaternion_boat()` does
not rotate or reseed P. The exact-real `anorm < 1e-8` and `norm_axis < 1e-8`
branches are materialized. Native shipping regression now exercises both the
24-7-25 retained-yaw reset and a nonzero sub-1e-8 near-parallel axis; universal
binary32/libm/cutoff/nonfinite correspondence remains fail-closed.

`finite_source_bound_live_word.py` carries the source continuation, concrete
bias history, sensor history identities and persistent runtime configuration
alongside the actual Live product. Its IMU call appends the next source segment
before running the filter; magnetic and hold calls preserve the same endpoint.
The inherited source-product fixture now uses the actual `J0/J1/J2` fields.

`finite_source_continuation.py` checks physical p/v/a/S vector caps, the coupled
three-axis moment IQC, a necessary rotation chord bound, physical sampled rate,
and the analytic BIAS envelopes. It carries the same actual bias factor over
consecutive 5 ms segments and preserves exact outward contract endpoints.
`finite_brmm_moment_prefix.py` derives the prefix moments and their energy budget
by exact Gramian/kernel concatenation, without independent source boxes or S
restarts. See [`ou3-alt-source-continuation.md`](ou3-alt-source-continuation.md)
for the real-arithmetic projection/induction proof and its limits.

The legacy `Qualified*` names do NOT make arbitrary runtime objects admitted
COMPLETE-BRMM/BIAS histories. Source labels are ancestry, not admission. The
fresh Live origin itself is checked directly against the carried root, physical
vector caps and BIAS envelope, so an asynchronous magnetic call may occur before
IMU transition 1 without fabricating a predecessor segment or advancing source
time/ordinal.

`finite_complete_brmm_restriction.py` now closes the universal physical-source
restriction that this graph needs. COMPLETE-BRMM's primary definition is one
bounded physical p/v/a/attitude history with a uniformly bounded centered
primitive; spectral and shaping-state constructions are only sufficient
certificate methods. Restricting any admitted primary history to the 5 ms grid
therefore yields the finite p/v/a/S caps, one-origin S bound, exact primitive
recurrence, coupled moment IQC and necessary rotation-chord bound. No common
numerical generator representation is required or assumed. Runtime token
matching still does not prove membership of an arbitrary object; the remaining
frequency/frontend, full BIAS-history and deployment-arithmetic relations stay
open.

The older `phase1_closure.py` is now explicitly a source/Jacobian ledger, not a
finite-storage gate. Its old pointwise Jacobian cocycle is useful ancestry but
cannot authorize storage; `assert_finite_storage_master` remains the only gate.
This prevents the old Phase-1 labels from sending ALT into a premature rho/metric
search.

## Next work

Do not start storage/rho search. The next finite-event obligations are to bind
the machine band's vertical input and the sigma stillness operands to the same
persistent guarded/private-Mahony machine history, then propagate the existing
joint24/21x21 prediction displacement through pending covariance floors,
per-track due/not-due S service, actual machine R_S, accelerometer correction
and reset. Retain both track histories and all cross-covariances.

The WPE raw-period/log/exp, band/statistics/tuner libm, Q-axis operation-level
roundoff, Eigen branches and nonfinite/threshold correspondence remain native
deployment obligations. Bound their supplies on the admitted BRMM+BIAS+ISS
family; do not replace same-history products by independent parameter boxes or
use a successful component regression as all-history membership.

Universal startup capture, every-prefix retention, compatible joint24 storage
and an ultimate bound remain open. Indefinite machine execution additionally
has the retained floating-clock stall and unchecked signed magnetic-counter
lifetime obstructions. Tiling finite words must preserve the one-time Live/S
origin and all frontend/covariance/bias history; no wordwise reinitialization.

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`,
`storage_search_allowed=false`.
