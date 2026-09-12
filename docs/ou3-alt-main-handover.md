# OU-III ALT continuation entry point

Read `AGENTS.md`, the normative [`ou3-alt-proof-plan.md`](ou3-alt-proof-plan.md),
[`ou3-alt-contraction-handover.md`](ou3-alt-contraction-handover.md), the ALT
section of `ou3-proof-research-state.md`, and
[`ou3-alt-live-magnetic-word.md`](ou3-alt-live-magnetic-word.md).
Continue an explicitly named active PR; after merge, start a new PR from latest
`main`. The original P2/P3/P4/P5 route remains independently continuable through
`ou3-brmm-main-handover.md`; `P3=1e-18` is frozen.

## Current checkpoint

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

Do not start storage/rho search. The primary physical-history -> finite-word
restriction is now closed without a common generator assumption. The immediate
blocker is to bind the remaining physical frequency relation, full BIAS
generating histories, and finite estimator-owned coefficient/product graphs
(frontend/tuner/guard/Racc, prediction Q/F, S service, accel/mag
innovation/K/Joseph/reset decisions) to that same admitted history at every
literal branch. The magnetic sample-zero endpoint and all later represented
magnetic endpoints already consume the same carried root/reference object. Then
close startup/ungauged and deployment-arithmetic paths.

For deployment correspondence, finish sqrt/acos/asin/atan2/AngleAxis and
normalization/cutoff/nonfinite branches, floating clocks, Eigen decisions, and
the finite-width signed magnetic counter. Only after the finite product graph is
source-uniform and `assert_finite_storage_master` accepts it may common joint24
storage search begin, followed by every-prefix retention and the ultimate bound.
Traces, random seeds, frozen gains, independent boxes, covariance-consistency
entry assumptions and convenient error sets do not replace these obligations.

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.
No certified rho, storage basin, ultimate bound or end-to-end theorem is claimed.
