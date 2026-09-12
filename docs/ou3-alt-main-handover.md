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
The exact same-history magnetic discrepancy is retained for storage work.

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
not rotate or reseed P. Deployment libm/normalization/cutoff/nonfinite details
remain fail-closed.

## Next work

Do not start storage/rho search. Finish the remaining literal IMU/source/BIAS,
startup/ungauged and deployment-arithmetic paths in the complete source-uniform
600-step word. For the tilt edge, close the near-parallel cutoff and target
binary32/libm correspondence rather than replacing them with an ideal branch.
Close remaining libm/casts/Eigen/nonfinite branches and lifetime clock/counter
arithmetic, including the shipping signed magnetic count.

Only after `assert_finite_storage_master` accepts that complete representation
may high-precision/common joint24 storage search begin, followed by every-prefix
retention and the ultimate bound. Traces, random seeds, frozen gains, independent
boxes, covariance-consistency entry assumptions and convenient error sets do
not replace these obligations.

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.
No certified rho, storage basin, ultimate bound or end-to-end theorem is claimed.
