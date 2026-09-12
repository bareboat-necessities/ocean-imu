# OU-III parallel ALT proof handover

## Current continuation checkpoint

Continue an explicitly named active PR. After it is merged, sync latest `main`
and create a new PR from that head. Do not replace the independent proof track.

Read, in this order:

1. `AGENTS.md`;
2. `docs/ou3-alt-proof-plan.md`;
3. this handover;
4. ALT section of `docs/ou3-proof-research-state.md`;
5. `docs/ou3-alt-contraction.md`;
6. `docs/ou3-alt-finite-measurement-proof.md`;
7. `docs/ou3-alt-finite-core-composition.md`;
8. `docs/ou3-alt-runtime-primitives.md`;
9. `docs/ou3-alt-mahony-binary32.md`;
10. `docs/ou3-alt-live-magnetic-word.md`;
11. `docs/ou3-brmm-main-handover.md` for the independent original proof track.

The original P2/P3/P4/P5 route remains independently continuable. Do not weaken, delete, rewrite, or make it depend on ALT. `P3=1e-18` remains frozen.

## Immutable ALT scope

ALT is the zero-wind-heel theorem branch. `wind_heel_rad_==0` from construction onward; no `update_wind_heel()` event is admitted. Shipping already defaults to zero heel, so this is not a filter change. Every theorem-facing IMU event also rejects a nonidentity de-heel map.

Retain joint24 `z=(c,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta)`, full 21-state covariance, all motion/bias cross terms, one persistent Live S origin and corrected COMPLETE-BRMM including acceleration `<=8.8 m/s^2` and all-time centered primitive `D_S<=1100 m*s`. Preserve same-signal WPE/bandpass/sigma/tau/T_S ancestry, staged tuner state and separate BIAS0/BIAS1/BIAS2 histories.

## Retained finite runtime and startup relations

The conditional finite graph already materializes prediction, covariance prediction, accel/mag/S measurement updates, Joseph/reset relations, accel-bias projection, periodic a_w covariance synchronization, held-bias semantics, runtime tuner application, frontend/private-Mahony/WPE/band/sigma/stillness state, startup magnetic acquisition/gauge, hard startup handoff, the first actual Live sample, H18/A21 control, and successive IMU/magnetic/hold substitution. These are finite identities and conditional relations, not source-uniform certification.

### First actual Live sample is rooted in startup

`finite_startup_first_live_step.py` removes the synthetic Live root at the first represented sample. The first prediction/S-service/accelerometer/tuner-WPE prefix is invoked directly from the startup-produced H18 state, with the raw sensor packet and physical segment required to start at exactly that same fresh physical endpoint.

This closes the former event-boundary gap between startup and the first represented Live IMU sample. Negative tests use a valid sensor packet at the wrong physical endpoint so the same-history ancestry guard itself is exercised.

### H18/A21 control semantics

`MAG-CALL-SCHEDULE-v1` requires the first post-Live magnetometer call within 40 ms and subsequent gaps <=40 ms. Shipping counts attempted post-delay `updateMag()` calls independently of innovation acceptance, so the 250-count internal unlock guard clears within 10 s of gauged Live and the strict `>1 s` condition is met by some call within that deadline.
The 250th call itself can occur too early; see the corrected timing proof in
`ou3-alt-live-magnetic-word.md`.

Do **not** assert eventual A21 for arbitrary external hold histories. The theorem hybrid language contains both branches:

- no hold: qualified unlock applies the exact H18->A21 BA variance-floor edge;
- hold active: internal lock clears but H18 may persist indefinitely;
- first later hold release while Live applies the same exact H18->A21 edge;
- asserting hold from A21 returns to H18 and zeros BA cross-covariances exactly as shipping does.

## Successive Live magnetic composition

`finite_live_interleave.py` carries the startup bridge, tilt watchdog, magnetic
calibration/control and schedule-prefix state through successive IMU, magnetic
and hold events. Every magnetic call reads the current physical endpoint and
private observer directly from that product state. Its full covariance/mean
successor feeds the next IMU; calibration, tuner, clocks and Live S origin cannot
be restarted at a represented boundary.

`finite_live_magnetic_word.py` composes source qualification, raw continuous
accumulation, acquisition refinement, same-mean reference/yaw writes, hold
release, continuous hard-iron application and the same corrected measurement in
shipping order. The actual default continuous estimator remains enabled. Its
memory begins before startup gravity admission, persists through goLive, and
supplies both its fit and reference correction from one sufficient-statistic
state. Both refinement and continuous reference writes can occur in one call.

`finite_continuous_mag_runtime.py` retains finite-real accumulator, information,
regularized solve, bias/residual gates, and failed-application clock/anchor
side effects. The supplying note proves the weighted-rotation variance identity
and finite-real all-history bounds: raw <=82 uT, accepted/applied continuous
offset <=28.7 uT, active reference/corrected sample <=110.7 uT. The resulting
coarse discrepancy ceiling is 221.4 uT, NOT an asserted useful ISS margin.
The exact vector discrepancy remains in the event graph.

The timing proof uses `a+max((n-1)g,1+g)` rather than claiming that the n-th
call itself must satisfy the strict one-second guard. The default internal
unlock deadline remains ten seconds, under locally finite calls continuing
through unbounded physical time. No eventual A21 is inferred under arbitrary
external hold. Finite prefix timing checks do not establish infinite coverage.

## Live tilt reset is now same-operand finite-real

`finite_tilt_reset_runtime.py` removes the free watchdog angle and free final
preserve-yaw quaternion from the theorem-facing Live edge. The watchdog predicate
is derived from the actual post-accelerometer nominal attitude. A firing reset
derives the old yaw from that same predecessor, derives accel-only tilt from the
same guarded accelerometer already consumed by the sample, reconstructs shipping
pitch/roll, and composes yaw-pitch-roll.

The covariance order now mirrors shipping exactly. `initialize_from_acc()` first
installs the accel-only qref and immediately calls
`set_accel_only_attitude_covariance_()`, so its anisotropic yaw axis is world-down
expressed in that accel-only intermediate body frame. Only afterwards does
`initialize_from_acc_preserve_yaw()` build the yaw-restored quaternion and call
`set_quaternion_boat()`. That setter changes qref and zeros the attitude-error
bookkeeping but does **not** rotate or reseed covariance. The finite graph carries
the intermediate covariance axis explicitly; a rational 24-7-25 regression also
checks that the ideal yaw restore leaves the gravity axis unchanged while
preserving the actual shipping evaluation order for later binary32 proof.

`finite_live_tilt_prefix.step_from_shipping_operands` owns this edge, and
`finite_live_interleave.imu_step` rejects attempts to inject `tilt_deg` or an
independent reset output. The >70-degree comparison is resolved with a rigorous
rational enclosure of the real threshold corresponding to
`acos(cos_tilt)*57.295779513f`; the tiny enclosure boundary remains fail-closed.
Near-parallel accelerometer cutoff behavior, nonfinite fallback, sqrt/acos/asin/
atan2/AngleAxis/normalization rounding and target compiler/libm correspondence
remain deployment obligations. They are not converted into theorem exclusions.

## Next proof target

The decisive missing object remains the complete source-uniform finite 600-step
shipping word. Continue from the product composer, not disconnected snapshots:

1. close deployment correspondence for the now-bound Live tilt reset and startup
   attitude/magnetic arithmetic, including normalization cutoffs, libm,
   float/double casts, Eigen decisions and nonfinite outcomes;
2. finish the remaining startup/ungauged and magnetic source/history paths;
3. attach corrected COMPLETE-BRMM, BIAS0/1/2 and all declared disturbances to
   every literal IMU/magnetic/hybrid prefix and deployment lifetime arithmetic;
4. only after the finite-master guard accepts the complete source-uniform
   relation, attempt common joint24 storage, every-prefix retention and ultimate
   bounds.

The shipping `int mag_updates_applied_` keeps incrementing after unlock. Its
finite-width lifetime behavior and floating clocks are explicit open deployment
obligations. Do not replace either with an unbounded Python arithmetic claim.
Do not silently disable calibration or shrink the physical domain to get PASS.

## Validation boundary

The pre-tilt continuation focused finite-identity selection had passed 470 tests,
including calibration equations, rejection side effects, source/clock continuity,
yaw writes, H18/A21 control and IMU->mag->IMU covariance substitution. The
focused workflow includes `test_finite_tilt_reset_runtime`; latest-head CI must
be read separately before claiming it passes. Native shipping correspondence on
the prior tested tree passed with bit-identical observed/plain sample states and
unchanged tracked headers. These are implementation/algebra checks, not
source-uniform stability evidence.

The full inherited suite is known not to be green on the shared baseline because
of the existing continuous-Mahony/source prerequisites. Do not weaken those
prerequisites to make ALT green. The prior local default `make all` environment
also lacked `/usr/include/eigen3`; no full-build PASS is inferred from a focused
native correspondence check.

## Gate state

`ALT_LIVE_PASS=false`

`ALT_STARTUP_PASS=false`

`ALT_END_TO_END_PASS=false`

Original P4/P5 remain independent and unpromoted. No certified common storage, `rho`, retained basin, ultimate bound, or end-to-end theorem is claimed at this handoff.
