# OU-III parallel ALT proof handover

## Current continuation checkpoint

Continue an explicitly named active PR. After it is merged, sync latest `main`
and create a new PR from that head. Do not replace the independent proof track.

Read, in this order:

1. `AGENTS.md`;
2. `docs/ou3-alt-main-handover.md`;
3. this file;
4. `docs/ou3-alt-proof-plan.md`;
5. ALT section of `docs/ou3-proof-research-state.md`;
6. `docs/ou3-alt-finite-measurement-proof.md`;
7. `docs/ou3-alt-finite-core-composition.md`;
8. `docs/ou3-alt-runtime-primitives.md`;
9. `docs/ou3-alt-mahony-binary32.md`;
10. `docs/ou3-alt-live-magnetic-word.md`;
11. `docs/ou3-brmm-main-handover.md` for the independent original proof track.

The original P2/P3/P4/P5 route remains independently continuable. Do not weaken, delete, rewrite, or make it depend on ALT. `P3=1e-18` remains frozen.

## Certified ALT deployment scope: wind heel excluded

ALT explicitly excludes the optional wind-heel retarget feature. Certified histories require

- `wind_heel_rad_ == 0` from construction onward; and
- zero calls to `update_wind_heel()`.

Shipping is unchanged and already defaults heel to zero. Hence B'=B throughout the certified history, de-heeling is identity, and no wind-heel/body-frame retarget event belongs to the ALT hybrid language. `deployment_scope.py` and the storage guards fail closed if this scope is lost. Do not spend ALT proof effort on wind heel unless the theorem scope is explicitly widened in a future PR.

Retain joint24 `z=(c,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta)`, full 21-state covariance, all motion/bias cross terms, one persistent Live S origin and corrected COMPLETE-BRMM including acceleration `<=8.8 m/s^2` and all-time centered primitive `D_S<=1100 m*s`. Preserve same-signal WPE/bandpass/sigma/tau/T_S ancestry, staged tuner state and separate BIAS0/BIAS1/BIAS2 histories.

## Retained finite runtime and startup relations

### Finite shipping runtime algebra

The ALT graph now contains exact/local finite descriptors for:

- physical attitude and translation prediction with same-history continuous physical increments;
- full 21-state covariance prediction with runtime-generated attitude F/Q, literal integrated-OU Qaxis branches and BA decay;
- pending a_w covariance floor, scheduler/S service and covariance hygiene branches;
- inverse-free measurement relations `K Sigma = N`, SafeLDLT first-attempt/retry/reject, Joseph/reset/projection;
- accelerometer vibration guard, applied R_acc and held-sample physical forcing;
- private Mahony, WPE, adaptive band/statistics, stillness and staged tuner commit;
- asynchronous magnetometer control and H18/A21 wrapper logic;
- first ordinary Live IMU prefix.

The initialized private-Mahony path also has a named exact binary32 graph under its stated arithmetic profile. Target compiler/libm/profile qualification remains open.

### Deterministic magnetic startup source

`MAG-BMM150-DET-v1` is the theorem-facing commissioned-installation source class:

- `20 <= ||B_W|| <= 75 uT`;
- horizontal field `>=15 uT`;
- body hard iron `<=5 uT`;
- deterministic residual `<=2 uT` per theorem sample.

Together with the declared startup gravity-direction error `<=0.02 rad`, the real-arithmetic startup argument gives total horizontal perturbation `<=8.5 uT`, `|sin(delta_yaw)| <= 17/30`, yaw error `<0.61 rad`, and total attitude error `<0.63 rad < pi/4`. No statistical `1/sqrt(N)` reduction is used.

### Exact fresh H18 joint24 entry

The gauged zero-heel handoff is composed through shipping `goLive` / `initialize_from_attitude` / `enterLive_`.

Important frame convention: shipping receives the accepted handoff as boat-to-world `q_BW`, while the internal MEKF nominal attitude is world-to-body `q_WB = conjugate(q_BW)` under the zero-heel scope. A nontrivial quaternion regression exists specifically to prevent accidentally reversing this convention.

`finite_startup_live_entry.py` seats `P_aw,aw` on the same committed `Sigma_aw`, clears every a_w cross-covariance, requires the committed Live `R_S`, and keeps accelerometer-bias learning disabled.

`finite_fresh_joint24_entry.py` derives the actual fresh state from the SAME physical `Reference` and estimator coordinates:

`c = Cayley(q_true_WB * conjugate(q_hat_WB))`,
`e_bg=b_g-b_g_hat`, `e_v=v-v_hat`, `e_p=p-p_hat`,
`e_S=S_centered-S_hat`, `e_aw=a-a_w_hat`,
`e_ba=beta-b_a_hat`, with final joint24 coordinates equal to the same true `beta`.

The one-time Live origin and fresh centered physical S=0 are enforced. There is no assumed covariance-consistency entry set and no independently chosen fresh-entry error box.

### Startup frontend memory is preserved through goLive

`finite_startup_live_runtime_bridge.py` closes the control/memory bridge from TunerReady to Live:

- Mahony/WPE/band/statistics/stillness/vibration-guard memory is preserved;
- only startup stage/clock changes at handoff;
- active tau/Sigma_aw/pseudo cadence/Live R_S are derived from the SAME carried `TuneState`;
- the persistent S scheduler is retargeted to the committed pseudo period;
- an online pending-tune bit is intentionally preserved into the first Live IMU boundary, matching shipping;
- startup is forbidden from carrying a periodic Live-only a_w floor request.

### First actual Live sample is now rooted in startup

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

The covariance ordering mirrors shipping: `initialize_from_acc()` first installs
the accel-only qref and calls `set_accel_only_attitude_covariance_()`, so the
anisotropic yaw axis is world-down expressed in that intermediate body frame.
Only afterwards does `initialize_from_acc_preserve_yaw()` restore yaw with
`set_quaternion_boat(q_new_bw)`. That setter changes qref and zeros attitude
error bookkeeping but does not rotate or reseed covariance. The finite graph
carries the intermediate axis explicitly. A nontrivial rational regression
checks the 24-7-25 tilt plus retained yaw and the resulting exact covariance.

The exact-real accelerometer branches now also use shipping's strict cutoffs:
`anorm < 1e-8` rejects initialization and `norm_axis < 1e-8` selects the
near-parallel/anti-parallel quaternion branch. The latter is represented by the
squared nonnegative comparison `axis_norm^2 < 1e-16`, avoiding a free sqrt.
Binary32 sqrt/comparison correspondence at these cutoffs remains open.

`finite_live_tilt_prefix.step_from_shipping_operands` owns this edge, and
`finite_live_interleave.imu_step` rejects attempts to inject `tilt_deg` or an
independent reset output. The >70-degree comparison is resolved with a rigorous
rational enclosure of the real threshold corresponding to
`acos(cos_tilt)*57.295779513f`; the tiny enclosure boundary remains fail-closed.
Nonfinite fallback, sqrt/acos/asin/atan2/AngleAxis/normalization rounding and
target compiler/libm correspondence remain deployment obligations. They are not
converted into theorem exclusions.

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
