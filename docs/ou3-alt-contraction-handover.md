# OU-III parallel ALT proof handover

## Resume point after PR #523

PR #523 is the completed handoff checkpoint for this ALT continuation and is intended to be merged into `main`. A new conversation should start a **new PR from latest `main`**, not continue the old branch.

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
10. `docs/ou3-brmm-main-handover.md` for the independent original proof track.

The original P2/P3/P4/P5 route remains independently continuable. Do not weaken, delete, rewrite, or make it depend on ALT. `P3=1e-18` remains frozen.

## Certified ALT deployment scope: wind heel excluded

ALT explicitly excludes the optional wind-heel retarget feature. Certified histories require

- `wind_heel_rad_ == 0` from construction onward; and
- zero calls to `update_wind_heel()`.

Shipping is unchanged and already defaults heel to zero. Hence B'=B throughout the certified history, de-heeling is identity, and no wind-heel/body-frame retarget event belongs to the ALT hybrid language. `deployment_scope.py` and the storage guards fail closed if this scope is lost. Do not spend ALT proof effort on wind heel unless the theorem scope is explicitly widened in a future PR.

Retain joint24 `z=(c,e_bg,e_v,e_p,e_S,e_aw,e_ba,beta)`, full 21-state covariance, all motion/bias cross terms, one persistent Live S origin and corrected COMPLETE-BRMM including acceleration `<=8.8 m/s^2` and all-time centered primitive `D_S<=1100 m*s`. Preserve same-signal WPE/bandpass/sigma/tau/T_S ancestry, staged tuner state and separate BIAS0/BIAS1/BIAS2 histories.

## What PR #523 closed

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

`MAG-CALL-SCHEDULE-v1` requires the first post-Live magnetometer call within 40 ms and subsequent gaps <=40 ms. Shipping counts attempted post-delay `updateMag()` calls independently of innovation acceptance, so the 250-count internal unlock guard clears within 10 s of gauged Live and the strict `>1 s` condition is automatically satisfied.

Do **not** assert eventual A21 for arbitrary external hold histories. The theorem hybrid language contains both branches:

- no hold: qualified unlock applies the exact H18->A21 BA variance-floor edge;
- hold active: internal lock clears but H18 may persist indefinitely;
- first later hold release while Live applies the same exact H18->A21 edge;
- asserting hold from A21 returns to H18 and zeros BA cross-covariances exactly as shipping does.

## Next proof target

The decisive remaining object is still the **complete source-uniform finite 600-step shipping word**. The next conversation should continue from latest `main` and attack, in this order:

1. compose successive Live prefixes with the literal asynchronous `updateMag()` interleave, including accepted/rejected magnetic measurements, reference refinement/generation changes and the wrapper count/time bookkeeping on the same physical history;
2. compose the firing Live tilt-reset branch (`initialize_from_acc_preserve_yaw`) into that same hybrid word instead of leaving it as `NotImplemented`/fail-closed;
3. finish source qualification for the remaining continuous-hard-iron/refinement paths actually retained in the declared theorem scope;
4. close deployment arithmetic: startup `atan2`, AngleAxis, quaternion normalization, handoff setters, binary32/binary64 clocks, exp/sqrt/acos/hypot/libm, Eigen solver/PSD/eigensolver and nonfinite branches;
5. bind corrected COMPLETE-BRMM plus BIAS0/1/2 and all declared sensor/model/finite-precision disturbances through every literal prefix;
6. only after `assert_finite_storage_master` accepts the complete source-uniform word, begin the high-precision feasibility/common joint24 storage search, then prove every-prefix retention and the disturbance-dependent ultimate bound.

Do not replace any of these with traces, more random seeds, frozen gains, independently boxed coefficients, covariance consistency, Jacobian products advertised as finite maps, or a convenient fresh-entry set.

## CI state at handoff

At the final pre-merge checkpoint, PR #523 was mergeable. The newest `ou3-alt-contraction` run for the final handover head had been queued and had not yet produced a result. Earlier focused **Finite physical identities (not stability qualification)** runs had passed exact algebra/anti-promotion checks, native shipping correspondence, rank-three benchmark and fail-closed gates.

Several independent original-proof/source workflows remained red on their already-open fail-closed blockers (including the continuous Mahony/source prerequisites). PR #523 did not weaken or bypass them. Do not interpret those failures as a new ALT runtime-algebra regression, and do not make them green by weakening the original proof.

## Gate state

`ALT_LIVE_PASS=false`

`ALT_STARTUP_PASS=false`

`ALT_END_TO_END_PASS=false`

Original P4/P5 remain independent and unpromoted. No certified common storage, `rho`, retained basin, ultimate bound, or end-to-end theorem is claimed at this handoff.
