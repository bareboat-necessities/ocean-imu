# OU-III ALT continuation entry point

## Start here after merged PR #523

The canonical continuation state is [`ou3-alt-contraction-handover.md`](ou3-alt-contraction-handover.md). A new conversation should sync latest `main` and open a **new PR** from it.

Read this file with `AGENTS.md`, the normative [`ou3-alt-proof-plan.md`](ou3-alt-proof-plan.md), the ALT section of `ou3-proof-research-state.md`, the finite identities in [`ou3-alt-finite-measurement-proof.md`](ou3-alt-finite-measurement-proof.md), [`ou3-alt-finite-core-composition.md`](ou3-alt-finite-core-composition.md), [`ou3-alt-runtime-primitives.md`](ou3-alt-runtime-primitives.md), and [`ou3-alt-mahony-binary32.md`](ou3-alt-mahony-binary32.md).

The original P2/P3/P4/P5 proof route remains separate through [`ou3-brmm-main-handover.md`](ou3-brmm-main-handover.md). Do not weaken it or make it depend on ALT. `P3=1e-18` remains frozen.

## Scope

ALT excludes wind heel. Certified histories require `wind_heel_rad_==0` from construction onward and no `update_wind_heel()` events. Shipping is unchanged. Under this scope B'=B and no wind-heel retarget belongs to the hybrid word.

## Current ALT checkpoint

PR #523 advanced ALT from disconnected finite primitives to an explicit startup-rooted runtime graph:

- physical prediction, full 21-state covariance, runtime attitude F/Q, integrated-OU Qaxis, BA decay, pending a_w floor, S scheduler/service, SafeLDLT branches and Joseph/reset are represented;
- private Mahony, WPE, adaptive band/statistics, stillness, vibration guard/R_acc and staged tuner memory are represented on the same sample ordering;
- deterministic magnetic startup admission/capture is declared and the real-arithmetic startup attitude bound is `<0.63 rad < pi/4`;
- the gauged zero-heel `goLive` handoff produces an exact fresh H18 joint24 `finite_core.State` instead of an assumed error box;
- the TunerReady frontend/tuner/guard memory is preserved through `goLive`, with the same carried `TuneState` supplying active tau/Sigma_aw/pseudo period/Live R_S;
- the first actual Live prediction/S/accelerometer/tuner-WPE prefix now starts directly from that startup-produced H18 state; the former synthetic Live root is removed;
- H18/A21 hybrid control is explicit: internal unlock is forced by the qualified mag-call schedule, while arbitrary external hold may keep H18 indefinitely; release gives the literal H18->A21 edge.

The handoff quaternion convention is important: shipping accepts boat-to-world `q_BW`, while the MEKF internal nominal attitude is world-to-body `q_WB=conjugate(q_BW)` under zero heel. Keep the nontrivial regression that enforces this.

## Next work

Do **not** start storage/rho search yet. The complete source-uniform 600-step shipping word is still missing.

Next, compose successive Live prefixes with asynchronous magnetometer interleaving, reference refinement/generation changes, accepted/rejected mag updates, and wrapper count/time bookkeeping on the same physical history. Then compose the firing Live tilt-reset branch (`initialize_from_acc_preserve_yaw`), finish remaining magnetic/hard-iron source paths retained in scope, and close deployment finite precision/libm/Eigen solver and hygiene branches.

After every literal prefix is attached to corrected COMPLETE-BRMM, BIAS0/1/2 and declared disturbances, allow `assert_finite_storage_master` to decide whether the representation is complete enough for the first high-precision/common joint24 storage search.

Traces, random seeds, frozen gains, independent coefficient boxes, Jacobian products, covariance-consistency assumptions and convenient fresh-entry sets are not substitutes.

## Gate state

`ALT_LIVE_PASS=false`, `ALT_STARTUP_PASS=false`, `ALT_END_TO_END_PASS=false`.

No certified common storage, contraction factor, retained basin, ultimate bound, or end-to-end theorem is available yet.
