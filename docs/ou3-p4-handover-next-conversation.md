# OU-III P4 canonical handoff for the next conversation

## Status

P4 is **not proved**. Keep all of the following false until the complete theorem family below closes:

- `P4_MOTION_PASS = false`
- `P4_PASS = false`
- `P5_MAY_START = false`

Do not weaken the deployed filter, hard entry domain, quality gates, actual-applied `R_S` semantics, or canonical P3 `delta = 1e-18` to obtain closure.

## Permanent theorem-family contract

From this merge forward, a result may be considered a P4 proof **only** if it is source-uniform over the deployed `COMPLETE_BRMM_NORMAL_LIVE_WORD` and jointly retains every item below through endpoint and literal-every-prefix augmented LDLT:

1. **BRMM only.** The physical source must be the complete admitted BRMM source/history. Pinned traces, replay, finite harmonic substitutes, independent coefficient boxes, or simplified source families cannot promote P4.
2. **All accelerometer-bias families:** `BIAS0`, `BIAS1`, and `BIAS2`. Each is a separate mandatory physical-driver family. A proof for one family does not imply either of the others.
3. **Signal-derived tau and sigma.** The same admitted IMU/physical signal history must drive the WavePeriodEstimator and the period-scaled sigma estimator. `tau` and `sigma` may not be selected independently or swept as unrelated theorem coordinates.
4. **Interdependent adaptive map.** Preserve the deployed dependency chain

   `same signal -> period/frequency -> tau`

   `same signal + same frequency -> raw sigma -> effective sigma`

   `tau -> T_S`

   `(tau, sigma, T_S) -> R_S`

   including the actual anisotropic `R_S` used by each `S=0` Joseph event.
5. **Adaptation schedule.** Preserve candidate/active EMA state, staged commit behavior, scheduler elapsed state, commit timing, and the exact active schedule used by each Riccati/Joseph event. A proof that replaces this with arbitrary per-event `(tau,sigma,T_S,R_S)` choices is invalid.
6. **Same-history Kalman/reset graph.** Reachable shipping `P`, event geometry, measurement `R`, `K`, reset transport, tuner state, scheduler state, physical bias state/driver, and projection must remain tied to the same source branch/history.
7. **Finite precision.** Include the declared binary32 Kalman/reset/platform arithmetic enclosure as an additive ISS channel.
8. **Full hard entry set and every prefix.** Closure must hold on every admitted hard-entry radial segment, at the endpoint, and after every literal event prefix with compatible storage/LDLT retention.

The final fail-closed promotion gate now encodes this contract. In particular it requires the exact bias-family set `['BIAS0','BIAS1','BIAS2']`, COMPLETE BRMM, the same-signal adaptive source contract, the interdependent `(tau,sigma,T_S)->R_S` relation, and the adaptation/commit/scheduler history before `P4_MOTION_PASS` can ever become true.

## What is already available

The merged branch contains substantial reusable machinery: exact finite-angle chord/projection sectors; same-cell Joseph/reset signed identities; hard entry-set admission; BIAS1 conditional root/driver recurrence; compatible storage infrastructure; source-cover/event lineage APIs; joint WavePeriodEstimator/AdaptiveWaveBandPass estimator plumbing; BRMM attachment machinery; actual-applied `R_S` provenance; and binary32 additive ISS infrastructure.

These are prerequisites and partial lemmas, not P4 closure.

## Mandatory next proof work

The next conversation should continue only along this path:

1. Materialize source-uniform `BIAS0` and `BIAS2` physical-driver recurrences alongside the existing `BIAS1` recurrence; do not infer them from BIAS1.
2. Complete the same-signal BRMM estimator cover so the input history jointly generates frequency, tau, raw/effective sigma, `T_S`, and SpectralMSE `R_S` while preserving candidate/active EMA, staged commits, and scheduler semantics.
3. Emit correlated same-history source cells for every admitted BRMM continuation, bias family, and hard-entry radial segment.
4. Derive `K` only from the same reachable `P/H/R` cell and retain exact reset transport/projection coupling.
5. Form source-correlated endpoint and literal-every-prefix augmented matrices for all three bias families.
6. Run outward LDLT with finite-precision ISS charges and prove the full retained-coordinate budgets.
7. Only after **all three bias families plus the adaptive BRMM source cover** close may the final gate be changed to report their closure bits true.

## Invalid routes

Do not claim P4 from any of the following:

- BRMM replay or a single captured trajectory;
- BIAS1-only, BIAS0-only, or BIAS2-only closure;
- independent `(f,sigma)` or `(tau,sigma,T_S,R_S)` rectangles;
- oracle/fixed tuner coefficients detached from the input signal;
- an adaptation schedule chosen independently of the shipping EMA/commit/scheduler state;
- detached `K`/`P`/`H`/`R` boxes;
- endpoint-only contraction without every-prefix retention;
- covariance consistency substituted for the declared hard physical entry set;
- floating-point behavior omitted from the final theorem.

## Suggested prompt for the next conversation

> Continue OU-III P4 from current `main` using `docs/ou3-p4-handover-next-conversation.md` as the canonical contract. Only COMPLETE BRMM proofs covering BIAS0, BIAS1, and BIAS2 are admissible. Preserve the same input-signal ancestry that estimates tau and sigma, the dependent `tau -> T_S` and `(tau,sigma,T_S) -> R_S` mapping, candidate/active EMA + staged commit + scheduler behavior, same-history reachable Kalman/reset coefficients, and finite precision. Close endpoint and literal-every-prefix augmented LDLT for all three bias families. Keep P4/P5 false until the final fail-closed gate genuinely closes.
