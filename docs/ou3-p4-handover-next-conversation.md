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

The merged branch contains substantial reusable machinery: exact finite-angle chord/projection sectors; same-cell Joseph/reset signed identities; hard entry-set admission; compatible storage infrastructure; source-cover/event lineage APIs; joint WavePeriodEstimator/AdaptiveWaveBandPass estimator plumbing; BRMM attachment machinery; actual-applied `R_S` provenance; and binary32 additive ISS infrastructure.

All three physical-driver recurrences are now materialized, each from its own declared parameter box and its own module, and none inferred from another:

- `ou3_p4_bias0_family.py` — composite turn-on/thermal/strain/non-Gauss-Markov driver with a declared pathwise Gauss-Markov increment cap;
- `ou3_p4_bias1_family.py` — one-root one-parameter driver;
- `ou3_p4_bias2_family.py` — bounded-variation drift admitting `phi_true=1`, i.e. no relaxation root;
- `ou3_p4_bias_family_joint_iss_supply.py` — the joint `[e_b;b_true]` supply for all three, each pushed through the deployed 24-state event lift with the shared `w` column retained.

All three now also reach the same-history projection/Joseph graph on the exact-chord bridge, so the **bias-family half of the contract is closed**: `bias_family_source_uniform_same_history_closed` is true for BIAS0, BIAS1 and BIAS2. The lemma that makes this legitimate rather than an inheritance from BIAS1 is that the prerequisites see a family only through four things — admission, one retained physical bias history, the shared `w` column, and `|b_true|` — because the radial projection map `F_R(e,beta)=beta-Pi_R(beta-e)` contains no driver term and the Joseph/reset gains come from the reachable `P/H/R` cell. The three declared boxes share one true-bias envelope (spread 5.6e-17), so the fourth dependence is one number and the compactness bound `|e_b| <= R + B_true = .6252` is the same for all three.

A uniform BIAS2 separation constant `mu_sep` is **not** a prerequisite. The declared objective is bounded bias error plus regional practical ISS of the other 18, and its bias half comes from the closed radial projection sector, which holds for a non-relaxing truth (`phi_true = 1`) exactly as for a relaxing one. `mu_sep` would only sharpen the motion-channel gains; it stays unproved and is reported as an optional sharpener.

BIAS0 assembled-sensor qualification remains open, and it is deployment qualification rather than conditional mathematics.

## Mandatory next proof work

The next conversation should continue only along this path:

1. Replace the independent 300 m`*`s integral-displacement entry ball with a correlated one derived from the deployed S=0 regulation. With it as part of `L` the H18 prefix excursion reaches Cayley norm 4.5788 against the declared chart bound 1.0, so no chart-valid `L_chart` satisfies the theorem's own `Gamma*L+C_p < L_chart` and the hypothesis is unsatisfiable at that entry set. Without it the excursion is .9411 / .6537 and a chart-valid `L_chart` exists. Kinematics alone will not supply it: `S_next=S+dt*p+dt^2/2*v` is unleaked, so `|e_S|` reaches `20*T_handoff`, above 600 m`*`s at the declared live-entry timing floor.
2. Declare the chart-valid level `L_chart` numerically. `thm:brmm-bounded-bias-motion` already separates the entry level `L` from `L_chart` and asks for `Gamma*L+C_p < L_chart`, so an excursion above an entry radius is not a failure; leaving the chart is. Velocity needs roughly 4.0 (H18) and 4.9 (A21) times its entry radius, which is physically forced by gravity mis-resolution from the attitude ball and carries no chart constraint, so absorbing it is a declaration rather than a lemma. Attitude is the row that does carry one, and it fits at .9411 / .6537 once the independent integral ball is gone.
3. Complete the same-signal BRMM estimator cover so the input history jointly generates frequency, tau, raw/effective sigma, `T_S`, and SpectralMSE `R_S` while preserving candidate/active EMA, staged commits, and scheduler semantics.
4. Emit correlated same-history source cells for every admitted BRMM continuation, bias family, and hard-entry radial segment.
5. Derive `K` only from the same reachable `P/H/R` cell and retain exact reset transport/projection coupling.
6. Form source-correlated endpoint and literal-every-prefix augmented matrices for all three bias families. Run these per family: BIAS0 carries the largest admitted supply and BIAS2 the non-relaxing truth, so a budget that survives BIAS1 says nothing about either.
7. Run outward LDLT with finite-precision ISS charges and prove the full retained-coordinate budgets.
8. Only after **all three bias families plus the adaptive BRMM source cover** close may the final gate be changed to report their closure bits true.

## The single remaining object

Every one of the six blockers the gate still reports reduces to one missing object: the **source-uniform COMPLETE BRMM cover**. The correction/reset domain, the endpoint and every-prefix augmented LDLT and the every-prefix hard-domain retention each consume a source-uniform cell family that does not exist yet, so none of them can close before it. The cover itself needs the estimator-owned transition operator materialized over every admitted BRMM continuation, every hard-entry radial segment and every correlated Joseph cell, with the coefficient image proved inside the theorem's target cell.

The captured word cannot supply it, and the contract says so explicitly: `point_trace_can_promote_source_uniform_cover` and `trajectory_replay_or_pinned_generator_may_establish_uniform_cover` are both false. One word is a diagnostic; the cover is a statement about a continuum.

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
