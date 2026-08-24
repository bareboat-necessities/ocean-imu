# OU-III implementation stability proof-completion plan

## Objective

Prove stability of the **actual deployed OU-III implementation**, including the initialization sequence before the MEKF enters Live mode. The target is not an existential local result of the form “there exists a sufficiently small neighborhood.” The target is a numerical, machine-checked certificate for the source-realizable hybrid implementation.

The proof target is:

> Under an explicit configured-runtime and marine operating envelope, the deployed startup reset/Mahony observer reaches a certified Live handoff set; the exact held-bias Live mode (H), the exact held-to-active transition, the exact active-bias Live mode (A), magnetic lock/refinement, tilt reset/re-lock/cooldown, and periodic covariance synchronization remain inside source-node Lyapunov funnels; the exact nonlinear source-word maps enter a certified invariant inner funnel in finite time and thereafter satisfy the practical-ISS theorem. The stochastic layer must separately satisfy the paper's localized mean-square/concentration theorem.

“Actual implementation” means `SeaStateFusionFilter_OU_III` together with `Kalman3D_Wave_OU_III`, using the shipping update order, gains, adaptation law, clamps, covariance operations, reset maps, and startup state machine.

## Existing evidence retained from `main`

The existing tuning, eight-sea replay, multi-seed validation, exact executed-map replay, RMS gates, and evidence/provenance machinery stay in place. They remain important for performance regression, falsification, and locating active constraints.

They are **not** reimplemented in this PR and are **not** substituted for source-complete theorem bounds. A replay counterexample can invalidate a proposed bound; replay success alone cannot promote a theorem claim.

## Proof scope and unavoidable qualifications

The final claim is conditional only where the physics makes an unconditional theorem impossible:

1. **Configured runtime timing.** The public API accepts arbitrary positive finite `dt`; until code enforces a bounded interval, the numerical theorem is scoped to the source-defined configured scheduler.
2. **Heading persistent excitation.** Full-heading convergence requires recurring accepted, non-collinear gravity/magnetic (or equivalent heading) information. Permanent rejection/collinearity is unobservable. Gravity-only operation is treated on the yaw quotient.
3. **Physical startup/motion/noise envelope.** Specific force, initial gyro bias, heading-reference error, and stochastic noise must be bounded by explicit numerical theorem inputs. These are deployment assumptions, never fitted from the eight trajectories.
4. **SO(3) topology.** The smooth observer cannot have a globally asymptotically stable attitude equilibrium. However, the implemented first-accelerometer reset discards prior attitude and maps every stored prior attitude into a physical hemisphere whenever the non-gravitational specific force is smaller than gravity. The proof is therefore global with respect to the discarded pre-startup estimate and regional/almost-global with respect to the physical measurement ambiguity described in the paper.

## Non-negotiable proof rules

- No proof-specific filter, retuning, Schmidt restriction, fixed gain, disabled cross-gain, or altered covariance.
- Keep the complete `S=0` Kalman gain, including `S -> attitude` coupling.
- Bind each proof primitive to source and invalidate the certificate when source changes.
- Use outward-rounded/validated arithmetic for every continuous-source bound used for promotion.
- Cover accepted/rejected measurement branches that satisfy the theorem hypotheses, not only favorable replay words.
- Keep H (18-state) and A (21-state) fixed-dimensional words separate; H→A is an explicit dimension-changing hybrid jump.
- Include quaternion injection and left-error covariance reset in every applicable source word.
- Initialization is part of the final theorem composition.
- The final result must contain numerical margins and finite capture counts/times. `exists epsilon > 0` is not a completion criterion.

## Execution plan

### P0 — Freeze the exact implemented hybrid system

Create/extend a source-derived manifest containing:

- H/A state dimensions and order;
- prediction and process-noise construction;
- `S=0`, accelerometer, and asynchronous magnetometer update order;
- Joseph covariance updates and complete Kalman gains;
- quaternion injection and left-error reset;
- `tau`, `sigma_aw`, `R_S`, pseudo-update cadence and source clamps;
- startup Mahony gains, gravity gate, aligned-branch test, hold times, rate veto and timeout;
- `goLive` state/covariance initialization and bias-held entry;
- H→A bias release;
- magnetic lock/refinement/re-gauge;
- tilt reset/re-lock/cooldown;
- periodic `a_w` covariance synchronization;
- configured runtime sampling contract.

**PASS:** the proof manifest is regenerated from the source tree and semantic parity tests pin every operation used in the theorem.

### P1 — Numerical startup/reset certificate before Live

Instantiate the paper's Mahony/reset theorem with the deployed gains `2 k_P = 0.2`, `2 k_I = 0.02` and validated arithmetic.

Prove numerically:

1. the first accepted accelerometer reset is independent of stored prior attitude and maps into the desired gravity hemisphere for the declared `|a_f| < g` envelope;
2. positive startup Lyapunov/chart margins for the declared initial gyro-bias and disturbance bounds;
3. the normal gravity-gate tilt bound from the implemented `0.075` sine gate plus the declared measured-gravity error;
4. the finite normal-gate comparison time `T_Q` when its robust disturbance floor permits it;
5. the implemented timeout path at 150 s, including the aligned-branch constraint;
6. the exact `goLive` map into the H-mode source coordinates, including seeded covariance, linear-block enable, current OU/R_S state, and held accelerometer bias;
7. full-heading and yaw-quotient startup branches separately.

**PASS:** every source-reachable accepted startup/handoff branch has a finite numerical initial Lyapunov level strictly inside a certified H-mode capture funnel. No assumption that the handoff is “small enough” is permitted without a computed inequality.

### P2 — Complete source-word language

Make the normal-Live source language finite-window complete relative to explicit theorem hypotheses.

- Supply an explicit finite vector-PE recurrence window as a deployment theorem input.
- Preserve arbitrary accepted/rejected branches between required PE packets.
- Combine it with the rigorous pseudo-measurement firing-gap bound.
- Tile every fixed-mode normal-Live execution by bounded H/A words.

**PASS:** every normal-Live execution satisfying the declared PE/timing hypotheses belongs to the machine-defined language; no favorable-subset selection is possible.

### P3 — Validated H/A Riccati/information-word certificate — COMPLETE

P3 is now closed by two machine-checked layers.

First, the outward-rounded direct matrix backend certifies the source-uniform H/A generalized endpoint inequality with the full coupled `[v,p,S,a_w]` translation block. Four `S=0` firings provide the complete translation observability qualification, with the validated integer spread search choosing the strongest rigorous admissible four-firing information bound. The three-firing `[v,p,S]` construction plus stable `a_w` is used only to sharpen a covariance upper bound; it cannot by itself satisfy the P3 observability qualification. The H and A endpoint margins are both strictly positive, and the arbitrary former `1e-18` threshold is not a theorem condition.

Second, `tools/ou3_p3_word_algebra.py` binds the endpoint comparison to every fixed-dimensional normal-Live covariance operation in the shipping source. For each prefix it proves the exact decomposition

`P_s = Phi_s P_0 Phi_s^T + Omega_s`, with `Omega_s >= 0`.

Prediction, accepted Joseph corrections, rejected/not-due identity branches, the left-error covariance reset, and the periodic PSD `a_w` covariance increment all have the common affine-PSD form

`P+ = A P A^T + B`, with `B >= 0`.

Thus for every `0 < delta < 1`, an established information margin propagates exactly through any subsequent source branch:

`Omega+ - delta P+ = A (Omega - delta P) A^T + (1-delta) B >= 0`.

The implemented left-error reset has attitude block `G = I + 0.5 [dtheta]_x` and

`det(G) = 1 + ||dtheta||^2 / 4 >= 1`,

so reset congruence is nonsingular for every finite injected correction and does not require an additional small-angle hypothesis at P3. The same covariance decomposition gives, by Schur complement,

`Phi_s^T P_s^-1 Phi_s <= P_0^-1`,

hence the exact source-uniform worst prefix information gain is `1.0` rather than a replay estimate.

The final producer `tools/ou3_explicit_information_word_certificate.py` emits `P3_IMPLEMENTATION_WORD_CERTIFICATE=PASS` only when the direct H/A interval inequalities, optimized four-S qualification, current source-word language, exact Live-operation algebra, reset nonsingularity, branch coverage, and unit prefix-information bound all validate together.

**PASS:** strict source-uniform H and A information contraction with explicit positive endpoint margins and an exact finite prefix information-gain bound. P4 is the next unresolved normal-Live obligation.

### P4 — Exact nonlinear SO(3) word certificate

Lift P3 to the exact nonlinear implementation rather than a C1 existence argument.

Include:

- exact Rodrigues/quaternion attitude correction;
- actual vector measurement residual geometry;
- complete `S -> attitude` cross-gain;
- exact left-error reset;
- all continuous source parameters and accepted branch guards.

Use adaptive validated subdivision to prove, for an explicit `0 < theta_* < pi` and explicit source-node level,

- `mu_W > 0` for every H/A source word;
- every within-word prefix stays chart/domain safe;
- the endpoint remains inside the destination source-node funnel.

**PASS:** nonzero explicit nonlinear H/A funnel levels with positive `mu_W` and prefix-safety margins. This is the quantitative normal-Live certificate required by the paper; a local existence theorem is insufficient.

### P5 — Initialization-to-inner-funnel finite capture

Propagate the P1 startup handoff level through the certified H-mode words from P4.

Compute the source-node recurrence and a finite integer capture bound `N_H`, including the actual early covariance/tuner staging reachable immediately after `goLive`.

**PASS:** every certified normal or timeout handoff reaches the inner H funnel in finite words/time without leaving a certified prefix-safe domain.

### P6 — Prove every implemented hybrid jump

For each source-reachable hard event, validate the exact jump image and a strict destination inequality:

1. startup handoff;
2. H→A accelerometer-bias release (18→21 dimensions, including new-coordinate energy);
3. magnetic lock;
4. magnetic refinement/re-gauge;
5. tilt reset;
6. tilt re-lock;
7. cooldown re-entry;
8. periodic `a_w` covariance synchronization (nonexpansive PSD/Loewner proof, with strict decrease supplied by surrounding words).

After any jump that lands in an outer funnel, recompute finite recapture to the inner funnel.

**PASS:** every source-reachable jump lands strictly inside a certified destination funnel and the hybrid execution cannot escape through an unproved branch.

### P7 — Gravity-only quotient route

Machine-certify the no-magnetometer route on the yaw quotient using the same startup, translation, nonlinear, prefix and hybrid machinery, but without making an absolute-yaw convergence claim.

**PASS:** regional practical stability modulo yaw for the source-complete gravity-only implementation.

### P8 — Stochastic localized theorem

Using the same certified source words/funnels from P4-P6, validate the paper's Gaussian/localized drift and finite-horizon concentration constants from primitive source/noise bounds.

**PASS:** positive stochastic drift margin, localization radius inside the deterministic funnel, and recomputed finite-horizon failure probability below the declared theorem budget.

### P9 — Independent final implementation-stability composition gate

The final gate must regenerate source contracts and recompute all composition arithmetic. It must not trust an upstream `PASS` field.

The final status may be `PASS_IMPLEMENTATION_STABLE` only if all of the following are true simultaneously:

- source/implementation parity P0;
- startup/reset/handoff P1;
- source-word completeness P2;
- H linear certificate P3;
- A linear certificate P3;
- H nonlinear/prefix certificate P4;
- A nonlinear/prefix certificate P4;
- finite initialization capture P5;
- all hybrid jumps P6;
- gravity-only quotient result P7 when that mode is claimed;
- stochastic theorem P8 when the stochastic claim is requested;
- existing `main` performance/regression/evidence gates remain green.

Anything less is reported by its first failed proof obligation, never as “stable by existence.”

## PR completion criterion

This PR is complete only when it either:

1. produces `PASS_IMPLEMENTATION_STABLE` with the numerical margins/capture bounds above; or
2. identifies a mathematically genuine obstruction in the implemented filter/theorem assumptions, with the first failing validated inequality and no weakening of the proof gate.

A tiny unspecified local basin, replay-only success, or an analytical existence theorem is explicitly not an acceptable completion result.