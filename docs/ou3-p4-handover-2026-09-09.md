# OU-III P4 handover — 2026-09-09

This file is the canonical continuation note for the next proof conversation after PR #509 is merged to `main`.

## Status at handover

P4 is **not closed**.

The branch deliberately keeps:

- `P4_MOTION_PASS = false`
- `P4_PASS = false`
- `P5_MAY_START = false`

No deployed filter code, physical domain, quality gate, or canonical P3 threshold `delta = 1e-18` was weakened for proof convenience.

PR #509 establishes substantial infrastructure and several real lemmas, but the source-uniform nonlinear theorem still lacks a complete same-history source cover and therefore lacks the final endpoint/every-prefix augmented LDLT and retained-domain certificate.

## What is genuinely established

### 1. Joint held-bias/source recurrence

The P4 architecture no longer treats accelerometer bias as an independent per-sample disturbance. It carries the physical true-bias root and driver together with corrected bias error and motion state.

For the driven conditional BIAS1 source,

`beta_i = phi_true * beta_(i-1) + w_i`

and the **same** driver increment enters corrected bias error:

`e_b,i = phi_hat * e_b,(i-1) + (phi_true - phi_hat) * beta_(i-1) + w_i`.

The 24-D joint recurrence, active projection, same-history finite factors, and compatible storage have passed on the nonzero-driver attached execution. Canonical evidence:

`reports/results/rao_stability/joint-nonzero-bias-driver.json`

The point run is evidence for algebra/connectivity only; it is not source-uniform P4.

### 2. Exact nonlinear radial projection sector

For radial projection `Pi_R` and

`F_R(e,beta) = beta - Pi_R(beta-e)`, the branch proves globally

`||Delta F_R||^2 <= ||Delta e||^2 + ||Delta beta||^2`

including saturated, unsaturated, and boundary-crossing branches. This avoids a frozen projection multiplier or saturation-pattern enumeration.

Canonical producer/result:

- `tools/stability/ou3_projection_sector.py`
- `reports/results/rao_stability/projection-sector.json`

Finite-precision projection is handled separately as an additive arithmetic channel.

### 3. Exact/reduced Joseph-reset algebra and prefix correction-history domination

The branch contains:

- reduced same-cell Joseph/reset signed identities;
- exact covariance-frame reset relations;
- exact reset graph/IQC helpers;
- every-prefix Joseph covariance-history domination;
- source-uniform finite-angle H18/A21 strict information blocks;
- exact chord/vector finite-angle machinery;
- explicit binary32 arithmetic ISS infrastructure.

These pieces are designed to preserve same-event/cross-term correlation and must not be replaced by detached rowwise gain bounds.

### 4. Hard entry set, BIAS1 conditional admission, and P3 premise separation

The proof uses the full declared hard entry set rather than shrinking an ellipsoid until a bound passes.

BIAS1 conditional root/driver admission and the canonical P3 execution premises are represented explicitly. Deployment/source qualification remains separate from the conditional mathematical theorem.

### 5. Source-cover infrastructure

The branch now has theorem-facing same-history source-cover APIs and interval transition operators for:

- shipping covariance prediction;
- Joseph covariance update;
- immediate covariance reset;
- tuner candidate/active EMA/commit state;
- pseudo-update scheduler branch preservation;
- raw/effective sigma distinction;
- WavePeriodEstimator high-pass/leaky-integrator/EW-moment state recurrence;
- correlated downstream target mapping once the upstream estimator coordinates are supplied.

Important files include:

- `tools/stability/ou3_p4_complete_brmm_source_cover_contract.py`
- `tools/stability/ou3_p4_complete_brmm_source_cover_transition.py`
- `tools/stability/ou3_p4_complete_brmm_adaptive_transition.py`
- `tools/stability/ou3_p4_complete_brmm_target_cell.py`
- `tools/stability/ou3_p4_wave_period_interval_transition.py`
- `tools/stability/ou3_p4_complete_brmm_frontend_geometry_transition.py`
- `tools/stability/ou3_p4_complete_brmm_differential_events.py`
- `tools/stability/ou3_p4_complete_brmm_finite_map_mean_value.py`

## Critical correction for the next conversation: f and sigma are NOT independent

This is the most important handover point.

Do **not** continue P4 using a rectangle in `(f, sigma)` or by treating period and sigma as independent tuner coordinates.

The deployed tuning architecture uses one physical signal history to generate both quantities:

1. `WavePeriodEstimator` filters the vertical-acceleration signal, maintains correlated EW velocity/elevation moments, and derives

   `T_z = 2*pi*sqrt(m0/m2)`

   from the moment ratio; its canonical state then gives exactly

   `f = 1/T_z`.

2. The **same** estimated wave frequency scales `AdaptiveWaveBandPass`, so the sigma-band filter coefficients depend on that same `f`.

3. The sigma statistic is estimated from the output variance of that period-scaled band after the propagated white-noise contribution is removed.

4. Only then are the deployed downstream coefficients derived:

   `f -> tau`,

   `(same signal, same f) -> sigma`,

   `tau -> T_S`,

   `(tau, sigma, T_S) -> R_S` via the deployed SpectralMSE map.

The theorem-facing coefficient-source obligation is therefore a **joint estimator-history relation**. The broad clamp rectangle remains only a sanity/admissibility superset and is explicitly non-promoting.

The branch was changed immediately before handover so that:

- `coefficient_target_inclusion_closed = false` until this joint signal/estimator relation is materialized;
- CI must reject any future proof that promotes from independent `(f,sigma)` target boxes;
- independent `R_S` and pseudo-cadence selection remain forbidden.

This correction supersedes earlier branch text that said clamp-defined coefficient inclusion was closed.

## Exact next mathematical task

Build a **joint same-signal estimator cell** and its branch-preserving transition.

It must carry, from one common physical/IMU signal history:

### Wave-period side

- two high-pass stages;
- leaky velocity/elevation states;
- EW weights, means, and second moments;
- the correlation needed to evaluate positive velocity/elevation variances without dividing independent interval boxes;
- valid/invalid raw-moment-ratio branch;
- canonical log-period state;
- one-way usable-period latch and invalid-update hold semantics;
- exact reciprocal relation `f = exp(-logT)`.

### Sigma side

- the same source sample entering the period-scaled `AdaptiveWaveBandPass`;
- band filter state and time-varying coefficients driven by the **same estimated f**;
- its propagated white-noise covariance state `p00,p01,p11`;
- the tuner variance accumulation/statistic used for sigma;
- the measurement-noise subtraction/floor exactly as shipping implements it;
- raw sigma target and later effective OU sigma floor as distinct coordinates.

### Downstream same-cell images

From that joint estimator state only:

- `tau_target(f)`;
- pseudo cadence `T_S(tau)`;
- SpectralMSE `R_S(tau,sigma,T_S)`;
- candidate/active EMA and staged-commit recurrence;
- actual applied anisotropic `R_S` for each S=0 event.

No future contraction experiment should run on a coefficient family that breaks those relations.

## After the joint estimator relation is materialized

Then return to the source-uniform P4 path:

1. emit correlated same-history `SourceCoverCell` sequences over all admitted BRMM continuations and every hard-entry radial segment;
2. keep reachable shipping Riccati `P`, event geometry, `R`, tuner state, scheduler state, BIAS1 true bias, and projection tied to the same branch;
3. derive `K` only from the same `P/H/R` cell;
4. use the exact/reduced Joseph-reset identities and finite-angle chord sector already in the branch;
5. form the source-correlated endpoint and literal-every-prefix augmented matrices;
6. run outward LDLT with the full finite-precision ISS charge included;
7. prove explicit retained coordinate budgets on the full declared hard entry set;
8. only if all of those pass may `P4_MOTION_PASS` or `P4_PASS` change.

## Dead ends / routes not to repeat

Do not spend time on the following unless new structure materially changes them:

- isolated spectral radius of a captured word;
- one independently selected Lyapunov metric per word without consecutive compatibility;
- detached per-sample bias-energy ports;
- marginal `Pbar` + hard residual boxes for reset-domain closure;
- rowwise independent `K` boxes for contraction/reset;
- scalar nonlinear eta budgeting against the tiny global Riccati/P3 margin;
- independent `(f,sigma)` tuner rectangles as a theorem family;
- covariance ellipsoids used as if they were hard physical entry sets;
- replay/pinned RAO history used as source-uniform proof.

The branch contains explicit diagnostics demonstrating why several of these lose too much correlation.

## Useful quantitative landmarks

These are diagnostics/lemmas, not P4 closure:

- nonzero-driver joint point critical entry levels: H18 about `4.42138 sigma`, A21 about `111.9058 sigma`;
- compatible-storage sample ratios on that point: H18 about `0.9999986022`, A21 about `0.9997983430`;
- finite-angle chart chord factor at 45 deg: `cos^2(22.5 deg) ~= 0.853553`;
- detached marginal same-cell reset-radius route can explode to roughly `9.62e7` on S=0 and is formally fenced off as non-promoting;
- old scalar moving-Riccati injection margin is about `2.13e-35` and is diagnostic only, not the active P3 theorem.

## CI state at handover

The handover was prepared while a large batch of PR checks was still queued/pending after the final same-signal guardrail commits. Do not infer overall green solely from this document.

When starting the next conversation from `main`, first inspect the current workflow state for the merge commit and rerun/fix only genuine regressions. The proof flags must remain false regardless of CI until the mathematical source-cover/LDLT obligations above are satisfied.

## Suggested first prompt for the next conversation

> Continue OU-III P4 from current `main` using `docs/ou3-p4-handover-2026-09-09.md` as the canonical handoff. Do not treat `(f,sigma)` as independent. Materialize the same-signal WavePeriodEstimator + AdaptiveWaveBandPass + sigma-statistic transition, preserve the joint estimator relationship through `(f,tau,sigma,T_S,R_S)`, then continue the same-history source cover and endpoint/every-prefix augmented LDLT. Keep P4/P5 false unless the full certificate actually closes.
