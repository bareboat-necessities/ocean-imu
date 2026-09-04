# OU-III proof research state

This is the short current-state ledger required by the root `AGENTS.md`. Replace stale research history rather than accumulating it here.

## Current theorem target

Certify OU-III for **perturbations of physically admissible oscillatory sea states**, not for every motion admitted by independent instantaneous source bounds.

Use a fixed `M_max = 3` directional spectral family

`E(omega,theta) = sum_{r=1}^3 S_J(omega; H_r,T_p,r,gamma_r) D_r(theta; beta_r,s_r)`

with inactive slots represented by `H_r = 0`, PM included by `gamma_r = 1`, and modal energy coupled by `H_s^2 = sum H_r^2`. Wind sea, swell, crossing seas and three-system multimodal seas are therefore one fixed-dimensional class.

`T_p` is a physical spectral parameter; it is not the deployed tuner period. The shipping WavePeriodEstimator estimates a zero-crossing/moment period after its own front-end dynamics. For multimodal seas the moments combine before the ratio is taken. SEA0 must therefore certify a physical-spectrum -> response -> estimator relation rather than equating `T_p` with the tuner period or averaging modal periods.

The SEA3 language must be source generated and satisfy

`L_actual_sea subset Lhat_SEA3 subset L_current_source`.

It refines, rather than redefines, `OU3_P2_CORRELATED_STAGE_TRANSFER_V1`. Exact WavePeriodEstimator state, fixed-prior startup branch, variance/parameter EMAs, stage/commit history and pseudo-update scheduler remain in the source path. No replay fitting or filter tuning is allowed.

## Controlling inequality

Three active spectral partitions need not have a common period. The controlling P4 object is a finite **sea-window** ratio, not a one-sample or one-cycle ratio:

`rho_W = V_after(F_W(x)) / V_before(x)`.

The target post-capture inequality is

`V_{k+N_W} <= rho V_k + gamma_s D_s,k + gamma_n D_n,k`, with `rho < 1`.

The window `T_W = N_W h` must be chosen from a certified finite-window information condition such as

`O_{k,N_W} >= alpha_O I`, `alpha_O > 0`,

on every admissible SEA3 window. It must not be selected from replay observations. P5 remains a separate finite-capture obligation from the declared 45-degree entrance into the P4 invariant funnel.

## Current proof architecture

1. **SEA0 -- three-mode directional-sea admissibility.** Certify physical `(H_r,T_p,r,gamma_r,beta_r,s_r)` domains and rate bounds, total-energy coupling, directional vessel/IMU response, response-weighted moment/estimator relation, and a source-generated oscillator/shaping or hard finite-window IQC enclosure.
2. **P1 -- operating branch.** Keep the existing startup/Normal-Live/hybrid conditions and bind them to SEA0. Every valid Normal-Live IMU sample executes the accelerometer update; no accelerometer-rejection branch is admitted.
3. **P2 -- SEA3 reachable language.** Preserve phase/frequency continuity, direction/cross-axis correlation, multimodal energy coupling, the fixed-prior -> estimator takeover, tuner state, stage/commit history and scheduler phase. Machine-check both language inclusions above.
4. **P3 -- finite-window linear certificate.** On the same SEA3 histories, certify H=18 and A=21 finite-window information/covariance bounds and strict canonical contraction. The canonical usefulness gate remains exactly `1e-18`.
5. **P4 -- nonlinear lifted dissipation.** Reuse the exact complete-word/Cayley/source-operation machinery on SEA3-complete windows and establish a useful invariant inner funnel. No pointwise decrease is required inside an oscillation.
6. **P5 -- finite capture.** Prove finite entry from the 45-degree entrance into that funnel, then compose consecutive SEA3 windows.

The existing arbitrary implementation-source route remains a stronger diagnostic envelope. It is not deleted and a failure there is not silently relabeled as a SEA3 pass.

## Current deployed theorem envelope

The proof must source/parity-check the tightened shipping values:

- `0.03 <= f_tune <= 1.2 Hz`;
- `0.02 <= tau <= 12 s`;
- `sigma_aw <= 4 m/s^2`;
- `0.15 <= R_S <= 100`;
- `T_S <= 0.150 s` with the nominal 200 Hz sample lower guard;
- dynamic EMA horizon `<= 35 s`;
- no accelerometer rejection in the declared Normal-Live branch.

The pseudo-scheduler starvation defect is no longer current: `set_pseudo_update_period_s()` uses the progress-preserving retarget helper, and the retained scheduler certificate gives the 30-sample maximum recurrence implied by the 150 ms clamp.

## Retained SEA0 evidence

### Spectral-moment bridge

`tools/ou3_sea3_spectral_moment_bridge.py` is replay free and non-promoting. It currently establishes:

- PM (`gamma = 1`) surface-elevation `Tz/Tp = ((5/4)*pi)^(-1/4) = 0.710370680986...`;
- for uncorrelated normalized modal partitions, `1/Tz_mix^2 = sum_r w_r/Tz_r^2`, `w_r = H_r^2/sum_j H_j^2`;
- ideal unbanded JONSWAP acceleration variance is not a valid finite theorem quantity because `S_eta ~ omega^-5` and `omega^4 S_eta ~ omega^-1`.

Its current `1 <= gamma <= 7` continuum screen reports

`0.709119591169 <= Tz_eta/Tp <= 0.829083684991`.

That shape screen is still a padded numerical screen rather than an interval-integration promotion and cannot prune P2. It is retained to size/falsify the response-weighted construction.

### WavePeriodEstimator front-end and startup causality

`tools/ou3_sea3_wave_period_frontend.py` is the second retained SEA0 subcertificate. It keeps the validated single-frequency discrete front-end identity and now also source-certifies the startup frequency-source split.

For the committed 5 ms schedule and 0.03--1.2 Hz tuning channel, the validated steady sinusoidal front-end distortion remains below about 59 ppm in period. This is not the current limiter.

The more important source-language result is:

- before `WavePeriodEstimator::getFrequencyHz()` is finite, the tuner uses the fixed `TUNE_FREQ_PRIOR_HZ = 0.2 Hz`;
- `update_tuner(..., tuner_frequency_hz_())` executes before the current sample calls `wave_period_.update(...)`, so the tuner consumes the previous sample's period-estimator state;
- when a period value first becomes finite, it can affect the tuner no earlier than the following valid sample;
- takeover uses the first finite positive estimator frequency and **does not wait for `WavePeriodEstimator::isReady()`**;
- filter `TunerReady` is based on `SeaStateAutoTuner` readiness, not WavePeriodEstimator readiness, so Live handoff may occur while the fixed-prior branch is still active;
- with the source default 0.02 Hz leak corner, no WavePeriodEstimator moment is accepted before the mandatory `6/lambda` integrator settling interval, about 47.75 s. This is a lower bound on first possible finite period, not a claim that a valid variance ratio must exist at exactly that time.

Therefore SEA3/P2 must carry at least two tuning-frequency source modes -- fixed prior and estimator-driven -- plus the one-sample takeover edge. `TunerReady` must never be used as shorthand for a settled sea-period estimate.

## Repository cleanup retained on this branch

The proof tree has been flattened rather than expanded:

- `ou3-proof.yml` is the sole general OU-III proof workflow;
- the temporary SEA0 workflow, one-shot cap-sync workflow, duplicate new-clamp P4 workflow and old whole-word probe workflow are removed;
- the superseded four-max and ordered-witness P3 diagnostics and tests are removed;
- the exploratory `ou3_sea3_initial_estimates.py` helper, JSON artifact and dedicated test are removed after their useful conclusions were absorbed into retained theorem/source contracts;
- `test_ou3_proof_cleanup.py` forbids these retired routes from returning.

Distinct studies such as the lever-arm workflow remain separate because they are not alternate proof pipelines.

## Current limiter

The raw spectral frequency-shape bridge, discrete WavePeriodEstimator sinusoidal warping, scheduler recurrence, and startup prior/takeover causality are no longer the immediate proof blockers.

The next missing soundness link is the **directional vessel/IMU response over the deployed finite band**, followed by the finite multimodal moment/log-period dynamics after estimator takeover. Until that exists, the spectral screen cannot soundly restrict the 800-state implementation P2 language.

A Gaussian directional spectrum by itself is not an infinite-time deterministic pointwise bound. The deterministic theorem must use a hard finite-window oscillator/IQC enclosure; a stochastic sea-realization statement is a later corollary.

## Next PR / next falsifiable experiment

This work is intentionally the boundary of the current PR. After it is merged, start a fresh PR from updated `main` for the response-weighted SEA3 -> P2 step:

1. declare and reference-justify a provisional three-partition directional JONSWAP/PM parameter/rate box and a conservative directional vessel/IMU response family;
2. form response-weighted elevation/velocity/acceleration spectral moments over the deployed finite band while preserving cross-axis/directional coupling;
3. propagate those moments through the retained fixed-prior/estimator WavePeriodEstimator source split, finite moment EMAs and log-period state;
4. drive the exact shipping variance/tau/sigma/R_S EMAs, stage/commit logic and pseudo scheduler, and mechanically test `Lhat_SEA3 subset L_current_source`;
5. only if that inclusion holds, compute high-precision complete-window H/A information/contraction ratios on boundary cells with essentially no interval pessimism;
6. report worst H/A ratio, limiting sea parameters/directions/phase state, maximizing error direction, operation-by-operation margin consumption, and distance to `rho = 1`.

Interpretation follows the root research protocol: clear `rho < 1` justifies rigorous SEA3 enclosure work; near-one requires a metric/theorem review; `rho > 1` falsifies this SEA3 P4 formulation rather than authorizing blind subdivision.

## Retained facts

- Canonical P3 usefulness threshold is exactly `1e-18`.
- Physical P2 tuner partition remains 800 states; arbitrary Cartesian tau/sigma/R_S switching remains forbidden.
- Same-history source correlation and exact operation order remain mandatory.
- H=18 and A=21 are both required.
- No replay fitting, operating-domain shrink, gate tuning or deployed-filter tuning is allowed to make the proof pass.
- Lever arm remains disabled and the vibration guard remains on its dormant/transparent proof branch.
- P4 cannot promote before canonical P3; P5 cannot promote before strict canonical P4 contraction.
- The theorem is regional for the current fixed 45-degree entrance. A semiglobal corollary requires certificates for arbitrary prescribed compact entrance regions.

## DEAD_ENDS / SHELVED

Do not resume these without a new mathematical fact:

- independent Cartesian tau/sigma/R_S extrema;
- arbitrary per-sample acceleration boxes as the intended marine theorem source;
- replay-selected sea/source words as a certificate domain;
- equating `T_p` with the deployed `T_z` estimator state or averaging modal periods;
- using unbanded JONSWAP acceleration variance as a finite theorem quantity;
- treating `TunerReady` as proof that WavePeriodEstimator is ready;
- fixed/common-period reasoning that assumes three modal peaks are commensurate;
- blind subdivision, scalar-norm tightening, coefficient tuning, domain shrink or gate tuning;
- additional P4 micro-certificates before a SEA3 source language and high-precision complete-window diagnostic are sound.
