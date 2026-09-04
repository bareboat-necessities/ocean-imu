# OU-III proof research state

This is the short current-state ledger required by the root `AGENTS.md`. Replace stale research history rather than accumulating it here.

## Current hypothesis

The stability theorem should be certified for **perturbations of physically admissible oscillatory sea states**, not for every motion admitted by independent instantaneous source bounds.

Use a fixed `M_max = 3` directional spectral family

`E(omega,theta) = sum_{r=1}^3 S_J(omega; H_r,T_p,r,gamma_r) D_r(theta; beta_r,s_r)`

with inactive slots represented by `H_r = 0`, PM included by `gamma_r = 1`, and modal energy coupled by `H_s^2 = sum H_r^2`. Wind sea, swell, crossing seas and three-system multimodal seas are therefore one fixed-dimensional class.

`T_p` is a physical spectral parameter; it is not the deployed tuner period. The shipping WavePeriodEstimator's canonical state is the zero-crossing/moment period `T_z = 2*pi*sqrt(m0/m2)` after its response/front-end dynamics. For multimodal seas the moments combine before the ratio is taken. SEA0 must therefore certify a `T_p -> T_z` relation through the declared directional response and exact estimator rather than equating the two periods or averaging modal periods.

The SEA3 language must be source generated and satisfy

`L_actual_sea subset Lhat_SEA3 subset L_current_source`.

It refines, rather than redefines, `OU3_P2_CORRELATED_STAGE_TRANSFER_V1`; exact WavePeriodEstimator, variance/parameter EMAs, stage/commit history and pseudo-update scheduler remain in the source path. No replay fitting or filter tuning is allowed.

## Controlling inequality

Three active spectral partitions need not have a common period. The controlling P4 object is therefore a finite **sea-window** ratio, not a one-sample or one-cycle ratio:

`rho_W = V_after(F_W(x)) / V_before(x)`.

The target post-capture inequality is

`V_{k+N_W} <= rho V_k + gamma_s D_s,k + gamma_n D_n,k`, with `rho < 1`.

The window `T_W = N_W h` must be chosen from a certified finite-window information condition such as

`O_{k,N_W} >= alpha_O I`, `alpha_O > 0`,

on every admissible SEA3 window. It must not be chosen from replay observations. P5 remains a separate finite-capture obligation from the declared 45-degree entrance into the P4 invariant funnel.

## Current proof architecture

1. **SEA0 -- three-mode directional-sea admissibility.** Certify physical `(H_r,T_p,r,gamma_r,beta_r,s_r)` domains and rate bounds, total-energy coupling, directional vessel/IMU response, the induced `T_p -> T_z` moment/estimator relation, and a source-generated oscillator/shaping or hard finite-window IQC enclosure.
2. **P1 -- operating branch.** Keep the existing startup/Normal-Live/hybrid conditions and bind them to SEA0. Every valid Normal-Live IMU sample executes the accelerometer update; no accelerometer-rejection branch is admitted.
3. **P2 -- SEA3 reachable language.** Add enough sea memory to preserve phase/frequency continuity, direction/cross-axis correlation, multimodal energy coupling, tuner state, stage/commit history and scheduler phase. Machine-check both language inclusions above.
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

The earlier pseudo-scheduler starvation defect is no longer the current limiter: `set_pseudo_update_period_s()` now calls the progress-preserving retarget helper. The SEA3 proof must certify the scheduler as it exists now rather than carrying the obsolete `fmod` failure forward.

## Evidence

`tools/ou3_sea3_spectral_moment_bridge.py` now provides the first replay-free SEA0 subcertificate. It does **not** promote full SEA0 or any P2--P5 stage.

Analytical results already fixed by this subcertificate:

- PM (`gamma = 1`) surface-elevation period ratio is exactly `Tz/Tp = ((5/4)*pi)^(-1/4) = 0.710370680986...`;
- for uncorrelated normalized modal partitions, `1/Tz_mix^2 = sum_r w_r/Tz_r^2`, `w_r = H_r^2/sum_j H_j^2`; modal periods are therefore combined through energy-weighted inverse squares, not an arithmetic average;
- ideal unbanded JONSWAP acceleration variance is not finite because `S_eta ~ omega^-5` and `omega^4 S_eta ~ omega^-1`; the vessel/IMU response and deployed wave band are mandatory before acceleration statistics can enter the theorem.

The current non-promoting numerical shape screen uses `1 <= gamma <= 7`, the `0.07/0.09` JONSWAP peak widths, 240 gamma cells, and conservative quadrature/tail pads. It reports the surface-elevation screening enclosure

`0.709119591169 <= Tz_eta/Tp <= 0.829083684991`.

That continuum screen is **not** a rigorous interval-integration certificate and cannot prune P2. Its role is to size and falsify the next response-weighted construction. The machine contract explicitly forbids substituting surface `Tz_eta` for the shipping tuner `Tz`.

`tools/ou3_sea3_initial_estimates.py` separately reads the current source clamps. Its artifact remains `exploratory_non_promoting`, source generated and replay free. At `h = 5 ms` and committed tuning frequency `0.03..1.2 Hz`:

- an individual oscillator advances `0.054..2.16 deg` per IMU sample;
- over the maximum 150 ms S-update gap it advances `1.62..64.8 deg`;
- the maximum nominal S-update gap is 30 IMU samples;
- the committed tune-period envelope is `0.833333..33.333333 s`;
- using the WavePeriodEstimator four-period / 20--180 s moment rule only as a sizing proxy gives `20..133.333333 s`, or about `4,000..26,667` IMU samples;
- 2/4/6 oscillator pairs per each of three partitions would contribute 12/24/36 sea-shaping states before tuner/scheduler augmentation.

These are sizing numbers only. The tune-frequency envelope is not asserted to be a physical `T_p` range.

## Current limiter

The raw JONSWAP/PM **frequency-shape bridge is no longer the immediate limiter**. The next missing soundness link is the directional vessel/IMU response plus deployed finite band, followed by exact propagation through the WavePeriodEstimator's leaky filters, finite moment EMAs and log-period state.

Until that response-plus-estimator enclosure exists, the spectral `Tz_eta/Tp` screen cannot restrict the 800-state P2 tuner language. After it exists, the next machine obligation is to prove

`Lhat_SEA3 subset L_current_source`

with the exact tuner/scheduler history retained. Only then is a new SEA3 P3/P4 feasibility diagnostic meaningful.

A Gaussian spectrum by itself is not an infinite-time deterministic pointwise bound. The deterministic theorem must use a hard finite-window oscillator/IQC enclosure; a stochastic sea-realization statement is a later corollary.

## Alternatives

1. **Finite oscillator/shaping bank.** Build an outward frequency/direction factor for each partition and retain phase state explicitly; propagate the RAO output through the exact estimator and tuner.
2. **Hard dynamic IQC outer enclosure.** Replace a large oscillator bank by a finite-window hard IQC after the response-weighted spectral envelope is known.
3. **Parameter-dependent lifted metric.** If a common source metric remains too conservative after the SEA3 language is sound, use `M(zeta)` with certified sea/rate transitions rather than further scalar subdivision.

The first two are alternative SEA0 representations; they must enclose the same declared physical sea class.

## Next falsifiable experiment

Build a **response-weighted SEA0 feasibility diagnostic** before any new interval P3/P4 work:

1. declare and reference-justify a provisional three-partition directional JONSWAP/PM parameter/rate box and a conservative directional vessel/IMU response family;
2. form response-weighted elevation/velocity/acceleration spectral moments over the deployed finite band, preserving cross-axis/directional coupling;
3. propagate those sources through the exact WavePeriodEstimator front end, including the two high-pass stages, leaky integrations, moment EMAs, log-period state, startup/ready condition and clamps;
4. drive the exact shipping variance/tau/sigma/R_S EMAs, stage/commit logic and pseudo scheduler, and mechanically test `Lhat_SEA3 subset L_current_source`;
5. if the inclusion holds, compute high-precision complete-window H/A information/contraction ratios on representative boundary cells with essentially no interval pessimism;
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
- fixed/common-period reasoning that assumes three modal peaks are commensurate;
- blind subdivision, scalar-norm tightening, coefficient tuning, domain shrink or gate tuning;
- additional P4 micro-certificates before a SEA3 source language and high-precision complete-window diagnostic are sound.
