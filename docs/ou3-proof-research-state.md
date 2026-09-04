# OU-III proof research state

This file is the current proof ledger. It records only the active theorem, completed certificate facts, open obligations, and the next falsifiable step. Historical routes belong in Git history, not here.

## Active theorem

The target is **uniform regional practical stability with finite capture for perturbations of physically admissible multimodal directional sea states**.

The sea class has fixed dimension with

`M_max = 3`

and directional spectrum

`E(omega,theta) = sum_{r=1}^3 S_J(omega; H_r,T_p,r,gamma_r) D_r(theta; beta_r,s_r)`.

Inactive partitions use `H_r = 0`; PM is the `gamma_r = 1` boundary; modal energies obey

`H_s^2 = sum_r H_r^2`.

The theorem covers one-, two-, and three-system seas, including wind sea + swell and crossing seas. It does **not** prove stability for every arbitrary bounded acceleration sequence.

The physical peak periods `T_p,r` are not identified with the deployed tuner period. The shipping `WavePeriodEstimator` produces a zero-crossing/moment period after its own leaky front end and finite-memory statistics. For multimodal seas, spectral moments combine before the period ratio is formed.

The source language to be certified is

`L_actual_sea subset Lhat_SEA3 subset L_current_source`,

where `L_current_source` is the existing implementation-correlated P2 language. SEA3 refines that language; it does not redefine the frozen P2 contract.

Because multiple active spectral partitions need not be commensurate, the stability object is a finite **sea window**, not a common wave cycle. The required nonlinear endpoint inequality is

`V_{k+N_W} <= rho V_k + gamma_s D_s,k + gamma_n D_n,k`, `rho < 1`,

and P4 must also certify a uniform within-window prefix bound

`V_{k+l} <= kappa_V V_k + kappa_s D_s,k + kappa_n D_n,k`, `0 <= l < N_W`.

The prefix inequality is required for the theorem's full-sample practical-ISS claim; chart/source containment alone is insufficient.

## Current proof architecture

1. **SEA0 — directional sea admissibility.** Physical three-partition parameter/rate domain, directional vessel/IMU response, finite-band moment enclosure, source-generated oscillator/IQC enclosure, and exact estimator/tuner source dynamics.
2. **P1 — operating branch.** Startup/Normal-Live/hybrid assumptions, including no rejected accelerometer branch in declared Normal Live.
3. **P2 — SEA3 reachable language.** Preserve sea phase/frequency continuity, direction/cross-axis correlation, modal-energy coupling, prior/period-estimator source mode, tuner memory, stage/commit history, and pseudo-scheduler phase. Prove both language inclusions.
4. **P3 — finite-window linear certificate.** H=18 and A=21 information/covariance bounds and strict canonical contraction on the same SEA3 histories. The canonical usefulness threshold remains exactly `1e-18`.
5. **P4 — nonlinear lifted dissipation.** Prove endpoint contraction, quantitative prefix gain, and an invariant inner funnel for exact shipping complete words.
6. **P5 — finite capture.** Prove entry from the declared 45-degree entrance into the P4 funnel, then compose consecutive SEA3 windows.

The existing arbitrary implementation-source P2/P3 route remains a stronger diagnostic envelope and is not relabeled as the SEA3 theorem.

## Shipping theorem envelope

Current source/parity values used by the proof are:

- `imu_dt = 0.005 s`;
- `0.03 <= f_tune <= 1.2 Hz`;
- `0.02 <= tau <= 12 s`;
- `sigma_aw <= 4 m/s^2`;
- `0.15 <= R_S <= 100`;
- `0.005 <= T_S <= 0.150 s`;
- dynamic EMA horizon `<= 35 s`;
- every valid Normal-Live IMU sample executes the accelerometer update;
- accelerometer rejection is outside the declared Normal-Live theorem branch;
- lever arm remains disabled;
- the accelerometer vibration guard is restricted to its dormant/transparent proof branch.

The pseudo-update scheduler uses the progress-preserving retarget helper. The retained scheduler certificate gives a maximum 30-sample recurrence at the 150 ms ceiling; the former starvation/fmod route is obsolete.

## Completed SEA0 subcertificates

### Spectral-moment bridge

`tools/ou3_sea3_spectral_moment_bridge.py` is replay free and non-promoting.

Established analytical facts:

- PM (`gamma = 1`) surface-elevation period ratio
  `Tz/Tp = ((5/4)*pi)^(-1/4) = 0.710370680986...`;
- for uncorrelated normalized partitions,
  `1/Tz_mix^2 = sum_r w_r/Tz_r^2`,
  with `w_r = H_r^2 / sum_j H_j^2`;
- ideal unbanded JONSWAP acceleration variance is not a finite theorem quantity because `S_eta ~ omega^-5` while acceleration weighting gives `omega^4 S_eta ~ omega^-1`.

The current padded numerical JONSWAP shape screen over `1 <= gamma <= 7` gives

`0.709119591169 <= Tz_eta/Tp <= 0.829083684991`.

That continuum result is a feasibility screen, not an interval-integration promotion and not a P2 pruning certificate.

### WavePeriodEstimator front end and startup causality

`tools/ou3_sea3_wave_period_frontend.py` is replay free and non-promoting.

Completed facts:

- validated steady single-frequency discrete front-end period distortion is below about 59 ppm over the committed 5 ms / 0.03--1.2 Hz channel;
- before `WavePeriodEstimator::getFrequencyHz()` is finite, the tuner uses fixed `TUNE_FREQ_PRIOR_HZ = 0.2 Hz`;
- `update_tuner(..., tuner_frequency_hz_())` executes before the current sample calls `wave_period_.update(...)`;
- a newly finite WavePeriodEstimator value can therefore affect the tuner no earlier than the following valid sample;
- takeover uses the first finite positive estimator frequency and does **not** wait for `WavePeriodEstimator::isReady()`;
- filter `TunerReady` is based on `SeaStateAutoTuner` readiness, not WavePeriodEstimator readiness;
- with the default `0.02 Hz` leak corner, no WavePeriodEstimator moment can be accepted before `6/lambda ~= 47.75 s`; this is a lower bound on first possible finite period, not an exact readiness time.

Therefore P1/P2 must retain two frequency-source modes — fixed prior and estimator-driven — plus the one-sample takeover edge. `TunerReady` must not be interpreted as a settled sea-period estimate.

## Current certificate status

- **SEA0:** partial. Spectral shape/moment bridge and estimator startup/front-end subcertificates exist. Directional response, finite-band multimodal moments, hard finite-window sea enclosure, and full finite-memory estimator/tuner reachability remain open.
- **P1:** existing operating-domain assumptions retained; SEA3-specific binding must include the prior/estimator source split above.
- **P2:** current implementation-correlated interface remains valid; SEA3-to-P2 inclusion is not yet proved.
- **P3:** existing canonical implementation-source route remains diagnostic. No SEA3 H/A finite-window contraction margin is promoted yet.
- **P4:** no SEA3 nonlinear endpoint/prefix certificate is promoted yet.
- **P5:** no SEA3 finite-capture certificate is promoted yet.

No new SEA0 artifact promotes P2, P3, P4, or P5.

## Current blocker

The immediate missing soundness link is the **directional vessel/IMU response over the deployed finite wave band**.

That response must preserve directional/cross-axis coupling and turn the three-partition spectrum into bounded response-weighted elevation/velocity/acceleration moments. Those moments must then be propagated through the already-certified fixed-prior/estimator source split, finite WavePeriodEstimator moments/log-period state, exact variance/tau/sigma/R_S adaptation, stage/commit logic, and pseudo scheduler.

Until that construction exists, the physical sea spectrum cannot soundly prune the existing 800-state P2 tuner language.

A Gaussian directional spectrum alone is not an infinite-time deterministic pointwise bound. The deterministic theorem needs a hard finite-window oscillator/IQC enclosure; a stochastic sea-realization statement is a later corollary.

## Next proof increment

Start from current `main` after this theorem/SEA0-foundation increment is merged.

1. Declare and justify a provisional three-partition directional JONSWAP/PM parameter and rate domain.
2. Define a conservative directional vessel/IMU response family over the deployed finite band.
3. Compute response-weighted matrix spectral moments while preserving cross-axis/directional coupling.
4. Propagate them through the retained fixed-prior/estimator source split and finite WavePeriodEstimator/log-period dynamics.
5. Drive the exact shipping variance/tau/sigma/R_S adaptation, stage/commit logic, and pseudo scheduler.
6. Mechanically prove `Lhat_SEA3 subset L_current_source`.
7. Only after that inclusion passes, run high-precision complete-window H/A feasibility diagnostics with minimal enclosure pessimism.
8. Report worst H/A ratio, limiting sea parameters/directions/phase state, maximizing error direction, operation-by-operation margin consumption, and distance to `rho = 1`.

Interpretation rule: clear `rho < 1` justifies rigorous SEA3 enclosure work; near-one requires theorem/metric review; `rho > 1` falsifies the proposed SEA3 contraction formulation rather than authorizing blind subdivision.

## Retired / forbidden proof shortcuts

Do not reintroduce without a new mathematical reason:

- independent Cartesian `tau/sigma/R_S` extrema;
- arbitrary per-sample acceleration boxes as the intended marine theorem source;
- replay-selected sea/source words as a certificate domain;
- equating `T_p` with deployed `T_z` or averaging modal periods;
- unbanded JONSWAP acceleration variance;
- treating `TunerReady` as WavePeriodEstimator readiness;
- common-period/Floquet reasoning that assumes three modal peaks are commensurate;
- blind subdivision, scalar-norm tightening, gate tuning, operating-domain shrink, or deployed-filter tuning to obtain a pass;
- parallel fallback P3/P4 workflows or diagnostic producers after their result has been absorbed into the canonical path.
