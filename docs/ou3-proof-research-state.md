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

The existing arbitrary implementation-source P2/P3 route remains a stronger implementation envelope and is not relabeled as a separately constructed SEA3 certificate.

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

## Completed SEA0 / source-bridge subcertificates

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

Therefore P1/P2 retain two frequency-source modes — fixed prior and estimator-driven — plus the one-sample takeover edge. `TunerReady` is not interpreted as a settled sea-period estimate.

### Uniform directional RAO-family enclosure

`tools/ou3_sea3_directional_response_domain.json` and
`tools/ou3_sea3_directional_p2_ha_feasibility.py` quantify the CoG translational response over a continuum rather than selecting one hull.

Let

`G = ess sup_(omega,theta) ||h(omega,theta)||_2`.

The response theorem holds **for every finite `G >= 0`**, with arbitrary complex phase, arbitrary frequency dependence inside the deployed finite response band, arbitrary heading dependence, and arbitrary cross-axis coupling. A six-DOF parent RAO is allowed; `h` is its translational CoG projection on the theorem's zero-lever-arm branch.

No nominal RAO, finite RAO grid, or fixed numeric RAO gain cap is part of the theorem. The matrix outer product `h h*` is retained before the PSD/trace enclosure. Uniformly,

`tr M_disp <= G^2 H_s^2 / 16`,

`tr M_vel <= omega_hi^2 G^2 H_s^2 / 16`,

`tr M_acc <= omega_hi^4 G^2 H_s^2 / 16`.

These are symbolic `G^2` inequalities, not sampled gain cases.

### SEA3-to-current-P2 right inclusion

The same bridge mechanically proves

`Lhat_SEA3 subset L_current_source`

for every Normal-Live SEA3 realization generated by any finite RAO in the quantified family. The inclusion is independent of `G`: shipping clamps the period-derived frequency, `tau`, `sigma_aw`, and `R_S` before the exact EMA/stage/commit language already over-approximated by the frozen P2 source graph.

This is intentionally a **non-pruning inclusion**. It proves sound containment in the existing 800-state P2 partition; it does not yet declare any of those P2 cells unreachable.

If the unique canonical full-P2 H/A P3 certificate passes, that stronger full-source result is inherited by the entire SEA3 RAO family by set inclusion. If full-P2 P3 fails, SEA3 remains inconclusive until source-faithful SEA3 pruning is constructed. The canonical usefulness threshold remains exactly `1e-18`.

## Current certificate status

- **SEA0:** partial. Surface spectral bridge, estimator startup/front-end subcertificate, and a uniform finite-band complex RAO-family matrix-moment theorem now exist. The hard finite-window physical sea realization/IQC enclosure and the left inclusion `L_actual_sea subset Lhat_SEA3` remain open.
- **P1:** existing operating-domain assumptions retained. A physical sea/RAO pair must satisfy the declared Normal-Live hard response/source branch; no Gaussian pointwise bound is inferred from the spectrum alone.
- **P2:** the right inclusion `Lhat_SEA3 subset L_current_source` is now mechanically closed for the full quantified finite RAO family. The result is non-pruning.
- **P3:** H=18/A=21 feasibility is evaluated through the unchanged unique canonical full-P2 route. A full-P2 PASS transfers uniformly to the SEA3 family; a full-P2 FAIL is only inconclusive for SEA3.
- **P4:** no SEA3 nonlinear endpoint/prefix certificate is promoted yet.
- **P5:** no SEA3 finite-capture certificate is promoted yet.

No SEA3 bridge changes the filter, shrinks the declared operating domain, changes the P3 gate, or creates new P3/P4 promotion authority.

## Current blocker

The immediate remaining soundness link is now the **left physical inclusion** rather than a missing single-hull RAO.

SEA0 must construct a hard finite-window oscillator/IQC or equivalent deterministic realization enclosure for the three-partition directional sea and prove that the admitted physical sea/RAO pairs satisfy the existing P1 Normal-Live hard response/source bounds. A Gaussian directional spectrum alone is not an infinite-time deterministic pointwise bound; a stochastic sea-realization statement remains a later corollary.

If canonical full-P2 H/A is already green, RAO-specific P2 pruning is not required to establish P3 for the SEA3 subset. If it is not green, then the response-weighted finite WavePeriodEstimator/log-period/variance/tuner state must be propagated far enough to prune source histories before rerunning the same frozen H/A theorem interface.

## Next proof increment

1. Complete a replay-free hard finite-window oscillator/IQC realization enclosure for the three-partition directional sea family.
2. Bind the physical sea parameter/rate domain and the quantified RAO family to the existing P1 Normal-Live hard response/source limits, closing `L_actual_sea subset Lhat_SEA3`.
3. Read the unchanged canonical full-P2 H/A result. If it passes, record the inherited uniform SEA3 P3 margin; if it does not, propagate the exact finite-memory WavePeriodEstimator/log-period/variance/tuner state to obtain source-faithful SEA3 pruning.
4. Only if pruning is required, rerun the unchanged H/A theorem interface on the narrower SEA3 history language and report limiting sea parameters, RAO direction/phase state, maximizing error direction, operation-by-operation margin consumption, and distance to the canonical gate.
5. Continue to SEA3 P4 endpoint/prefix dissipation and P5 finite capture only after the P3 status is unambiguous.

## Retired / forbidden proof shortcuts

Do not reintroduce without a new mathematical reason:

- one nominal hull RAO, one fixed numerical RAO gain, or a finite sampled RAO catalog as a universal response proof;
- independent Cartesian `tau/sigma/R_S` extrema;
- arbitrary per-sample acceleration boxes as the intended marine theorem source;
- replay-selected sea/source words as a certificate domain;
- equating `T_p` with deployed `T_z` or averaging modal periods;
- unbanded JONSWAP acceleration variance;
- treating `TunerReady` as WavePeriodEstimator readiness;
- common-period/Floquet reasoning that assumes three modal peaks are commensurate;
- blind subdivision, scalar-norm tightening, gate tuning, operating-domain shrink, or deployed-filter tuning to obtain a pass;
- parallel fallback P3/P4 workflows or alternate diagnostic producers after their result has been absorbed into the canonical path.
