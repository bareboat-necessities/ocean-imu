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
- `sigma_aw <= 4 m/s^2`, as the saturating commit `sigma_aw = min(0.9 sigma_a, 4)` rather than an assumed sea bound;
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

### Directional vessel/IMU response and response-weighted moments

`tools/ou3_sea3_directional_response_moments.py` is replay free and non-promoting.

The response family is an envelope, not a hull. With the lever arm disabled,
angular motion produces no specific force at the IMU, so only translational
response enters and the admissible family is bounded componentwise by

`|G_j(omega,theta)| <= rho(omega) omega^2 u_j(theta)`,
`u = [1, |cos theta|, |sin theta|]`,
`rho(omega) = rho_max min{1, (omega_L/omega)^p}`.

The roll-off is finite-waterline Froude-Krylov length averaging. Declared:
`rho_max = 2`, `p = 1`, `L_ref = 4 m`, hence `omega_L = 3.925 rad/s`. A shorter
reference hull raises `omega_L` and therefore assumes less.

Established analytical facts:

- **Modulus majorant.** For `|G_j| <= m_j` and any complex `v`,
  `v* Phi_u v <= |v|^T M |v|` with `M = int m m^T E dtheta`, so
  `lambda_max(Phi_u) <= lambda_max(M)`. The enclosure never reduces to
  independent per-axis scalars.
- **Directional factorization.** `E = sum_r S_r(omega) D_r(theta)` makes every
  response-weighted moment a sum of scalar frequency integrals times constant
  directional Gram matrices `Q_r(beta_r,s_r)`. Those are strongly non-diagonal
  (normalized off-diagonals to 0.998) and direction selective: a head sea loads
  surge (0.81 against 0.19), a beam sea reverses it.
- **Exact leak inversion.** With `d mu = |G_eta|^2 S_up d omega`, the deployed
  `omega^2 = (sigma_v/sigma_eta)^2 - lambda^2` equals `E_mu[omega^2]` exactly,
  for any input spectrum. It is an identity, not a narrow-band approximation.
- **Deployed period is high-passed.** Equivalently
  `d mu = rho^2 (omega^2/(omega^2+lambda^2))^4 S_eta d omega`: the deployed
  `T_z` is the moment period of the true elevation spectrum after a
  fourth-order high-pass at `lambda` and the response weight `rho^2`.
- **Mixture convexity.** `omega_est^2` is a convex combination over partitions,
  so the extreme deployed periods of the whole three-partition class occur on
  one-partition seas. The fifteen-parameter period enclosure collapses to a
  `(T_p, gamma)` screen.
- **Moment finiteness ladder.** Proxy elevation/velocity and band-passed
  acceleration moments are finite for any `p >= 0`; the raw acceleration moment
  needs `p > 0`; an acceleration moment of order `n` needs `n < 2p`. The
  deployed path needs no order above zero, so `p = 1` suffices.

Declared sea domain (SEA3-D): `T_p in [2.5, 20] s`, `gamma in [1, 7]`,
`H_s <= 8.5 m` (matching the declared operating domain, not narrower),
spreading exponent `s in [1, 25]`, `beta` free, and significant steepness
`<= 0.10` imposed **per partition and on the combined sea**. The steepness
ceiling is physical admissibility: without the per-partition form the search
returns a 0.28-steepness chop riding on a long swell, which is not a sea.

Numerical screens over that domain (fixed-grid quadrature with analytical tail
pads, not validated interval integration, therefore non-promoting):

- `0.656419 <= T_z^src/T_p <= 0.886099`, wider on both sides than the surface
  screen, because the high-pass shortens long swell and the roll-off lengthens
  short chop;
- `0.924 <= T_z^src/T_z,eta <= 1.115`; the surface ratio may not be substituted;
- induced tuning channel `0.0644 .. 0.5051 Hz`, inside the committed
  `0.03 .. 1.2 Hz`;
- worst admissible band-passed `sigma_a = 4.660 m/s^2`; at each sea's own
  deployed period only `3.905`, at the fixed 0.2 Hz startup prior `4.056`.

## Findings that change the next step

### Sampling supplies no band limit (obligation, not a failure)

The deployed chain is sampled at 5 ms, so any claim about what it observes must
bound the folded tail. With a flat response the specific-force density behaves
as `A/omega`, every decade above Nyquist folds in the same `A ln 10`, and the
total diverges. This is not a formal nicety: at the worst admissible partition
three unmodelled decades already fold in `83 m^2/s^4`, i.e. `9.14 m/s^2`, more
than the whole sigma clamp interval. Under the declared roll-off the same fold
is `2.36e-4 m^2/s^4`.

Therefore a certified vessel/IMU high-frequency roll-off is a **mandatory SEA0
theorem assumption**. Sampling cannot be offered in its place.

### The sigma clamp rail is reachable — failure analysis

1. **Failed inequality.** The declared sea domain was expected to map into the
   interior of the sigma source coordinate. It does not: the worst admissible
   sea gives `sigma_a = 4.660 m/s^2` against the saturation point
   `sigma_a^max/c_sigma = 4/0.9 = 4.444 m/s^2`. The witness is a crossing sea,
   `1.83 m` at `T_p = 4.59 s` on `8.18 m` at `T_p = 10.9 s`, `H_s = 8.38 m`,
   each partition below its own breaking steepness.
2. **Failure class.** Not a theorem failure, not an enclosure failure, not an
   implementation defect. It is a **source-language structure failure**:
   shipping code commits `sigma_aw = min(0.9 sigma_a, 4 m/s^2)`, so the
   coordinate is contained by saturation and no sea can violate it, but on the
   rail `sigma_aw` is constant and the sea-to-source map is not injective.
3. **What it invalidates.** The working hypothesis that a physical sea domain
   plus a directional response prunes the P2 tuner language *uniformly*. On the
   saturated cell SEA3 carries no more information than the frozen contract
   already does, so no pruning is available there.
4. **What it does not invalidate.** The SEA3 theorem formulation, the response
   family, the modulus majorant, the directional factorization, the leak
   inversion identity, the mixture convexity reduction, or the period channel —
   which stays comfortably interior to the committed tuning channel and is
   where the pruning leverage actually is.
5. **Current limiting quantity.** `sigma_a^worst / (sigma_a^max/c_sigma) =
   1.049`. It scales linearly in both `rho_max` and the steepness ceiling: the
   rail becomes unreachable at `rho_max <= 1.907` or steepness `<= 0.0954`.
6. **Sharpest structural fact.** Neither frequency-source mode reaches the rail
   at rest: `3.905` at the settled deployed period and `4.056` at the fixed
   prior are both below `4.444`. The rail is reached only at intermediate band
   references such as `0.258 Hz`, i.e. on the **estimator-lag branch between**
   the two source modes the startup subcertificate established. The two SEA0
   subcertificates therefore meet exactly here.

## Current certificate status

- **SEA0:** partial. Spectral shape/moment bridge, estimator startup/front-end, and directional response/response-weighted moment subcertificates exist. The hard finite-window oscillator/IQC sea enclosure, the exact discrete estimator/tuner reachability, and the rate relation `R_lambda` remain open.
- **P1:** existing operating-domain assumptions retained; SEA3-specific binding must include the prior/estimator source split above.
- **P2:** current implementation-correlated interface remains valid; SEA3-to-P2 inclusion is not yet proved.
- **P3:** existing canonical implementation-source route remains diagnostic. No SEA3 H/A finite-window contraction margin is promoted yet.
- **P4:** no SEA3 nonlinear endpoint/prefix certificate is promoted yet.
- **P5:** no SEA3 finite-capture certificate is promoted yet.

No new SEA0 artifact promotes P2, P3, P4, or P5.

## Current blocker

The directional response is now declared and its response-weighted moments are
constructed, so the previous blocker is discharged at the continuous-time,
steady-state level. Two things now stand between SEA3 and a usable P2 pruning.

**The sigma clamp rail.** SEA3 prunes the tuner language only where the
sea-to-source map is injective. It is not injective on the saturated
`sigma_aw = 4 m/s^2` cell, which the declared sea domain reaches. Either the
rail is put out of reach by a certified sharper response envelope, or the
saturated cell is carried explicitly and unpruned in the SEA3 P2 language.
Choosing the first for convenience is forbidden: `rho_max` and the steepness
ceiling are physical declarations and moving them to obtain a result is
operating-domain tuning.

**Discrete propagation.** The channel weights above are continuous-time
steady-state transfer functions. The deployed path is the exact discrete
WavePeriodEstimator moment EMAs, its canonical log-period state and settling
gate, the `AdaptiveWaveBandPass` recursion with moving corners, the debiased
tuner EMAs, stage/commit logic, and the pseudo scheduler. Nothing may be
promoted across that gap by analogy.

A Gaussian directional spectrum alone is still not an infinite-time
deterministic pointwise bound. The deterministic theorem needs a hard
finite-window oscillator/IQC enclosure; a stochastic sea-realization statement
is a later corollary.

## Next proof increment

Start from current `main` after this response/moment increment is merged. The
research question is the **source-language inclusion**, and it is one PR.

1. Decide the sigma clamp rail: certify a translational-RAO bound that puts
   `rho_max` below 1.907, or declare the saturated cell as an unpruned SEA3 P2
   state. Record which, and why, before any further construction.
2. Replace the continuous-time channel weights by the exact discrete
   WavePeriodEstimator moments, log-period state, settling gate, and the
   `AdaptiveWaveBandPass` recursion with moving corners.
3. Declare and justify the sea rate relation `R_lambda`; the estimator's finite
   memory only bounds source motion once the sea's own rate is bounded.
4. Drive the exact shipping variance/tau/sigma/R_S adaptation, stage/commit
   logic, and pseudo scheduler from those discrete sources.
5. Mechanically prove `Lhat_SEA3 subset L_current_source`, reporting which P2
   cells SEA3 actually removes and which it cannot.
6. Only after that inclusion passes, run high-precision complete-window H/A
   feasibility diagnostics with minimal enclosure pessimism.
7. Report worst H/A ratio, limiting sea parameters/directions/phase state,
   maximizing error direction, operation-by-operation margin consumption, and
   distance to `rho = 1`.

Interpretation rule: clear `rho < 1` justifies rigorous SEA3 enclosure work; near-one requires theorem/metric review; `rho > 1` falsifies the proposed SEA3 contraction formulation rather than authorizing blind subdivision.

If step 5 removes no P2 cell of consequence, that falsifies the pruning premise
of the SEA3 route itself, and the architecture review is due then rather than
after another enclosure attempt.

## Retired / forbidden proof shortcuts

Do not reintroduce without a new mathematical reason:

- independent Cartesian `tau/sigma/R_S` extrema;
- arbitrary per-sample acceleration boxes as the intended marine theorem source;
- replay-selected sea/source words as a certificate domain;
- equating `T_p` with deployed `T_z` or averaging modal periods;
- unbanded JONSWAP acceleration variance;
- offering the 5 ms sampling rate in place of a certified response roll-off;
- a single-level (whole-sea only) steepness ceiling, which admits breaking partitions;
- substituting the surface `Tz_eta` or a settled-run diagnostic for a bound that must also cover the fixed prior and the estimator-lag branch;
- per-axis scalarization of the directional source spectrum;
- treating `TunerReady` as WavePeriodEstimator readiness;
- common-period/Floquet reasoning that assumes three modal peaks are commensurate;
- blind subdivision, scalar-norm tightening, gate tuning, operating-domain shrink, or deployed-filter tuning to obtain a pass;
- parallel fallback P3/P4 workflows or diagnostic producers after their result has been absorbed into the canonical path.
