# Minimal (27-state) realization of the shared-modes QMEKF (experimental)

**Status: experimental, simulation only, not a default.**  The shipping
filter stays `SeaStateFusionFilter_OU_III`.  The 33-state shared-modes QMEKF
([kalman-carma-shared-modes.md](kalman-carma-shared-modes.md),
`KalmanCarmaSharedModes<T>`) is unchanged and bit-identical
(`carma_shared_modes-test` [2]); this is a research alternative to it.

* Motion block: `src/kalman_carma_shared_modes/CarmaSharedModesMinimalMotion.h`,
  `CarmaMinimalSharedModeMotion<T, R>`, six states per axis,
  `R in {Cascade (default), Quotient, QuotientHold}`.
* QMEKF core: `KalmanCarmaSharedModesMinimal.h`,
  `KalmanCarmaSharedModesMinimal<T = float, R = Cascade>`, the shared-modes
  core with this block (27 error states).
* Wrapper: `SeaStateFusionFilter_CarmaSharedModesMinimal.h`,
  `SeaStateFusion_CarmaSharedModesMinimal<R>`, the CARMA3 schedule.
* Tests and tools: `tests/kalman_carma_shared_modes/carma_shared_modes_minimal-test.cpp`
  (run by `make all`), the bench names `CSM27`, `CSM27Q`, `CSM27H`, the
  replay variants `W3D_CSM_VARIANT=CSM27...`, the host timing tool
  `kalman_carma_shared_modes_minimal-timing` (on demand), the ESP32-S3
  probe `tests/kalman_carma_shared_modes/arduino/csm_size_probe`.
* Study: `tools/carma_shared_modes_minimal_analysis.py` (mathematics),
  `tools/carma_shared_modes_minimal_study.py` (Monte Carlo, robustness,
  traces, dataset, timing, size, report) ->
  `reports/results/kalman_carma_shared_modes_minimal/`.
* Article: `doc/kalman_carma_shared_modes/kalman_carma_shared_modes-minimal.tex`.

**Findings** (details in section 6; every number from the committed study data):

* **Equivalence.**  Exactly two directions per axis are unobservable, for
  every excitation and attitude, so a 27-state filter suffices.  At fixed
  parameters, and through every amplitude, coupling, confidence and fallback
  change, the 27-state QMEKF *is* the 33-state one (estimates, innovations,
  covariance).  Frequency adaptation moves the unobservable subspace
  (0.45 deg per percent of kappa), so exact equivalence under adaptation is
  impossible.
* **The removed states carry no useful information.**  The quotient with a
  regression lift (`CSM27Q`) stays within 0.12 % of the 33-state filter on
  all 26 simulated seas x 30 seeds, and within 0.02 % on 24 of them.  NEES
  agrees to a median ratio of 0.9999.
* **Recommended realization: the cascade (`CSM27`, the default).**  With the
  coupling off it is the CARMA3 baseline A41 exactly, adaptation included.
  It is the best conditioned of the realizations (correlation condition 5.5
  against 147).  Against the 33-state filter its 3-D error is better on 20
  of 26 seas (+0.1 to +0.9 %) and heave better by at least 1 % on 16.  It is
  worse on two seas: heave RAO x0.5 (-2.2 %) and coloured accelerometer noise
  (-1.6 %).  Consistency, wrong-coupling safety and fallback are unchanged.
* **Dataset replay:** gated by its own bars in `run_tests.sh`.  Against
  TFG-CARMA3 it has 7-9 % lower 3-D error and 20-25 % lower heave; the heave
  advantage is the excitation mode (with Evidence excitation TFG-CARMA3 is
  within 0.04 percentage points of Hs).  TFG-CARMA3 is better on pitch (4 %) and on the 3-D
  biases (2 %).
* **ESP32-S3:** -42 % multiply-adds per sample (13 686 -> 7 992), -51 % per
  model commit, -3.2 kB (-21 %) filter object, -592 B largest stack frame,
  +10.6 kB (+5.4 %) code.
* **Dimension:** the minimal Markov order of the three-axis law is
  `4 rank(Sigma) + 2 rank(Im Sigma)` (16 motion states, 25 in all).  That
  realization depends on the excitation and is not a better operating
  point.
* **Unexpected:** the shared excitation mapping over-couples when the hull's
  heave-to-surge ratio departs from the prior (applied coherence 0.38 for a
  measured 0.87).  With the heave RAO halved, every shared-mode filter, 33 or
  27 states, is 35 % worse than A41.  This is outside this change; see
  section 6.

## 1. Why two directions per axis are unobservable

Per axis the 33-state block carries the complexified C41J mode
`zeta = B/(P+ D) xi` as eight real states `[Re(p, v, a, l), Im(p, v, a, l)]`,
and every accelerometer row reads `Re a`.  Each root `r` of the real
quadratic `D` (-7.93 and -1.25 rad/s at kappa = 1) is an eigenvalue of both
`A_c` and `conj(A_c)`, so the realified `A8` has it twice with a
two-dimensional eigenspace, which one scalar output cannot resolve.  Exactly:

    U = {0} x (span_C{w(r3), w(r4)} intersected with R^4),   w(r) = [B(r), r B(r), r^2 B(r), 1]

is the unobservable subspace: two real directions per axis, entirely in the
imaginary half, independent of the excitation `Sigma`, of the attitude and
of every QMEKF operation.  Its annihilator is the pair of rows whose
polynomial is `D(s)(alpha + beta s)`; they read `Im m`, `Im m'` of the
D-free mode `m = xi / P+`.  Two of the four imaginary directions stay
observable for every `Sigma` (a correction to the companion article's
observability paragraph, which calls `Im zeta` unobservable for real
`Sigma`).

## 2. Equivalence

* **Fixed parameters, any excitation sequence.**  `U` is invariant and
  annihilated by every measurement row (any rotation), the reset touches
  only attitude rows, and the fallback shrink
  `P_hz <- E + alpha (P_hz - E)` commutes with the projection.  The
  27-state QMEKF therefore gives the same attitude, bias, `p`, `v`, `a`,
  innovations and covariance as the 33-state filter, through amplitude
  commits, estimated coupling, confidence and fallback (test [3]:
  5-7e-12 relative in double precision over 400 s).
* **Three-dimensional law.**  Per axis `Re zeta = B/P * eta`,
  `eta = V [xi_r; xi_i]`, `V = [N_r, N_i] / D` with `P+ = N_r + i N_i`;
  `|D|^2 = |N_r|^2 + |N_i|^2` makes `V` lossless (co-inner).  Every
  realization of this `1 x 2` transfer function with the noise blocks
  `g2 [Re(Sigma_ab) I + Im(Sigma_ab) J2]` gives the same spectral density
  matrix: co-spectra, quadrature spectra, coherences and phases (test [2]:
  transfer functions 3e-14, stationary law and lag cross-covariance 9e-11).
* **Hermitian structure.**  `J` does not descend to the quotient, so the
  27-state block discretizes two intensities:
  `Q_ab = Re(Sigma_ab) Q_R + Im(Sigma_ab) Q_I`, `Q_R` symmetric from
  `g2 G G^T`, `Q_I` antisymmetric from `g2 G J2 G^T`.

## 3. Adaptation

The CARMA3 schedule commits a frequency scale `kappa` and amplitudes every
0.1 s, and the coupling changes `Sigma`.  A commit replaces `F` and `Q` and
keeps the state, so the adaptive model depends on the realization.

* **Excitation changes** (amplitudes, coupling, confidence, shrink) leave
  `U` and `F` alone: exact (above).
* **Frequency changes** are not.  Every realization scales exactly,
  `A(kappa) = kappa S A(1) S^-1` with a diagonal `S`, but
  `U(kappa) = S U(1)` moves: 0.45 deg per percent of `kappa`, 10 deg at
  +-25 %.  After a commit the 33-state filter's component in the old `U`,
  which it carries with a data-dependent posterior, leaks into the new
  observable coordinates.  No 27-state filter carries it, so none is
  exactly equivalent under frequency adaptation (motion-only experiment:
  `|p_33 - p_27Q| = 3e-5 m`, not rounding).
* **Quotient**: at a kappa commit the missing component is replaced by its
  stationary regression on the quotient (`L = P8 Pi^T (Pi P8 Pi^T)^-1`), and
  `x <- Pi(k_new) L(k_old) x` on every axis (`p, v, a` exactly continuous).
  With the scaling this is one double factorization at construction.
* **QuotientHold** holds the quotient coordinates (ablation).
* **Cascade** (default): the CARMA3 baseline's own `[p, v, a, l]` driven at
  `a` by `eta = xi_r + C_w w`, and the two colouring states `w` of `V` in
  input-normal coordinates (`A_w + A_w^T = -B_w B_w^T`, `C_w = -B_w e_1^T`,
  `B_w` upper triangular).  A commit holds `[p, v, a, l]` (as A41 does) and
  the normalized `w / sqrt(g2 sigma_a^2)`.  Because `V` is lossless with a
  constant Gramian, `eta` stays exactly white under any `kappa(t)`,
  `sigma(t)`: **with the coupling off the cascade is the CARMA3 baseline A41
  in law, adaptation included** (test [4]: 1.6e-13 relative; the 33-state
  filter and the quotient differ from A41 by up to 22 % under the same
  schedule).  This removes the "realization effect" documented for
  `CSM_OFF` in [kalman-carma-shared-modes.md](kalman-carma-shared-modes.md)
  section 2 by construction.  With the coupling on it is a third adaptive
  model, judged by the experiments below.

## 4. Is 27 the right dimension?

The minimal Markov order of the three-axis displacement law is
`4 rank(Sigma) + 2 rank(Im Sigma)` (residues of the spectral density: at the
roots of `D` the real part cancels and the residue is proportional to
`Im Sigma`, a real antisymmetric matrix of rank 0 or 2; Hankel singular
values confirm 16, 16, 12, 12, 8 for full Hermitian, quadrature, real,
diagonal and rank-one excitations).

* Real (in-phase) coupling or none: 12 motion states, i.e. the CARMA3
  baseline or Experiment A.
* Quadrature coupling (linear waves): 16 motion states, a 25-state filter,
  with the extra four in the vertical plane of propagation.
* 27 is the minimal realization that keeps the shared dynamics
  `F = I3 (x) F6` for every excitation.

A 25-state filter exists (C41J on every axis plus two colouring copies in
the plane of `range(Im Sigma)`), but its coordinates depend on `Sigma`: they
would move at every coupling, confidence and heading update, which loses
the exactness under excitation changes.  It is undefined at
`Im Sigma = 0`, which the confidence crosses at every warm-up and fallback,
and it loses the Kronecker structure.  All that for two states, about 14 %
of the 27-state covariance work.  Not implemented.

## 5. Implementation notes

* Paired discretization `discretize_lti_pair` / `StepDiscretizationPair`:
  `discretize_lti`'s scaled series, doublings and jitter series for
  `(Q_R, Q_I)`.
* Per-axis re-coordinatization hook `take_axis_recoordinatization` in the
  CARMA3 core (`Kalman3D_Wave_CARMA3::recoordinate_motion_`, applied as
  `S P S^T` right after a commit).
* Exact structural zeros.  Round-off residue of `Pi A8 L` (down to 1e-19)
  defeated the discretization's structural scaling and produced denormals,
  a 7x slower refresh on x86.  The cascade's `B_w` is exactly upper
  triangular and `C_w = -B_w e_1^T` exactly.
* Closed-form 2 x 2 algebra (Lyapunov, Cholesky, Givens) in the cascade's
  construction.  Eigen's generic decompositions added 44 kB of code on the
  ESP32-S3.
* Conditioning (stationary correlation per axis): 147 for the 33-state
  block, 128 for the quotient, 5.5 for the cascade (= A41: the colouring
  states are white and uncorrelated with `[p, v, a, l]`).

## 6. Experiments

Paired Monte Carlo (`kalman_carma_shared_modes-bench`, carma3_bench engine:
random-phase directional seas through the benchmark hull's RAOs, one IMU and
magnetometer stream per row shared by every filter).  26 seas in five groups
(fixed; coupling strength and phase: spreading, heading, heave RAO x0.5 /
x1.5; gradual frequency change; rapid period, direction, height and
onset/cessation changes; disturbances: turns, manoeuvring, coloured and
variance-cycling accelerometer noise, +-25 % sample jitter, magnetic
disturbances with calibration residuals, accelerometer-bias drift).  30 seeds
each, 1200 s at 200 Hz, scored from 300 s, paired bootstrap over seeds.
Filters: `A41@C`, `CSM` (33), `CSM27` (cascade), `CSM27Q` (quotient + lift),
`CSM27H` (quotient held), `CSM:off`.  Raw rows:
`reports/results/kalman_carma_shared_modes_minimal/raw_runs.csv`; tables
and figures: the article (`csm27-*.tex-part`, figures in `figures/`).

Paired 3-D RMS change against CSM (%, + = better; 95 % interval half-width
about 0.05-0.5):

| group | CSM27 heave | CSM27 3-D | CSM27Q 3-D | CSM27H 3-D | A41 3-D |
|---|---|---|---|---|---|
| fixed (9) | +0.63 | +0.32 | +0.002 | -0.25 | -4.85 |
| coupling (6) | +0.63 | -0.09 | +0.021 | -0.01 | +3.62 |
| gradual (1) | +0.89 | +0.16 | +0.001 | -0.18 | -5.45 |
| rapid (3) | +0.42 | +0.03 | +0.004 | -0.27 | -3.42 |
| disturbance (7) | +0.69 | -0.01 | -0.002 | -0.15 | -2.29 |

* **Removed directions.**  CSM27Q's largest |change| is 0.12 % (heave RAO
  x0.5); the transient envelopes around the frequency, direction and height
  changes show `|p_27Q - p_33|` three to four orders below the error.
  Holding the quotient instead of lifting (CSM27H) costs up to 1.25 %.
* **Cascade.**  It is better where the complexified realization's adaptation
  cost is largest (regular wave heave +3.3 %, short wind sea +2.2 %, spread
  JONSWAP +1.5 %; CSM:off is a median 1.2 % worse than A41 in heave).  It is
  worse in horizontal RMS on JONSWAP-type seas (-0.1 to -0.45 %), in heave
  on swell, the abrupt change, bimodal and PM (-0.2 to -1.4 %), and in 3-D on
  heave RAO x0.5 and coloured accelerometer noise.  The commit-rule ablation
  (`CSM27/rawhold`: hold the colouring states rather than the normalized
  ones) does not remove these losses.  They come from the C41J coordinates
  themselves, which also produce the gains.
* **Consistency.**  Median position NEES 1.670 (CSM), 1.652 (CSM27), 1.670
  (CSM27Q); identical 95 % coverage (0.914 / 0.945 / 0.912); no failure in
  4680 filter runs.
* **Wrong coupling** (external prior reversed, turned, phase-shifted,
  doubled; 3 seas x 20 seeds).  With the confidence on, every wrong prior is
  withdrawn and both filters sit at A41's error.  Without it, the doubled
  prior is the worst case for both (57 mm, NEES 14).  The 27-state filter is
  within 1.6 % of the 33-state one in every cell.
* **Dataset replay and gates.**  CSM27 replays the eight scored records with
  the seeds and adapter of every CARMA filter and is gated by its own bars
  (worst default-seed value plus about 0.5 %;
  `kCsm27RegressionBars` in `kalman_carma_shared_modes-sim.cpp`):

| gate (worst of 8) | CSM27 | bar | TFG-CARMA3 (T41, Legacy) | T41 bar | T41, Evidence | CSM (33) |
|---|---|---|---|---|---|---|
| heave % Hs, JONSWAP | 2.455 | 2.467 | 3.057 | 3.072 | 2.459 | 2.518 |
| heave % Hs, PM-Stokes | 2.403 | 2.416 | 3.215 | 3.231 | 2.44 | 2.458 |
| 3-D %, JONSWAP | 4.41 | 4.433 | 4.865 | 4.89 | 4.907 | 4.397 |
| 3-D %, PM-Stokes | 4.879 | 4.904 | 5.252 | 5.278 | 5.277 | 4.849 |
| roll (deg) | 0.3571 | 0.359 | 0.3551 | 0.3569 | 0.3555 | 0.3572 |
| pitch (deg) | 0.222 | 0.2231 | 0.2128 | 0.2138 | 0.2131 | 0.2234 |
| yaw (deg) | 0.7556 | 0.7595 | 0.7623 | 0.7662 | 0.7631 | 0.7513 |
| accel z bias % | 4.164 | 4.186 | 4.185 | 4.206 | 4.173 | 4.166 |
| accel 3-D bias % | 83.56 | 83.99 | 82.12 | 82.53 | 82.2 | 83.75 |
| gyro 3-D bias % | 12.64 | 12.71 | 12.33 | 12.4 | 12.34 | 12.65 |

* **Over-coupling of the excitation mapping (unexpected; not changed
  here).**  On the heave RAO x0.5 sea, A41 and CSM:off are 35.4 %
  [32.6, 37.6] better in 3-D than CSM, and CSM's heave error exceeds A41's by
  more than 20 % on 28 of 30 seeds.  The measured multiple coherence is the
  same as on the unscaled sea (0.87), but the applied coherence rises from
  0.13 to 0.38.  The mapping keeps the measured regression `S_hz / S_zz`,
  the vertical prior and the horizontal prior, so the applied coherence is
  `gamma^2 (S_hh / S_zz) / (sigma_h^2 / sigma_z^2)`.  The consistency test
  compares two data-driven regressions and cannot see this.  It affects the
  33- and 27-state filters alike.  A candidate fix is to cap the applied
  coherence at the measured one; that would change the 33-state filter and
  needs its own validation.

## 7. Embedded efficiency

| filter | states | P (B) | multiply-adds / sample | per 0.1 s commit | ESP32-S3 filter object (B) | code (B) | largest frame (B) | update chain (B) |
|---|---|---|---|---|---|---|---|---|
| A41 | 21 | 1764 | 4,248 | 4,608 | 6136 | 173,375 | 1328 | 4832 |
| CSM | 33 | 4356 | 13,686 | 28,672 | 15776 | 196,266 | 2480 | 4896 |
| CSM27 | 27 | 2916 | 7,992 | 13,932 | 12528 | 206,920 | 1888 | 4864 |
| CSM27Q | 27 | 2916 | 7,992 | 19,440 | 13600 | 235,752 | 4864 | 4864 |

* Multiply-adds are the analytic counts of the CARMA3 core's structured
  implementation (`carma3_analysis.cost_model`, extended: Kronecker motion
  block, paired `Q_R`/`Q_I` additions, commit discretization and
  re-coordinatization).
* ESP32-S3 sizes: Espressif GCC `esp-14.2.0_20241119` (archive SHA-256
  matches Espressif's release checksum), `-Os`, arduino-esp32 flags, object
  of the complete wrapper (`csm_size_probe/build.sh`).  The filter object
  holds `P` (729 floats against 1089) and the motion discretization (paired
  6-state series against one 8-state series).  The wrapper itself holds its
  core through a pointer, hence equal `sizeof_top`.  The extra code is the
  paired discretization and the closed-form construction; generic Eigen
  decompositions there had cost 44 kB.
* Per commit the 33-state block also solves, in double precision (software
  on the ESP32-S3), a 4 x 4 Lyapunov equation by doubling and a complex 4 x 4
  inverse.  The cascade needs neither: its `kappa` scaling is closed form.
HOST_TIMING

## 8. Recommendation

* Use the 27-state filter wherever the 33-state shared-modes filter would be
  used: it is equivalent at fixed parameters and preserves the
  three-dimensional law, cross-spectra, Hermitian structure, `p, v, a`
  estimates and covariance consistency.  The removed states carry no
  information during adaptation, and it costs 42 % fewer multiply-adds per
  sample and 21 % less filter memory, at 5 % more code.
* Prefer the **cascade** (default).  Its coupling-off model is the CARMA3
  baseline exactly, it is the best conditioned, its commits are cheapest,
  and it is better on 20 of 26 seas.  Use the **quotient** when bit-level
  continuity with the validated 33-state filter matters more than those
  properties: it reproduces it to 0.12 % at the same per-sample cost, with
  more commit work and 29 kB more code than the cascade.
* Neither 25 states (excitation-dependent coordinates) nor 21 (no
  quadrature coupling) is a better operating point.
* The over-coupling of the excitation mapping, which both filters share, is
  the larger effect and the next experiment: cap the applied coherence at
  the measured one and re-run the heave RAO x0.5 and x1.5 seas together with
  the 26-sea set.
