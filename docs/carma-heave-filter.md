# Experimental CARMA heave filter versus OU-III

**Status: experimental research comparison.**  Nothing here changes the
shipping OU-III estimator, its tuning, gates, magnetic handling or the OU-III
stability study.  `Kalman3D_Wave_OU_III` and `SeaStateFusionFilter_OU_III`
remain the deployed, authoritative implementation; the CARMA code lives in
`src/kalman_carma/` and is used only by tests and the benchmark.  The
CARMA prior is a statistical model, not a law imposed on the vessel's
physical history, and no result below is a stability or admissibility
statement (see the last section of
`doc/kalman_ou_iii/stochastic-marine-motion.tex`).

__RESULTS_SUMMARY__

## 1. Model

The vertical displacement `h` (NED `p_z`, about the local equilibrium) is a
Gaussian CARMA process with transfer function `H(s) = g B(s)/P(s)` (monic
`B`, `P`) and two-sided PSD

    S_h(w) = g^2 |B(iw)|^2 / |P(iw)|^2,      Var h = (1/2pi) int S_h(w) dw.

| Model | B(s) | P(s) | states |
|---|---|---|---|
| CARMA(3,0) | 1 | (s+lambda) L1(s) | 3 |
| CARMA(4,1) | s + w_z0 | L1(s) L2(s) | 4 |
| CARMA(5,2) | s^2 + 2 zeta_z w_z s + w_z^2 | (s+lambda) L1(s) L2(s) | 5 |

with `L_j(s) = s^2 + 2 zeta_j w_j s + w_j^2` and `w = 2 pi f`.

**Relative degree and acceleration.**  `deg P - deg B = 3` for all three,
so `S_h = O(w^-6)`, `int w^4 S_h dw < inf`, and `h`, `v = h'` and `a = h''`
exist in mean square with finite variance; `a` is continuous but not
differentiable (its derivative carries the white driving noise).  JONSWAP has
an `f^-5` displacement tail, so none of these models is asymptotically
equivalent to it.

**Parameters.**  `tools/carma_fit.py` minimizes the RMS log-spectral error
`sqrt(mean_f [d(f) - mean_f d]^2)`, `d = 10 log10(S_C/S_J)`, against JONSWAP
(`f_p` = 0.29 Hz, `gamma` = 3.3) on a 1 mHz grid over 0.17-1.30 Hz, i.e. shape
only (the driving intensity is set at run time).  Constraints: every pole
damping in [0.05, 1]; CARMA(5,2) `zeta_z >= 0.20`; CARMA(4,1) `f_z0 >= 0`.
Results are in `reports/results/carma_benchmark/carma_fit.json`.

| Model | parameters (Hz, damping) | RMS dB, 0.17-1.30 Hz | max dB | DC level re peak | `T_z f_p` |
|---|---|---|---|---|---|
| CARMA(5,2), preliminary (used) | f1 0.27812/0.09715, f2 0.44274/0.52299, f_lambda 1.86402, f_z 0.15105/0.20 | 0.84 | 9.7 | -23.2 dB | 0.783 |
| CARMA(5,2), refit, same constraints | f1 0.27678/0.09708, f2 0.44283/0.49116, f_lambda 2.24923, f_z 0.15229/0.20 | 0.83 | 9.4 | -23.2 dB | 0.777 |
| CARMA(5,2), unconstrained (`zeta_z >= 0.01`) | f1 0.28168/0.1479, f2 0.66801/0.86989, f_lambda 0.06859, f_z 0.16952/0.02709 | 0.34 | 1.1 | -2.2 dB | 0.895 |
| CARMA(4,1) | f1 0.306945/0.104191, f2 0.815318/0.478188, f_z0 0 | 1.86 | 17.5 | zero at DC | 0.782 |
| CARMA(3,0) | f1 0.386457/0.266692, f_lambda 5.0 (on bound) | 2.93 | 24.0 | -5.8 dB | 0.766 |
| JONSWAP | | | | | 0.779 |

The preliminary constrained CARMA(5,2) parameters are used unchanged: the
refit under the same constraints gains 0.006 dB.  The unconstrained fit
reproduces the warned-about behaviour: a narrow notch (`zeta_z` = 0.027) at
0.17 Hz on the band edge and a 0.069 Hz real pole that puts 23 % of the
displacement variance below 0.1 Hz.  It is reported, not used.  The
CARMA(4,1) optimum places its zero at DC; the CARMA(3,0) optimum sits on the
5 Hz real-pole bound (an oscillator-like `f^-4` slope across the band).

![PSD overlay](../reports/results/carma_benchmark/carma_psd_overlay.svg)

All curves use the same frequency axis and are each normalized to their own
peak, i.e. shapes only.  The fractional oscillator mentioned in the study
request is not in the repository and is not plotted.

## 2. State-space realization

`src/kalman_carma/CarmaVerticalModel.h` writes every model in the
coordinates

    x = [h, v, a, l_1 .. l_{N-3}],    dx = A x dt + g e_a dW,

so that `h = C x = x_0`, `v = C A x = x_1`, `a = C A^2 x = x_2` exactly (the
first two rows of `A` are `h' = v`, `v' = a`).  This is the observability
change of basis applied to the factored (cascade) realization; the
remaining coordinates are latent states of that cascade, not additional
physical variables.  For CARMA(5,2) the cascade is

    u' = -lambda u + g xi,     r'' + 2 zeta2 w2 r' + w2^2 r = u,
    h'' + 2 zeta1 w1 h' + w1^2 h = u + beta1 r' + beta0 r,
    beta1 = 2 zeta_z w_z - 2 zeta2 w2,  beta0 = w_z^2 - w2^2,

(`B/L2 = 1 + (beta1 s + beta0)/L2`), and eliminating `u` through
`a = u + beta1 r' + beta0 r - w1^2 h - 2 zeta1 w1 v` gives, with
`gamma = beta1 - lambda` and latent `(l1, l2) = (r, r')`,

    a'  = gamma w1^2 h + (2 gamma zeta1 w1 - w1^2) v + (gamma - 2 zeta1 w1) a
          + (-beta1 w2^2 - gamma beta0) l1 + (beta0 - 2 beta1 zeta2 w2 - gamma beta1) l2 + g xi
    l1' = l2
    l2' = w1^2 h + 2 zeta1 w1 v + a - w_z^2 l1 - 2 zeta_z w_z l2.

The latent pair therefore obeys the numerator ("zero") dynamics driven by the
input of the first oscillator; it is stable because `zeta_z > 0`.  At
`kappa` = 1 the CARMA(5,2) matrix is

    [  0        1        0        0        0     ]
    [  0        0        1        0        0     ]
    [-43.491   -7.889  -14.582  -77.805  -35.510 ]
    [  0        0        0        0        1     ]
    [  3.054    0.340    1       -0.901   -0.380 ]

with eigenvalues -11.712, -0.170 +- 1.739i, -1.455 +- 2.371i (the poles of
`P`).  CARMA(4,1) has one latent state `r` (the input of `L2`):
`r' = w1^2 h + 2 zeta1 w1 v + a - w_z0 r` and
`a' = d w1^2 h + (2 d zeta1 w1 - w1^2) v + (d - 2 zeta1 w1) a + (-w2^2 - d w_z0) r + g xi`,
`d = w_z0 - 2 zeta2 w2`.  CARMA(3,0) is the companion form in `[h, v, a]`:
`a' = -lambda w1^2 h - (2 lambda zeta1 w1 + w1^2) v - (lambda + 2 zeta1 w1) a + g xi`.

In these coordinates OU-III's vertical chain is the special case
`a' = -a/tau + noise` with no `h`, `v` feedback (and an extra integral `S`);
CARMA feeds displacement and velocity back into the jerk equation, which is
what makes its displacement stationary.

**Scaling.**  A fixed shape is moved by one frequency scale `kappa`
(every pole and zero times `kappa`).  With unit gain this gives exactly
`Var_kappa(a) = Var_1(a)/kappa` and `Var_kappa(h) = Var_1(h)/kappa^5`, so a
requested stationary acceleration standard deviation `sigma_a` is met with
`g^2 = sigma_a^2 kappa / Var_1(a)` without a run-time Lyapunov solve; the
displacement zero-crossing frequency is `kappa` times its base value
(0.3787, 0.3708 and 0.3704 Hz for CARMA(3,0), (4,1), (5,2)).

## 3. Exact discretization

For any step `h`,

    F = exp(A h),   Q = int_0^h exp(A t) W exp(A^T t) dt,   W = g^2 e_a e_a^T,

is computed by `discretize_lti()`: Taylor series of both integrands at
`h_s = h/2^s` with `||A h_s||_inf <= 1/4`, evaluated in the diagonally
scaled coordinates `x = D z`, `D = diag(h_s^{d_i})`, `d_i` the coupling
distance of state `i` from the noise input, followed by `s` exact doublings
`Q(2h) = Q(h) + F Q F^T`, `F(2h) = F^2`.  The scaling keeps the small entries
(`Q_hh ~ h^5`) at full relative precision and avoids subnormal intermediates;
every series term is a product, never the difference `P - F P F^T` of
stationary covariances, which cancels catastrophically in float at 5 ms.  `Q`
is a sum of positive semidefinite terms up to rounding.  The stationary
covariance is the same doubling carried until `F` has decayed.  Non-uniform
steps are supported directly (the filter caches `F`, `Q` for an unchanged
step and model; with a varying `dt` they are recomputed every sample).

`tests/kalman_carma/carma_model-test.cpp` (801 checks) verifies, for all three
orders and several `kappa`:

* the realization reproduces `g^2|B|^2/|P|^2` from the polynomial factors
  (relative error < 1e-9, 5 mHz-20 Hz), `C G = C A G = 0`, `C A^2 G = g`, and
  the eigenvalues equal the factor poles, all strictly stable;
* `F`, `Q` against a double-precision Van Loan matrix exponential (scaled
  coordinates) for steps 0.1 ms-50 ms, and against `P - F P F^T` for
  0.4-40 s: < 1e-9 relative to each entry's natural scale; the float
  implementation within 1e-5 (F, whitened) and 1e-4 (Q) at IMU-rate steps;
* `Q` positive semidefinite; the semigroup law for non-uniform steps
  `Q(h1+h2) = F(h2) Q(h1) F(h2)^T + Q(h2)`;
* the doubling stationary covariance against a Kronecker Lyapunov solve
  (< 1e-9) and against quadrature of `S_h`, `w^2 S_h`, `w^4 S_h` (< 1e-5);
* the discrete recursion `P <- F P F^T + Q` converging to the continuous
  stationary covariance (double < 1e-9, float < 2e-3);
* the `kappa` scaling law; a 10 000 s simulation staying bounded.

## 4. Integration into the INS

`Kalman3D_Wave_CARMA<float, N, with_vertical_integral>`:

| block | states | from |
|---|---|---|
| attitude error, gyro bias | 6 | OU-III, same transition and structured Qd |
| horizontal `[v_x, v_y, p_x, p_y, S_x, S_y, a_x, a_y]` | 8 | OU-III integrated OU chain and `S = 0` update, unchanged |
| vertical `[h, v_z, a_z, latent]` (+ `S_z`) | N (+1) | CARMA |
| accelerometer bias | 3 | OU-III residual-bias OU, gate and freeze |

CARMA(5,2) without `S_z` has 22 states (OU-III: 21).  The vertical
acceleration enters the attitude-aware accelerometer model exactly where
OU-III's `a_w,z` does,

    f_b = R_wb ([a_x, a_y, a_z]^T - g e_z) + b_a(T),

so the Jacobian with respect to the vertical block is `R_wb e_z e_a^T`, next
to `-[f_b]x` for attitude and `I` for the bias; every cross-covariance is
retained (the gain uses the full `P C^T`).  When bias learning is gated off
the bias uncertainty is still marginalized into the innovation covariance and
the bias rows are frozen, as in OU-III.  The position and velocity
pseudo-measurement APIs act on `h` and `v_z`.  Not carried over: wind heel,
IMU lever arm and the commanded-parameter random-walk filter (none enabled
by the deployed wrapper).  The block-diagonal linear transition is applied
structurally rather than formed densely; besides being cheaper, the dense
11 x 11 loop of the CARMA(3,0) build faulted under GCC 13 -O3 on x86-64 (an
aligned vector load on a misaligned column; clean under clang, ASan and
UBSan).

**`S_z`.**  The CARMA(5,2) prior has `S_h(0) = c0 > 0`, so `S_z = int h dt`
is not stationary under it: near DC `h` is white with level `c0` and `S_z`
is a random walk.  No zero pseudo-measurement is carried over.  The optional
variant propagates `S_z' = h` exactly with the CARMA states and applies
`S_z = 0` every `T_S` as an explicit regularizer of that random walk: with
information-rate-equivalent noise `r^2 T_S` it closes a first-order loop at
`w_R = sqrt(c0/(r^2 T_S))`; placing `w_R tau = kappa_R` (OU-III's pole
constant 0.3627) gives `r_Sz = sqrt(c0) tau / (kappa_R sqrt(T_S))`.  It is
tested separately (C52S).  CARMA(4,1)'s zero at DC makes `S_z` stationary
and needs no such update.

**Wrapper.**  `SeaStateFusionFilter_CARMA` reuses the shared front end,
adaptation mechanics, proxy startup and magnetic acquisition, and OU-III's
per-sample schedule, magnetometer-service bias gate, tilt watchdog and Racc
inflation; the horizontal chain keeps OU-III's SpectralMSE `r_S` law.  The
same measured operating point that tunes OU-III scales the CARMA shape:
`kappa = f_tune / f_z,model`, `f_tune = c_tau/(2 tau)` (the model's
zero-crossing frequency is placed on the measured one), and
`sigma_a = max(0.05 m/s^2, noise floor, sigma_aw / c_sigma)` (the measured
band acceleration RMS, with OU-III's floor).  These mappings were fixed
before any Monte Carlo run.  `setFixedVerticalCarma(1, 0)` holds the shape at
its 0.29 Hz design (C52F); `setCarmaSigmaFloor()` exists only for the
mechanism ablation C52N.

## 5. Stochastic-process validation

`carma_process-sim` generates 11.6 h (20 ms) and 400 s (1 ms) realizations
with the library's float model; `tools/carma_process_validation.py` checks
them against quantities derived independently in Python (companion
realization, scipy Lyapunov, polynomial PSDs).

| check (11.6 h at 20 ms, 400 s at 1 ms) | CARMA(3,0) | CARMA(4,1) | CARMA(5,2) | pass if |
|---|---|---|---|---|
| largest pole real part (rad/s) | -0.648 | -0.201 | -0.170 | < 0 |
| sample / Lyapunov variance, `h`, `v`, `a` | 0.992, 0.993, 0.999 | 1.006, 1.003, 0.999 | 1.009, 1.006, 1.002 | within 8 % |
| Welch vs analytic PSD, mean abs dB, `h`, `v`, `a` | 0.10, 0.10, 0.12 | 0.10, 0.11, 0.09 | 0.13, 0.13, 0.07 | < 0.35 dB |
| `S_v/S_h` vs `w^2`, `S_a/S_h` vs `w^4`, mean abs dB | 0.005, 0.031 | 0.017, 0.005 | 0.015, 0.017 | < 0.3 dB |
| `(dv/dt - a)` residual / Brownian prediction `g sqrt(dt/12)` | 0.998 | 1.001 | 1.000 | within 5 % |
| `E[da^2] / (g^2 dt)` | 0.984 | 1.000 | 0.994 | within 5 % |
| increment scaling exponent, `a` (1 = Brownian), `v` (2 = differentiable) | 0.98, 1.99 | 1.00, 2.00 | 0.99, 2.00 | within 0.1 |
| max `|h| / sigma_h` over 11.6 h | 4.58 | 4.33 | 4.49 | finite |

So `h` and `v` are differentiable, `a` is continuous with Brownian
increments (not differentiable), the spectra and variances are the designed
ones, and long runs stay bounded.  `reports/results/carma_benchmark/process_validation.json`
holds the full output.

![CARMA(5,2) realization against analytic spectra](../reports/results/carma_benchmark/carma52_process_psd.svg)

## 6. Benchmark protocol

**Truth.**  `tests/kalman_carma/CarmaBenchTruth.h`: random-phase synthesis
(160 jittered frequencies over 0.55-4.5 `f_p` x 5 directions, cos^2s
spreading) of the incident sea.  Vessel motion is the response of a small
rigid hull: each RAO is a second-order mass-spring-damper response
`wn^2 / (wn^2 - w^2 + 2 zeta wn i w)` to the wave excitation times the
hull-averaging factor `exp(-k l)` (wavenumber `k`, effective averaging
length `l`) with which a hull of finite draft and size filters short waves.
Heave: 0.55 Hz, `zeta` 0.35, `l` 0.5 m.  Roll and pitch are driven by the
beam and bow slopes: 0.40 Hz / 0.15 / 1.0 m and 0.60 Hz / 0.30 / 1.5 m, gains
0.8 and 1.0, plus heel/trim and a slow yaw wander; surge and sway follow the
orbital motion.  All channels are analytic harmonic sums, so position,
velocity, acceleration, attitude and rate come from one history.  No CARMA
model is used.  Incident elevation and vessel heave are reported
separately; errors are against vessel heave, normalized by the incident
`H_s` (the repository's convention).

The motion levels were calibrated against the pinned 28 ft vessel-RAO
dataset, before any filter was run on the new truth: vertical acceleration
0.3-1.5 m/s^2 RMS with peaks below 6.5 m/s^2, roll 3-5 deg and pitch 2-5 deg
RMS with peaks up to about 20 deg (30 deg in the large-rotation case), body
rates below 0.8 rad/s (1.3 rad/s), against 0.24-1.72 m/s^2, peaks up to
6.6 m/s^2, roll 1-8 deg RMS and rates up to 0.9 rad/s in the dataset.  Every
sea respects a realistic steepness, `H_s/L_p` <= 0.06 (deep-water peak
wavelength `L_p`).  The generated table "Realized motion levels" lists every
scenario (`carma_bench --level truth`).

**Sensors.**  The shared harness models (`make_imu_noise_model`,
`make_mag_noise_model`, the parameters of `process_wave_file_for_tracker`):
accelerometer 1.51 mg white, +-8 mg turn-on bias, 5e-4 m/s^2/sqrt(s) bias
walk; gyro 0.00157 rad/s, +-0.05 deg/s, 1e-5; magnetometer 0.8 uT at 25 Hz
with residual bias, scale, cross-axis and misalignment; 200 Hz IMU.  One
corrupted stream per (scenario, seed) is fed to every filter.

**Levels.**  Level 1 isolates the motion model: the true attitude levels the
measured specific force exactly, `y = a_z + b_z + n`, and reduced vertical
filters run on it -- OU-III's chain `[v, p, S, a_w] + b` (shared
`IntegratedOUChain`, SpectralMSE `r_S`, tau-scaled `S` cadence, shipping
constants) and the CARMA models `+ b`.  One measurement-only operating point
(the shipping `WavePeriodEstimator`, sigma band and `SeaStateAutoTuner` on
`y`) drives all of them; `--params oracle` instead takes `T_z` and `sigma_a`
from the realized truth.  Level 2 runs the complete shipping
`SeaStateFusion_OU_III` (unchanged, simulator configuration of
`kalman_ou_iii-sim` minus its dataset-specific direction/low-wave RAO
settings, which the firmware does not set) against `SeaStateFusion_CARMA`
with identical configuration.  The pinned vessel-RAO dataset
(oceanography-waves-lib v1.2.3) is replayed through the shared
`W3dSimulationRunner` with the unchanged `kalman_ou_iii-sim` and `carma-sim`
(same adapter configuration).

**Filters.**  OU3 (shipping OU-III), C30, C41, C52 (adaptive scale and
amplitude), C52S (with `S_z`), C52F (frequency scale fixed at the 0.29 Hz
design, amplitude adaptive), C52N (amplitude floor = noise floor, quiet
scenarios only, post-hoc mechanism ablation), C52Mv (Level 1 oracle-shape
diagnostic, see section 8).

**Scenarios.**  JONSWAP `f_p` 0.20 Hz (`H_s` 1, 2 m), 0.29 Hz (0.4, 0.8,
1.1 m; `gamma` 1.5/3.3/7 at 0.8 m), 0.45 Hz (0.2, 0.4 m) and a long 0.10 Hz,
4 m sea; Pierson-Moskowitz 0.15 Hz/2.5 m, 0.25 Hz/1 m, 0.40 Hz/0.4 m; tapered
broadband 0.10-0.70 Hz, 1 m; mixed swell (0.09 Hz, 1 m) and wind sea
(0.30 Hz, 0.8 m); bimodal (0.13 Hz/1.2 m + 0.42 Hz/0.4 m, 70 deg apart);
`f_p` rising 0.22 -> 0.36 Hz over 300-900 s (0.7 m); a 10 s change from
0.35 Hz/0.5 m to 0.18 Hz/2 m at 600 s; 5 cm `H_s`; still IMU (shared bias
model and post-calibration +-0.01 m/s^2); extra accelerometer drift (5 mg
ramp plus 2 mg oscillation); 15 deg heel and 5 deg trim; large rotation
(roll and pitch gains 2.0, +-15 deg yaw oscillation); +-25 % sample-interval
jitter; a no-RAO control (heave = incident elevation); 3600 s records.  100
seeds for the four primary cases (JONSWAP 0.29 Hz/0.8 m, PM 0.25 Hz/1 m,
mixed, changing `f_p`) at both levels, 30 for the rest, 20 for the 3600 s
records.

**Metrics.**  Scored on `t >= 300 s` of 1200 s records (and the last 60 s):
heave RMS, RMS/`H_s`, mean, peak; velocity and acceleration RMS; low-frequency
error (RMS of the 60 s moving average); linear drift trend; convergence time
(after Live, the 30 s rolling RMS stays below 1.5x the scored RMS); mean NIS
of the accelerometer innovation and the fraction inside its 95 % chi-square
interval; NEES of `[p_z, v_z]` (expected 2); failures (non-finite,
negative variance, `|e| > 10 H_s`); host update time and object size.
Paired statistics over seeds as described in the generated tables.

@@RESULTS_GRID@@

## 8. Mechanism: spectral shape, not observability

`tools/carma_mismatch_analysis.py` computes, for the true vessel-heave
spectrum of a case (the band-limited incident spectrum times the truth's
heave RAO) and the oracle operating point, the exact steady-state transfer
`T(f)` from the accelerometer to `h_hat` of each filter (scipy DARE on the
200 Hz discrete model; OU-III's `S = 0` update included at its
information-rate-equivalent variance) and splits the error variance into
signal distortion `int |1 + w^2 T|^2 S_h df` (what the prior does with the
true spectrum) and noise propagation `int |T|^2 N df` (measurement
observability):

| case | OU-III | CARMA(3,0) | CARMA(4,1) | CARMA(5,2), JONSWAP shape | CARMA(5,2), shape refitted to the truth |
|---|---|---|---|---|---|
| JONSWAP 0.29 Hz, 0.8 m | 13.4 (11.2 / 7.4) | 36.1 | 10.2 (9.7 / 3.2) | 13.7 (13.5 / 1.8) | 8.6 (7.4 / 4.4) |
| PM 0.25 Hz, 1 m | 17.3 (15.1 / 8.4) | 49.6 | 14.8 (14.3 / 3.7) | 21.3 (21.2 / 2.0) | 11.0 (9.9 / 5.0) |
| JONSWAP 0.20 Hz, 1 m | 19.6 (15.8 / 11.6) | 44.6 | 13.3 (12.1 / 5.6) | 16.6 (16.3 / 3.3) | 12.5 (11.1 / 5.6) |
| JONSWAP 0.45 Hz, 0.4 m | 6.6 (5.4 / 3.8) | 15.9 | 5.2 (4.9 / 1.6) | 7.0 (6.9 / 0.9) | 3.0 (2.5 / 1.7) |
| mixed swell + wind sea | 60.7 (59.0 / 14.4) | 128.0 | 93.7 (93.4 / 6.8) | 302.1 (302.0 / 3.8) | 129.8 (129.8 / 2.1) |
| bimodal | 31.6 (26.6 / 17.1) | 74.3 | 26.9 (25.6 / 8.4) | 51.3 (51.1 / 4.8) | 18.9 (17.3 / 7.6) |

Steady-state heave RMS error in mm (signal distortion / noise propagation).
Two effects compete.  OU-III's prior is nearly shape-free above its corner
`1/tau`, so its signal distortion is moderate everywhere, but it lets
accelerometer noise through at low frequency (7-17 mm, mostly below
0.1 Hz).  A stationary CARMA prior propagates two to four times less noise,
and wins whenever its shape is close enough to the vessel's heave spectrum;
when it is not, the filter attenuates or mis-phases the true motion wherever
the prior is wrong, and signal distortion dominates.  CARMA(4,1) (one zero
at DC, two resonances) sits on the right side of that balance for every
single-peaked sea and the bimodal case; the fixed JONSWAP-fitted CARMA(5,2)
is level with OU-III on JONSWAP and worse on PM, whose broader peak its
narrow first resonance (`zeta_1` = 0.097) does not cover.  For the mixed sea
the 0.09 Hz swell falls near CARMA(5,2)'s numerator zero once `kappa` scales
the shape to the wind sea: 289 mm RMS of its 302 mm error is below 0.1 Hz.
The constrained refit (last column) cannot represent two well-separated
peaks either and is still worse than OU-III there.

@@MATCHED@@

## 9. Does CARMA replace `S_z`?

@@SZ@@

## 10. Cost and ESP32-S3 build

@@COST@@

ESP32-S3 (`esp32:esp32:m5stack_atoms3`, arduino-esp32 3.3.11, xtensa GCC
14.2.0, `-Os`): every variant compiles and links, with no relocation or
literal-pool errors (the ESP32-S3 failure mode mentioned in the request).

| variant | program flash (bytes) | wrapper object | filter object | static RAM (globals) |
|---|---|---|---|---|
| OU-III (shipping baseline) | 683 246 | 3 032 | 5 692 | 22 400 |
| CARMA(5,2), 22 states | 671 226 | 3 000 | 5 224 | 22 400 |
| CARMA(4,1), 21 states | 670 118 | 3 000 | 4 760 | 22 400 |
| CARMA(5,2) + `S_z`, 23 states | 675 854 | 3 000 | 5 656 | 22 400 |

(object sizes from `sizeof` in the Xtensa build; the sketch allocates the
fusion object on the heap, so the static RAM reported by the linker is the
same for all four.)  Caveats: built against Eigen 3.4.0
through a local `ArduinoEigenDense.h` shim rather than the published
ArduinoEigen package, which the build host could not download; the M5Unified
library was staged without its AtomS3R sources; prototypes generated with
exuberant-ctags.  **No AtomS3R hardware was available: on-target update
time, stack use and numerical behaviour in float on the Xtensa FPU are not
verified.**  The sketch prints mean and maximum `micros()` per update when
run on hardware.

## 11. Engineering assessment

@@ASSESSMENT@@

## 12. Limitations and open items

* **Shape identification.**  The only CARMA(5,2) configuration that beats
  OU-III uses a shape refitted to the true spectrum.  An online identifier
  for the seven shape parameters, its excitation requirements and its effect
  on stability and covariance consistency are not studied.
* **Overconfidence under mismatch.**  With a wrong shape the CARMA filters'
  covariance does not cover their error (NEES up to about 70 against 2 on the
  mixed sea).  Nothing inflates it; this is reported, not fixed.
* **Quiet water.**  @@QUIET_LIMIT@@
* **Startup.**  OU-III and the adaptive CARMA filters share the startup and
  operating-point path; on 1-2 % of seeds both show a multi-metre startup
  excursion before Live settles (section 7).  It is not CARMA-specific and is
  not investigated here.
* **Not carried over**: wind-heel feed-forward, IMU lever arm and the
  parameter random-walk option of OU-III (all off in the deployed wrapper).
* **Truth model.**  One synthetic hull (heave RAO 0.55 Hz, `zeta` 0.35) and
  one pinned dataset.  The ranking of the priors depends on how far the
  vessel's heave spectrum is from the JONSWAP shape, so it can differ for
  other hulls.
* **No stability claim.**  None of this is part of, or evidence for, the
  OU-III stability study; the CARMA prior is a statistical model and no
  admissibility or practical-stability statement is made for it.
* **Hardware.**  ESP32-S3 compile and link only (section 10).

## Reproduction

    make -C tests/kalman_carma build test            # unit, process and smoke checks
    python3 tools/carma_fit.py                       # spectral fits
    (cd tests/kalman_carma && ./carma_process-sim process)
    python3 tools/carma_process_validation.py --input tests/kalman_carma/process \
        --check --out reports/results/carma_benchmark
    python3 tools/carma_benchmark.py run             # synthetic grid (resumable)
    make ensure-sim-data && make -C tests/kalman_ou_iii build
    python3 tools/carma_benchmark.py dataset         # pinned dataset replay
    python3 tools/carma_benchmark.py traces
    python3 tools/carma_mismatch_analysis.py
    python3 tools/carma_benchmark.py report

ESP32-S3 build: `tests/kalman_carma/arduino/carma_compile_check/` with
`-DCARMA_CHECK_VARIANT=0..3` (OU-III, CARMA(5,2), CARMA(4,1), CARMA(5,2)+`S_z`)
for `esp32:esp32:m5stack_atoms3:CDCOnBoot=cdc,USBMode=hwcdc`.  It is under
`tests/`, so the shipped-sketch CI jobs do not compile it.
