# Shared-modes QMEKF (`KalmanCarmaSharedModes`, experimental)

**Status: experimental, simulation only, not a default.**  The shipping
filter stays `SeaStateFusionFilter_OU_III`; the independent-axis CARMA3 QMEKF
(`src/kalman_carma/`, benchmark name `A41`) is unchanged and remains the
baseline.  Nothing here touches TFG or TFG-CARMA3.

This filter promotes Option C of the coupled-prior study
([carma-coupled-three-axis.md](carma-coupled-three-axis.md)) -- shared
directional stochastic modes -- into a standalone QMEKF with 27 error
states.  It adds the coupling confidence and fallback that the study showed
to be necessary.

* Core: `src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h`,
  `KalmanCarmaSharedModes<T = float>`, the QMEKF (scalar `T`, as every
  Kalman core in the repository), with `CarmaCouplingConfidence.h` and the
  motion block `CarmaSharedModesMotion.h` (`SharedModesMotion<T>`).
* Wrapper: `SeaStateFusionFilter_CarmaSharedModes.h`,
  `SeaStateFusion_CarmaSharedModes`, the unchanged CARMA3 schedule with the
  core as its inner filter (float, as every wrapper).
* Tests and simulators: `tests/kalman_carma_shared_modes/`
  (`carma_shared_modes-test`, the filter; `carma_shared_modes_realization-test`,
  the motion block; `kalman_carma_shared_modes-sim`, the dataset replay;
  `kalman_carma_shared_modes-bench`, the paired Monte Carlo;
  `kalman_carma_shared_modes-timing`, host core timing, on demand; the
  ESP32-S3 probe `arduino/csm_size_probe`).
* Study: `tools/carma_shared_modes_analysis.py` (mathematics),
  `tools/carma_shared_modes_study.py` (Monte Carlo, robustness, traces,
  dataset, gates, timing, size, report) ->
  `reports/results/kalman_carma_shared_modes/`; deterministic table
  `tools/carma_shared_modes_sim_table.py`.
* Articles: `doc/kalman_carma_shared_modes/`
  ([`kalman_carma_shared_modes.tex`](../doc/kalman_carma_shared_modes/kalman_carma_shared_modes.tex),
  the standalone presentation of the 27-state method: CARMA motivation,
  absence of an integral pseudo-measurement, spectral fitting, shared modes,
  state reduction, scaling, adaptation, proxy, period estimation, magnetic
  aiding and startup;
  `kalman_carma_shared_modes-27state.tex`, the numerical system-validation article;
  `kalman_carma_shared_modes-math.tex`, the model;
  `kalman_carma_shared_modes-minimal.tex`, the minimal realization and the
  numerical study; `kalman_carma_shared_modes-charts.tex`, the replay chart
  gallery).

**Findings** (details in section 6; every number from the committed study data):

* **27 states.**  The shared-mode law is one complexified CARMA(4,1) mode
  per axis; realified it has eight real states, and exactly two of them are
  unobservable for every excitation and attitude.  The filter carries the
  six observable ones (section 2).  It gives the same law of `p, v, a` on all
  three axes, every co- and quadrature spectrum included, and with the
  coupling off it is the CARMA3 baseline A41 in law, model commits included.
* **Accuracy** (paired Monte Carlo, 26 seas x 30 seeds, against A41 with the
  same Evidence excitation, `A41@C`): the 3-D position RMS is lower on 25 of
  26 seas, by 0.1 to 15.8 % (median 3.2 %; swell 15.8 %, regular wave 8.4 %,
  spread JONSWAP 3.3 %).  Heave is better on 24 seas (spread JONSWAP 5.9 %);
  it is worse under manoeuvring (-6.1 %, the confidence withdraws the
  coupling during turns; 3-D still +0.6 %) and on the one losing sea.
* **The one loss: over-coupling.**  With the hull's heave RAO halved the
  filter is 58 % worse than A41 in 3-D and overconfident (NEES 14): the
  excitation mapping keeps the measured regression with the prior's
  horizontal-to-vertical ratio, and the confidence cannot see it (section 6).
* **Consistency.**  Median position NEES 1.65 (A41 1.80), heave 95 %
  coverage 0.91 (0.92), no failed run in 1560.
* **Wrong coupling.**  Every deliberately wrong prior (reversed, turned
  90 deg, phase-shifted, doubled) is withdrawn by the confidence: the filter
  then has A41's error.  Without the confidence a doubled prior costs 122 %
  (57 against 26 mm, NEES 14).
* **Dataset replay** (eight scored records, gated by its own bars in
  `run_tests.sh`): against TFG-CARMA3, heave 20-25 % and 3-D 7-9 % lower
  (the heave advantage is the Evidence excitation; the 3-D one the
  coupling); pitch 4 % and the 3-D biases 2-3 % higher.
* **Embedded cost** (against A41): 7,992 instead of 4,248 multiply-adds per
  sample, 12.5 instead of 6.1 kB filter object, 33 kB more code on the
  ESP32-S3; host core time 1.6-1.7 times A41's.

## 1. What is reused, what is new

| part | source | status |
|---|---|---|
| QMEKF core: attitude error, gyro and accelerometer bias, Joseph updates, multiplicative reset, magnetometer | `kalman_carma/Kalman3D_Wave_CARMA3.h` (`Motion` parameter) | reused unchanged |
| coupling estimate and excitation mapping | `kalman_carma/CarmaCrossSpectrumEstimator.h` | reused unchanged |
| schedule: proxy startup, magnetic acquisition, tuner, 0.1 s commits | `kalman_carma/SeaStateFusionFilter_CARMA3.h` | reused unchanged |
| shared-mode motion block in its minimal realization (six states per axis, Hermitian excitation, paired discretization), with the short-horizon coupling monitor | `CarmaSharedModesMotion.h` | new |
| per-axis re-coordinatization at a model commit | `Kalman3D_Wave_CARMA3::recoordinate_motion_` (`take_axis_recoordinatization`) | new hook, instantiated only for a motion block |
| coupling confidence, hysteresis, smoothing | `CarmaCouplingConfidence.h` | new |
| PSD-preserving decoupling of the covariance | `Kalman3D_Wave_CARMA3::shrink_motion_vertical_coupling` | new member, instantiated only for a motion block |
| QMEKF core `KalmanCarmaSharedModes<T>`: confidence in `time_update` / `measurement_update_acc_only`, diagnostics, direction estimate | `KalmanCarmaSharedModes.h` | new (derives from the CARMA3 core) |
| float wrapper `SeaStateFusion_CarmaSharedModes` | `SeaStateFusionFilter_CarmaSharedModes.h` | new (the CARMA3 wrapper with the new core) |

With the confidence held at zero the filter is `CouplingMode::Off`, bit for
bit, and with the coupling off it is the CARMA3 baseline A41 in law, model
commits included (section 2).

## 2. State and model

Error state (27): `[dtheta(3), b_g(3), x(6), y(6), z(6), b_a(3)]`, NED,
attitude error in the body frame.

**The law.**  Per axis the motion is the real part of a complexified
CARMA(4,1) mode of the C41J shape,

    zeta = B(s) / (P+(s) D(s)) xi

* `P+` holds the upper-half-plane poles of the C41J polynomial `P`.
* `D` is the real Hurwitz quadratic with
  `|D(iw)|^2 = (|P+(iw)|^2 + |P+(-iw)|^2) / 2` (roots -7.93 and -1.25 rad/s
  at kappa = 1).
* `xi` is circular complex white noise.
* `Re zeta` has exactly the C41J law, and the noise enters `a` only
  (relative degree 3).

The axes share the dynamics and are driven with a Hermitian PSD excitation
`Sigma`.  The displacement cross-spectrum is
`S_ij(w) = (g^2/2) [Sigma_ij |H(iw)|^2 + conj(Sigma_ij) |H(-iw)|^2]`.  Its
co-spectrum is `Re Sigma` and its quadrature spectrum `Im Sigma`, which is the
90 deg x-z relation of linear waves, shifted by the heave RAO.  Derivations,
exact discretization, the complex-OU factorization (Fasen-Hartmann and
Scholz) and the MCARMA context (Stelzer's 2024 correction) are in the
mathematical article.

**The realization: six states per axis.**  Realified, the mode has eight
real states `[Re(p, v, a, l), Im(p, v, a, l)]`, and two of them are
unobservable for every excitation and attitude: each root `r` of `D` is an
eigenvalue of both the complex mode and its conjugate, so the realified
dynamics have it twice, and one scalar output cannot resolve it.  The
unobservable subspace

    U = {0} x (span_C{w(r3), w(r4)} intersected with R^4),   w(r) = [B(r), r B(r), r^2 B(r), 1]

lies in the imaginary half and does not depend on `Sigma`, the attitude or
any QMEKF operation.  The filter carries the observable six states.  With
`P+ = N_r + i N_i` (real polynomials),

    Re zeta = B/P * eta,   eta = V [xi_r; xi_i],   V = [N_r, N_i] / D,

and `|D|^2 = |N_r|^2 + |N_i|^2` makes `V` lossless (co-inner).  Per axis the
six states are the CARMA3 baseline's own `[p, v, a, l]`, driven at `a` by
`eta = xi_r + C_w w`, and the two states `w` of `V` in input-normal
coordinates (`A_w + A_w^T = -B_w B_w^T`, `C_w = -B_w e_1^T`, `B_w` upper
triangular).  `p, v, a` are states, so every accelerometer row is a unit
selector rotated by `R_wb`.  The axes share the dynamics, `F = I_3 (x) F_6`,
and the Hermitian excitation enters through two intensities:

    Q_ab = Re(Sigma_ab) Q_R + Im(Sigma_ab) Q_I

`Q_R` symmetric from `g2 G G^T`, `Q_I` antisymmetric from `g2 G J2 G^T`,
`J2 = [[0, -1], [1, 0]]` (a paired discretization of the CARMA3 core's
scaled series, doublings and jitter series).  Every quantity scales exactly
with `kappa` (`V_kappa(s) = V_1(s / kappa)`), so one closed-form 2 x 2
construction at `kappa = 1` is the only factorization.  This realization and
the realified mode give the same law of `p, v, a` on all three axes, every
co- and quadrature spectrum included, and at fixed parameters the same
estimates, innovations and covariance through every amplitude, coupling,
confidence and fallback change (the minimal-realization article, with
proofs).

**Excitation** (deployable).  Each sample, the filter complex-demodulates its
own world-frame acceleration `R_wb^T (f - b_a) + g e_z` at the operating
frequency `f_op`, then low-passes the envelope at `0.15 f_op` and averages
`A A^H` over 240 s.  The result is the band cross-spectrum `S_hat`.  The
excitation keeps the baseline diagonal and adds only the horizontal-vertical
terms:

    b = s gain c S_hat_hz,   c = sigma_z^2 / S_hat_zz,   s = max(0, 1 - 1/(N_eff gamma^2))

Here `gamma^2` is the debiased multiple coherence of z on (x, y).  If
`b^H H^-1 b` would exceed `rho sigma_z^2` (`rho = 0.8`,
`SharedModesConfig::cap_coherence`), `H` is enlarged along the coupling
(`CouplingStatus::cap_gain_x/y` report the variance factors),
`H <- H + ((lambda - 1) / q) b b^H` with `q = b^H H^-1 b` and
`lambda = q / (rho sigma_z^2)`.  This keeps the measured regression and the
cross-wave prior; `cap_coherence = 0` restores the earlier cross-term
shrinkage to `0.999 sigma_z^2`, which understates the regression under a weak
heave RAO.  The study bundle in `reports/results/kalman_carma_shared_modes/`
predates this policy; `reports/results/tfg_shared_modes/` compares both.  The excitation is
PSD by construction, equivariant under rotations about the vertical, and the
baseline diagonal without coupling evidence.  `gain` is the coupling
confidence `c in [0, 1]` of section 4.

**Sea-state adaptation.**  This is the CARMA3 wrapper's own adaptation,
unchanged.  Every 0.1 s it commits `kappa = f_tune / f_z,base` (one frequency
scale for every axis) and the isotropic amplitudes; a commit rebuilds `F`, `Q`
and the excitation.  The posterior covariance is reset to the stationary
prior only on entering Live, never during adaptation.  A commit holds
`[p, v, a, l]`, as the baseline does, and the normalized colouring
`w / s_a`, `s_a = sqrt(g2 sigma_a^2)` (`w <- (s_new / s_old) w` per axis, the
covariance through the same congruence).  `V` is lossless with a constant
Gramian, so `eta` stays exactly white under any change of `kappa` and
`sigma`: **with the coupling off the filter is the CARMA3 baseline A41 in
law, adaptation included** (`carma_shared_modes_realization-test` [3]:
4e-13 relative in double precision through 10 Hz kappa and amplitude
commits).

**Wave excitation: Evidence mode.**  The filter uses the CARMA3 wrapper's
Evidence excitation mode, the still-sea phantom fix that the AtomS3R sketches
deploy (noise-debiased band variance and a detection weight instead of the
Legacy 0.05 m/s^2 floor and stillness attenuation,
[carma-wave-excitation.md](carma-wave-excitation.md);
`SeaStateFusion_CarmaSharedModes::begin`, `kalman_carma_shared_modes-sim`).
Its regression bars are fitted in this mode, and the Monte Carlo baseline is
A41 in the same mode (`A41@C`).

## 3. Measurement update and reset

These are the baseline's, line for line: the accelerometer
`f = R_wb (a - g e_z) + b_a` (Joseph form, the motion columns are
`R_wb^T e_i` on the `a` selector of each axis), the magnetometer, the
multiplicative attitude reset (`dtheta -> 0`, `P <- G P G^T` on the attitude
block and its cross-covariances), and the bias states.  The motion block
enters the cross-covariances through the same congruence as the baseline's
CARMA blocks.

## 4. Coupling confidence and fallback

The research study found two kinds of failure:

* **Confident wrong coupling.**  On data generated by the model itself, a
  confidently reversed coupling gave NEES 380, against 6.8 for the correct
  coupling.
* **Deployable failures.**  In the deployable 50-seed comparison, the
  estimated coupling lost on manoeuvring (heave -12.7 %) and on wave
  onset/cessation (-2.9 % horizontal, -3.5 % heave), while helping on every
  steady directional sea.

`CarmaCouplingConfidence` turns the available evidence into one number
`c in [0, 1]`.

**Evidence.**

1. *Consistency of the prior with recent data.*  The prior's regression
   `beta_p = S_hz / S_zz` (the 240 s estimate or an external prior) is
   compared with that of a 60 s monitor estimate fed with the same
   world-frame accelerations:

       D^2 = sum_i |beta_f,i - beta_p,i|^2 / (v_f,i + v_p,i + (eps |beta|)^2),
       v_i = (1 - gamma_iz^2) S_ii / (N S_zz)

   `v_i` is the large-sample variance of a transfer-function estimate from
   `N` independent envelope samples.  With a consistent prior `E[D^2] ~ 2`
   (two complex components), and smaller still because the two estimates
   overlap.  The decision is inconsistent above 16 and consistent again below
   6 (hysteresis).
2. *Wave energy.*  The monitor's vertical band power is compared with the
   accelerometer noise floor of the same envelope filter,
   `sigma_n^2 dt pi b f / 2`.  There is no coupling at or below an SNR of 10
   and full coupling at or above 100, log-linear in between.
3. *Manoeuvring.*  The low-passed world-frame yaw rate (2 s) above 3 deg/s
   withdraws the coupling, and the withdrawal is held for 30 s afterwards
   while the monitor refills.
4. *Abrupt sea-state change.*  `|d ln f_op / dt|` (10 s) above 0.01/s
   withdraws it the same way, also with a 30 s hold.
5. *Innovation consistency.*  The coupling falls back when the 20 s mean
   accelerometer NIS exceeds twice its own 600 s mean, and recovers below
   1.5 times it.  The test is self-normalized because the configured R is
   inflated, so the absolute NIS level carries no information.
6. *Warm-up.*  There is no coupling until both the prior and the monitor
   exist (60 s for the monitor, 120 s for the estimate).

**Smoothing.**  The evidence `e` (the product of the gates, so 0 when the
prior is inconsistent) drives

    c <- c + (1 - exp(-dt/tau)) (e - c),   tau = 3 s falling, 60 s rising

so `|dc/dt| <= 1/tau` and `c` never steps (test [4]).

**Applying c.**  Every 0.1 s a changed `c` sets `gain = c` and the excitation
is rebuilt at the next step, with or without a model commit, so the
standalone core engages and withdraws the coupling on its own (test [9]).
`c = 0` is the independent-axis law exactly (test [2]).  The state is never
touched, so position, velocity and attitude are continuous through every
change (test [6]).

**Removing stale cross-covariance.**  When `c` falls, the covariance still
carries horizontal-vertical correlation that the reduced prior no longer
sustains.  Zeroing the `P_hz` blocks is not positive semidefinite in general:
both blocks also correlate with the attitude and bias states `r`.  The
fallback instead shrinks the part of `P_hz` that `r` does not explain:

    P_hz <- E + alpha (P_hz - E),   E = P_hr P_rr^-1 P_rz,   alpha = c_new / c_old

In the Schur complement `C = P_mm - P_mr P_rr^-1 P_rm` of the motion block,
this is `C <- alpha C + (1 - alpha) blockdiag(C_hh, C_zz)`, a convex
combination of PSD matrices, while `P_rr` and `P_rm` are unchanged.  `P`
therefore stays positive semidefinite (proof in the mathematical article,
test [3]), and the correlation induced by the shared attitude and bias errors
is kept.  A decline below 0.1 % per commit accumulates until it reaches
0.1 %, so slow declines also shrink (test [9]).  The cost is one 9 x 9 LDLT
and a 12 x 9 x 6 product, at most at 10 Hz and only while `c` falls.

**Attitude re-initialization.**  The band cross-spectrum, the monitor and the
confidence are world-frame statistics of the attitude they were measured
with.  Every `initialize_from_*` call therefore restarts them (a research hold
from `forceCouplingConfidence` survives); the coupling returns after the
usual warm-up.  A non-finite gyro sample is rejected by the core, and the
CARMA wrappers drop any sample with a non-finite rate or specific force.

**Recovery.**  Recovery is automatic.  Once the evidence returns, `c` rises
with `tau = 60 s` and the cross-covariance rebuilds through the dynamics and
the measurements (test [6]: a reversed prior corrected at 450 s returns to
`c > 0.99`).

**What the confidence cannot do.**  It cannot detect a wrong coupling that
the IMU data cannot contradict.  The consistency test compares the prior
with the same IMU's recent band statistics: a heading error, for example,
rotates both the prior and the monitor, so it shows up as heading error and
not as a coupling mismatch.  It also reacts only after the evidence has
accumulated (seconds to a minute), and it cannot see the over-coupling of the
excitation mapping (section 6).  The robustness results in section 6
quantify these limitations.

## 5. Public API

The core is a Kalman core like the others, templated on the scalar:

```cpp
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"
using ocean_imu::carma_shared_modes::KalmanCarmaSharedModes;   // <T = float>

KalmanCarmaSharedModes<float> f(sigma_a, sigma_g, sigma_m);    // the CARMA3 core's constructors
f.configureSharedModes(cfg);   // optional: mode, confidence thresholds, horizon, band
f.time_update(gyr, dt);        // as Kalman3D_Wave_CARMA3, plus the confidence step
f.measurement_update_acc_only(acc, temp_c);
f.measurement_update_mag_only(mag);
f.get_position(); f.get_velocity(); f.get_world_accel(); f.quaternion_boat();
f.coupling();       // mode, confidence c, state (warm-up/coupled/reduced/fallback),
                    // D^2, SNR, NIS ratio, applied coherence / phase, shrink events
f.waveDirection();  // direction of travel (deg from the filter's north), uncertainty, phase
f.health();         // NIS, covariance finite / symmetric / min diagonal
```

Behind the CARMA3 schedule (startup, magnetic north lock, tuner, model commits):

```cpp
#include "kalman_carma_shared_modes/SeaStateFusionFilter_CarmaSharedModes.h"
ocean_imu::carma_shared_modes::SeaStateFusion_CarmaSharedModes fusion;
fusion.begin(cfg);                          // as SeaStateFusion_CARMA3; shared-modes defaults
fusion.begin(cfg, shared_modes_cfg);        // ... or an explicit shared-modes configuration
fusion.update(dt, gyr, acc, temp_c); fusion.updateMag(mag);
fusion.raw().mekf().coupling();             // the core
```

`CouplingMode::Estimated` is the default and the only deployable mode.
`External` (with `setExternalCoupling`) is for research priors and the
robustness tests; `Off` gives independent axes; `forceCouplingConfidence` is
a test hook.  Ground truth and RAO knowledge never enter the classes.

## 6. Validation

`carma_shared_modes-test` (46 checks) and `carma_shared_modes_realization-test`
(18 checks), both run by `make all`:

| test | group | result |
|---|---|---|
| filter [1] | 27 states, bit-identical repeat runs and reset, finite symmetric covariance | yes |
| filter [2] | confidence disabled: the estimated coupling is applied; confidence 0 = Off | yes; bit-identical |
| filter [3] | fallback algebra: PSD, only P_hz changes, state untouched, alpha = 1 no-op, alpha = 0 leaves E | yes |
| filter [4] | confidence: consistent -> 1, reversed -> fallback < 0.05 in 10 s, rate bounds, hysteresis, recovery, manoeuvre hold, sea change, low energy, NIS ratio, warm-up | yes |
| filter [5] | direction of travel and phase from a regression | yes |
| filter [6] | synthetic directional sea: estimated coupling, reversed external prior with / without confidence, continuity through a forced fallback, recovery | see below |
| filter [7] | the core in double precision | finite, symmetric, coupling engages |
| realization [1] | 27 states, six per axis, Hurwitz, poles scale with kappa, colouring Gramian I | 2e-15; 3e-15 |
| realization [2] | against a test-local eight-state realified mode: noise -> (p, v, a) transfer functions; stationary three-axis law and lag cross-covariance (Hermitian excitation) | 3e-14; 3e-11 |
| realization [3] | coupling off = A41 through kappa and amplitude adaptation (double) | 4e-13 |
| realization [4] | p, v, a continuous through the core's commit path; transforms keep the p, v, a rows; PSD after 2000 commits | exact; exact; yes |
| realization [5] | float wrapper on a coupled sea: symmetric PSD covariance, coupling engages, position RMS below A41's | yes; 23.1 vs 27.5 mm |

Test 6 of the filter runs on an exact-kinematics directional sea, scored
after 300 s:

| run | horizontal RMS | heave RMS | NEES_p |
|---|---|---|---|
| independent axes | 23.1 mm | 13.1 mm | 1.92 |
| estimated coupling | 19.8 mm | 10.6 mm | 1.65 |
| reversed external prior, no confidence | 67.5 mm | 41.8 mm | 26.0 |
| reversed external prior, with confidence | 23.1 mm | 13.1 mm | 1.92 |

With the confidence, the reversed prior's run is the independent-axis result:
the coupling was never applied.

**Paired Monte Carlo** (`kalman_carma_shared_modes-bench`, carma3_bench
engine: random-phase directional seas through the benchmark hull's RAOs, one
IMU and magnetometer stream per row shared by every filter).  26 seas in five
groups (fixed; coupling strength and phase: spreading, heading, heave RAO
x0.5 / x1.5; gradual frequency change; rapid period, direction, height and
onset/cessation changes; disturbances: turns, manoeuvring, coloured and
variance-cycling accelerometer noise, +-25 % sample jitter, magnetic
disturbances with calibration residuals, accelerometer-bias drift).  30 seeds
each, 1200 s at 200 Hz, scored from 300 s, paired bootstrap over seeds.
Filters: `A41@C` (the CARMA3 baseline in the same Evidence excitation mode)
and `CSM`.  With the coupling off CSM is A41@C, so the paired difference is
the coupling.  Raw rows: `reports/results/kalman_carma_shared_modes/raw_runs.csv`;
tables and figures: `doc/kalman_carma_shared_modes/csm-*.tex-part` and
`figures/`.

Paired RMS change of CSM against A41@C by group (%, + = better; mean over
the seas of the group, range in brackets):

| group | heave | horizontal | 3-D |
|---|---|---|---|
| fixed (9) | +6.19 (+0.3 .. +17.4) | +4.11 (+0.1 .. +15.2) | +4.69 (+0.2 .. +15.8) |
| coupling (6) | -5.80 (-61.1 .. +7.9) | -7.90 (-57.1 .. +2.4) | -7.21 (-58.0 .. +3.9) |
| gradual (1) | +9.82 | +3.42 | +5.33 |
| rapid (3) | +5.69 (+4.7 .. +7.1) | +2.34 (+1.6 .. +3.7) | +3.32 (+2.6 .. +4.6) |
| disturbance (7) | +2.94 (-6.1 .. +5.9) | +1.45 (+0.1 .. +2.2) | +2.22 (+0.1 .. +3.4) |

The coupling group mean is dominated by heave RAO x0.5; its other five seas
are +1.6 to +4.0 % in 3-D.  Intervals are tight (95 % half-widths 0.1-0.5
points on the single-system seas).

* **Where the coupling helps most.**  The most coherent seas: swell (+15.8 %
  3-D, +17.5 % heave), the regular wave (+8.4 %, +13.2 %), PM (+5.2 %), the
  gradual frequency change (+5.3 %).  Least: the very low sea (+0.2 %) and
  coloured accelerometer noise (+0.1 %), whose 818 mm error is the shared
  low-frequency disturbance, not the wave model.
* **Manoeuvring.**  Heave -6.1 % [-7.1, -5.1] with 3-D +0.6 %: the confidence
  withdraws the coupling during the turns, and the withdrawal itself costs
  heave.
* **Transients.**  Around the frequency, period, direction and height
  changes and wave onset/cessation (10-seed ensembles, `traces/`), CSM's
  error envelope stays at or below A41's; during the direction change the
  two coincide (the coupling is withdrawn).
* **Consistency.**  Median position NEES 1.652 (A41 1.802), position-velocity
  NEES 5.80 (5.81), heave 95 % coverage 0.912 (0.918); no failure in 1560
  filter runs.
* **Wrong coupling** (external prior built from the truth, then reversed,
  turned 90 deg, phase-shifted 45 deg or doubled; 3 seas x 20 seeds).  With
  the confidence on, every wrong prior is withdrawn on every sea and CSM has
  exactly A41's error (25.7 mm on the spread JONSWAP); the true prior gives
  24.8 mm.  Without the confidence: doubled 57.0 mm (NEES 13.9), reversed
  29.4 mm, turned 27.1 mm, phase-shifted 24.9 mm; on the mixed sea reversed
  and turned priors cost 21-23 %.  A coupling frozen at 400 s in the keel
  frame through 90 deg turns: 25.4 mm (estimated: 25.2, A41: 25.7).
* **Over-coupling of the excitation mapping (heave RAO x0.5).**  CSM is
  58.0 % [51.8, 63.5] worse than A41 in 3-D and overconfident (median NEES
  14.3 against 2.8).  The measured multiple coherence is the same as on the
  unscaled sea (0.87), but the applied coherence rises from 0.13 to 0.38: the
  mapping keeps the measured regression `S_hz / S_zz`, the vertical prior and
  the horizontal prior, so the applied coherence is
  `gamma^2 (S_hh / S_zz) / (sigma_h^2 / sigma_z^2)`.  The consistency test
  compares two data-driven regressions and cannot see this.  A candidate fix
  is to cap the applied coherence at the measured one,
  `|b|^2 <= gamma^2 sigma_h^2 sigma_z^2`; it needs its own validation on the
  26 seas.

**Dataset replay and gates.**  `kalman_carma_shared_modes-sim` replays the
eight scored records of the pinned oceanography-waves-lib dataset with the
seeds, sensor model and adapter of every CARMA filter and is gated by CSM's
own bars (worst default-seed value plus about 0.5 %; `kCsmRegressionBars`):

| gate (worst of 8) | CSM | bar | TFG-CARMA3 (T41, Legacy) | T41, Evidence |
|---|---|---|---|---|
| heave % Hs, JONSWAP | 2.455 | 2.467 | 3.057 | 2.459 |
| heave % Hs, PM-Stokes | 2.403 | 2.416 | 3.215 | 2.44 |
| 3-D %, JONSWAP | 4.41 | 4.433 | 4.865 | 4.907 |
| 3-D %, PM-Stokes | 4.879 | 4.904 | 5.252 | 5.277 |
| roll (deg) | 0.3571 | 0.359 | 0.3551 | 0.3555 |
| pitch (deg) | 0.222 | 0.2231 | 0.2128 | 0.2131 |
| yaw (deg) | 0.7556 | 0.7595 | 0.7623 | 0.7631 |
| accel z bias % | 4.164 | 4.186 | 4.185 | 4.173 |
| accel 3-D bias % | 83.56 | 83.99 | 82.12 | 82.2 |
| gyro 3-D bias % | 12.64 | 12.71 | 12.33 | 12.34 |

Against CSM with the coupling off (`CSM_OFF`, which is A41 in Evidence mode)
the replay shows the Monte Carlo pattern on every record: heave 0.2-2.8 %
and 3-D 2-10 % lower, gyroscope bias equal, the 3-D accelerometer bias error
1.8-4.9 % higher on the seas of Hs 1.5 m and above (a cost of the coupling:
the two runs differ in nothing else).  The heave advantage over TFG-CARMA3 is the excitation mode
(with Evidence excitation TFG-CARMA3 is within 0.04 points of Hs); the 3-D
advantage is the coupling.  TFG-CARMA3's two-frame Lie group separates the
biases and pitch slightly better.

## 7. Embedded efficiency

| filter | states | P (B) | multiply-adds / sample | per 0.1 s commit | ESP32-S3 filter object (B) | code (B) | largest frame (B) | update chain (B) |
|---|---|---|---|---|---|---|---|---|
| A41 | 21 | 1764 | 4,248 | 4,608 | 6136 | 173,375 | 1328 | 4832 |
| CSM | 27 | 2916 | 7,992 | 19,440 | 12504 | 206,514 | 1888 | 4864 |

* Multiply-adds are the analytic counts of the CARMA3 core's structured
  implementation (`carma3_analysis.cost_model`, extended: Kronecker motion
  block, paired `Q_R`/`Q_I` additions, commit discretization and the dense
  per-axis re-coordinatization).
* ESP32-S3 sizes: Espressif GCC `esp-14.2.0_20241119`, `-Os`, arduino-esp32
  flags, object of the complete wrapper (`csm_size_probe/build.sh`).  The
  filter object holds `P` (729 floats) and the motion discretization (the
  paired 6-state series), the two coupling spectrum estimators and the
  confidence.  The wrapper holds its core through a pointer, hence equal
  wrapper sizes (3208 B).  The extra code is the paired discretization and
  the closed-form construction (generic Eigen decompositions there would cost
  44 kB).  The worst static call chain of one update is a lower bound
  (function pointers and library frames are not counted).
* **Host time** (x86-64 `-O3`, one job on an otherwise idle machine; ranking
  only, not an ESP32-S3 time).  The QMEKF core alone
  (`kalman_carma_shared_modes-timing`, 200 Hz, commits every 0.1 s, median of
  three runs): A41 2.18, CSM 3.69 us per sample; with x86 flush-to-zero 1.80
  and 2.92 (the float discretization series underflow into denormals, which
  x86 handles in microcode); without commits 1.97 and 3.03.  The complete
  wrapper (`kalman_carma_shared_modes-bench`, mean over three seeds): A41@C
  3.08, CSM 4.61 us, 99th percentile 5.7 and 13.1 us; the worst-case update
  is the commit step.
