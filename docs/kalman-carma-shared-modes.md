# Shared-modes QMEKF (`KalmanCarmaSharedModes`, experimental)

**Status: experimental, simulation only, not a default.**  The shipping
filter stays `SeaStateFusionFilter_OU_III`; the independent-axis CARMA3 QMEKF
(`src/kalman_carma/`, benchmark name `A41`) is unchanged and remains the
baseline.  Nothing here touches TFG or TFG-CARMA3.

This filter promotes Option C of the coupled-prior study
([carma-coupled-three-axis.md](carma-coupled-three-axis.md)) -- shared
directional stochastic modes -- into a standalone QMEKF.  It adds the
coupling confidence and fallback that the study showed to be necessary.

* Core: `src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h`,
  `KalmanCarmaSharedModes<T = float>`, the QMEKF (scalar `T`, as every
  Kalman core in the repository), with `CarmaCouplingConfidence.h` and
  `CarmaSharedModesMotion.h`.
* Wrapper: `SeaStateFusionFilter_CarmaSharedModes.h`,
  `SeaStateFusion_CarmaSharedModes`, the unchanged CARMA3 schedule with the
  core as its inner filter (float, as every wrapper).
* Tests and simulators: `tests/kalman_carma_shared_modes/`
  (`carma_shared_modes-test`; `kalman_carma_shared_modes-sim`, the dataset
  replay; `kalman_carma_shared_modes-bench`, the paired Monte Carlo).
* Validation: `tools/carma_shared_modes_validation.py` ->
  `reports/results/kalman_carma_shared_modes/`; deterministic table
  `tools/carma_shared_modes_sim_table.py`.
* Articles: `doc/kalman_carma_shared_modes/`
  (`kalman_carma_shared_modes-math.tex`, `kalman_carma_shared_modes-validation.tex`,
  and the replay chart gallery `kalman_carma_shared_modes-charts.tex`).

**Result (simulation only).**  Completed Monte Carlo: 18 seas x 50 seeds,
robustness 4 seas x 30 seeds, plus the dataset replay of the eight pinned
records.  Generated tables: `reports/results/kalman_carma_shared_modes/summary.md`.

* **Shared modes alone** (CSM against the CARMA3 baseline at the same
  Evidence excitation, A41@C): pooled over the 16 wave seas the displacement
  RMS is lower by **4.1 % (x), 1.8 % (y), 4.5 % (z)**, 95 % intervals within
  +-0.3 points.  CSM is significantly better on 15, 15 and 13 of the 16 seas.
  It is worse only on manoeuvring (-6.8 % heave), short-period wind waves
  (-0.8 % heave) and the 5 cm sea (-0.3 % x).  On the dataset replay its 3-D error is below the coupling-off
  variant's on all eight records (by 2-10 %).
* **Against the historical A41** (Legacy excitation): +9.3 / +5.6 / +14.1 %.
  Of that, Evidence excitation alone gives +5.4 / +3.8 / +10.0 %: it removes
  the calm-water phantom displacement (52-61 mm -> 0) but costs 25 %
  horizontally on a 5 cm sea.
* **Coupling confidence.**  It rejects every wrong prior injected (reversed,
  turned 90 deg, phase +45 deg, doubled amplitude): those runs equal the
  independent-axis fallback, where without the confidence they lose up to
  131 % with NEES up to 134.  It is conservative: it also withdraws a correct
  coupling on two-system seas and heading zigzags, so CSM keeps less gain than
  the research filter (SHM:est +6.0 / +2.6 / +6.0 % against A41 at the same
  Legacy excitation).  In exchange it halves the manoeuvring loss (-6.3
  against -12.7 % heave vs A41) and turns the wave-onset loss into a gain
  (+3.2 % vs A41@C, where SHM:est loses 3.5 % vs A41).
* **Consistency.**  CSM never failed (A41 tripped the bench's divergence check on seed 44 in
  calm water, still IMU and onset).  Its mean NEES is at or below A41@C's on
  every sea except wave onset (4.8-6.5 on the irregular single-system seas).
  In Evidence mode, though, the motion covariance collapses before waves
  arrive, so NEES is about 10^4 on wave onset for A41@C and CSM alike
  (section 8).
* **Cost.**  Host 6.6 us mean (32.9 us p99) per update against 3.2 (7.0) us.
  ESP32-S3: code 196 KB against 173 KB, QMEKF object 15.8 KB against
  6.1 KB, update-path stack 4.9 KB against 4.8 KB.
* **Status:** experimental, not a default; the CARMA3 baseline is unchanged.

## 1. What is reused, what is new

| part | source | status |
|---|---|---|
| QMEKF core: attitude error, gyro and accelerometer bias, Joseph updates, multiplicative reset, magnetometer | `kalman_carma/Kalman3D_Wave_CARMA3.h` (`Motion` parameter) | reused unchanged |
| shared-mode motion block (complexified CARMA(4,1) per axis, Hermitian excitation) | `kalman_carma/CarmaCoupledMotion.h` (`CarmaSharedModeMotion`) | reused unchanged (derived class) |
| coupling estimate and excitation mapping | `kalman_carma/CarmaCrossSpectrumEstimator.h` | reused unchanged |
| schedule: proxy startup, magnetic acquisition, tuner, 0.1 s commits | `kalman_carma/SeaStateFusionFilter_CARMA3.h` | reused unchanged |
| short-horizon monitor of the coupling | `CarmaSharedModesMotion.h` | new |
| coupling confidence, hysteresis, smoothing | `CarmaCouplingConfidence.h` | new |
| PSD-preserving decoupling of the covariance | `Kalman3D_Wave_CARMA3::shrink_motion_vertical_coupling` | new member, instantiated only for a motion block |
| QMEKF core `KalmanCarmaSharedModes<T>`: confidence in `time_update` / `measurement_update_acc_only`, diagnostics, direction estimate | `KalmanCarmaSharedModes.h` | new (derives from the CARMA3 core) |
| float wrapper `SeaStateFusion_CarmaSharedModes` | `SeaStateFusionFilter_CarmaSharedModes.h` | new (the CARMA3 wrapper with the new core) |

With the confidence disabled, the filter is the research filter `SHM:est`.
`carma_shared_modes-test` [2] finds them bit-identical over 400 s.  In the
benchmark they agree to float rounding (5e-5 relative): the two objects
differ in memory layout, which can change the summation order of vectorized
loops.  With the confidence held at zero it is `CouplingMode::Off`, bit for
bit, which is the CARMA(4,1) law on every axis.

## 2. State and model

Error state (33): `[dtheta(3), b_g(3), x(8), y(8), z(8), b_a(3)]`, NED,
attitude error in the body frame.  Per axis the complexified CARMA(4,1) mode
of the C41J shape is

    zeta = B(s) / (P+(s) D(s)) xi,   real state [Re(p, v, a, l), Im(p, v, a, l)]

* `P+` holds the upper-half-plane poles of the C41J polynomial `P`.
* `D` is the real Hurwitz quadratic with
  `|D(iw)|^2 = (|P+(iw)|^2 + |P+(-iw)|^2) / 2`.
* `xi` is circular complex white noise.
* `Re zeta` has exactly the C41J law, and the noise enters `a` only
  (relative degree 3).
* `p, v, a = Re(...)` are states, so every accelerometer row is a unit
  selector rotated by `R_wb`.

The axes share the dynamics, `F = I_3 (x) F_8`.  The axes are driven with a
Hermitian PSD excitation `Sigma`:

    Q_ab = Re(Sigma_ab) Q_8 + Im(Sigma_ab) J Q_8,   J = [[0, -I], [I, 0]]

The displacement cross-spectrum is
`S_ij(w) = (g^2/2) [Sigma_ij |H(iw)|^2 + conj(Sigma_ij) |H(-iw)|^2]`.  Its
co-spectrum is `Re Sigma` and its quadrature spectrum `Im Sigma`, which is the
90 deg x-z relation of linear waves, shifted by the heave RAO.  Derivations,
exact discretization, the complex-OU factorization (Fasen-Hartmann and Scholz)
and the MCARMA context (Stelzer's 2024 correction) are in the mathematical
article.

**Excitation** (deployable).  Each sample, the filter complex-demodulates its
own world-frame acceleration `R_wb^T (f - b_a) + g e_z` at the operating
frequency `f_op`, then low-passes the envelope at `0.15 f_op` and averages
`A A^H` over 240 s.  The result is the band cross-spectrum `S_hat`.  The
excitation keeps the baseline diagonal and adds only the horizontal-vertical
terms:

    b = s gain c S_hat_hz,   c = sigma_z^2 / S_hat_zz,   s = max(0, 1 - 1/(N_eff gamma^2))

Here `gamma^2` is the debiased multiple coherence of z on (x, y), and `b` is
scaled down if `b^H H^-1 b` would exceed `0.999 sigma_z^2`.  The excitation is
PSD by construction, equivariant under rotations about the vertical, and the
baseline diagonal without coupling evidence.  `gain` is the coupling
confidence `c in [0, 1]` of section 4.

**Sea-state adaptation.**  This is the CARMA3 wrapper's own adaptation,
unchanged.  Every 0.1 s it commits `kappa = f_tune / f_z,base` (one frequency
scale for every axis) and the isotropic amplitudes; a commit rebuilds `F`, `Q`
and the excitation.  The posterior covariance is reset to the stationary
prior only on entering Live, never during adaptation.

**Adaptation and the complexified realization.**  At fixed parameters the
complexified mode is the CARMA(4,1) law exactly.  Under adaptation it is not,
because matching spectra no longer imply matching laws once `sigma(t)` and
`kappa(t)` vary.  With coupling off (`CSM_OFF`) and the tuner frozen after
300 s, the dataset replay matches A41 to 0.03 %.  With adaptation it differs
by -1.6 % to +6 % in heave, and is worst on the two Hs 0.27 m records.  A
kappa commit deadband changes nothing (<= 0.2 %).  An amplitude deadband
narrows that gap on the dataset records, but on the screening seeds
1001-1020 it costs heave on most seas and up to 6.9 % horizontal on the very
low sea, so the pre-registered criterion rejected it (`screen_runs.csv`).
The default commits every amplitude, as the research filter does.

**Wave excitation: Evidence mode.**  The amplitude jumps come from the Legacy
excitation chain, which applies the 0.05 m/s^2 floor and the stillness
attenuation at 0.28 m/s^2 vertical RMS (docs/carma-wave-excitation.md).  On
the Hs 0.27 m records the sea sits near that threshold.  Evidence mode, the
still-sea phantom fix that the AtomS3R sketches deploy, replaces both with
the noise-debiased band variance and a detection weight.  It is therefore the
shared-modes filter's default (`SeaStateFusion_CarmaSharedModes::begin`,
`kalman_carma_shared_modes-sim`).  On the dataset replay with Evidence mode:

| record | heave Z % Hs, A41 Legacy / A41 Evidence / CSM | 3D %, A41 Evidence / CSM |
|---|---|---|
| JONSWAP 0.27 m | 3.06 / 2.46 / 2.52 | 4.06 / 4.03 |
| PM-Stokes 0.27 m | 3.21 / 2.41 / 2.46 | 4.08 / 4.07 |
| JONSWAP 1.5 m | 2.25 / 2.21 / 2.19 | 4.48 / 4.34 |
| PM-Stokes 8.5 m | 2.37 / 2.36 / 2.25 | 5.27 / 4.85 |

CSM's regression bars are fitted in this mode.  They are below the CARMA3
and TFG-CARMA3 bars (Legacy) for heave (2.531 / 2.470 against 3.072 / 3.231
% Hs), 3D (4.419 / 4.873 against 4.890 / 5.278) and yaw (0.755 against
0.766 deg).  They are above those bars for pitch (0.225 against 0.214 deg)
and for the accelerometer bias (84.2 against 82.4 / 82.5 %).  On the
gyroscope bias, CSM's bar (12.71 %) is below CARMA3's (12.78 %) and above
TFG-CARMA3's (12.40 %): the two-frame Lie group separates the gyroscope bias
better, and the QMEKF comparison is the relevant one here.

**Bias estimation against the CARMA3 QMEKF** (dataset replay, Evidence
mode, mean over the eight records):

| quantity | A41 | CSM_OFF | CSM |
|---|---|---|---|
| gyroscope bias, 3D RMS (mdeg/s) | 6.774 | -- | 6.765 |
| accelerometer bias, vertical RMS (mm/s^2) | 1.867 | -- | 1.862 |
| accelerometer bias, horizontal RMS (mm/s^2) | 49.20 | 49.19 | 50.43 |

The gyroscope and vertical accelerometer biases are equal to or slightly
better than CARMA3's on every record.  The horizontal accelerometer bias is
2-4 % worse on the seas of Hs 1.5 m and above and equal on the small seas.
That cost comes from the coupling, not the realization (`CSM_OFF` matches
A41), and is the price of its 7-24 % lower surge error.  Like for like, against A41 in Evidence mode, CSM's
small-sea heave is about 2 % worse, a residue of the realization effect
above, while its 3D error is lower on every record.  The Monte Carlo pairs
CSM with both A41 (Legacy, the historical baseline) and A41@C (Evidence).

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

**Applying c.**  At the wrapper's 0.1 s commit the excitation is rebuilt
with `gain = c`, so `c = 0` is the independent-axis law exactly (test [2]).
The state is never touched, so position, velocity and attitude are
continuous through every change (test [6]).

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
is kept.  The cost is one 9 x 9 LDLT and 2 kflop, at most at 10 Hz and only
while `c` falls.

**Recovery.**  Recovery is automatic.  Once the evidence returns, `c` rises
with `tau = 60 s` and the cross-covariance rebuilds through the dynamics and
the measurements (test [6]: a reversed prior corrected at 450 s returns to
`c > 0.99`).

**What the confidence cannot do.**  It cannot detect a wrong coupling that
the IMU data cannot contradict.  The consistency test compares the prior
with the same IMU's recent band statistics: a heading error, for example,
rotates both the prior and the monitor, so it shows up as heading error and
not as a coupling mismatch.  It also reacts only after the evidence has
accumulated (seconds to a minute).  The robustness results in section 6
quantify both limitations.

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

`carma_shared_modes-test` (43 checks, run by `make all`):

| # | group | result |
|---|---|---|
| 1 | 33 states, bit-identical repeat runs and reset, finite symmetric covariance | yes |
| 2 | confidence disabled = research SHM:est (400 s); confidence 0 = Off | bit-identical; bit-identical |
| 3 | fallback algebra: PSD, only P_hz changes, state untouched, alpha = 1 no-op, alpha = 0 leaves E | min corr. eig 0.015; yes; yes; yes; 6e-7 |
| 4 | confidence: consistent -> 1, reversed -> fallback < 0.05 in 10 s, rate bounds, hysteresis, recovery, manoeuvre hold, sea change, low energy, NIS ratio, warm-up | yes |
| 5 | direction of travel and phase from a regression | 4e-6 deg |
| 6 | synthetic directional sea: estimated coupling (c 0.99), reversed external prior with / without confidence, continuity through a forced fallback, recovery | see below |

Test 6 runs on an exact-kinematics directional sea, scored after 300 s.

| run | horizontal RMS | heave RMS | NEES_p |
|---|---|---|---|
| independent axes | 23.1 mm | 15.4 mm | 5.6 |
| estimated coupling | 19.8 mm | 13.2 mm | 6.1 |
| reversed external prior, no confidence | 119.9 mm | 42.9 mm | 92.8 |
| reversed external prior, with confidence | 23.1 mm | 15.4 mm | 5.6 |

With the confidence, the reversed prior's run is the independent-axis result:
the coupling was never applied.

`kalman_carma_shared_modes-sim` replays the eight scored records of the pinned
oceanography-waves-lib v1.2.3 dataset with regression bars (the CARMA3 replay
adapter with the shared-modes core); `kalman_carma_shared_modes-bench` runs
carma3_bench's Monte Carlo engine with this filter registered next to A41 and
SHM:est.

## 7. Results (completed Monte Carlo)

Every row runs one scenario and seed fed identically to every filter: A41
(the CARMA3 baseline, Legacy excitation), A41@C (the same with Evidence
excitation), SHM:est (the research filter, Legacy), CSM and CSM:off (both
Evidence).  Errors from 300 s; changes are `1 - mean RMS(candidate) / mean
RMS(base)` with 95 % paired-bootstrap intervals over seeds (positive =
better; `*` = interval includes zero).  Pooled figures are the mean over the
16 wave seas of the log RMS ratio, seeds resampled jointly.

### 7.1 Final comparison (50 seeds)

| sea | A41@C RMS x / y / z (mm) | CSM vs A41@C x / y / z (%) | CSM vs A41 x / y / z (%) | SHM:est vs A41 z (%) | NEES mean A41@C / CSM |
|---|---|---|---|---|---|
| calm water (afloat, no waves) | 0.07 / 0.06 / 0.06 um | -- (A41@C, CSM and CSM:off all below 0.1 um) | A41 61.2 / 60.8 / 52.2 mm -> 0 | -0.0* | 0.0 / 0.0 |
| still IMU | 0.07 / 0.06 / 0.06 um | -- (as calm water) | A41 59.1 / 58.3 / 52.2 mm -> 0 | -0.0* | 0.0 / 0.0 |
| regular monochromatic wave | 31.4 / 28.4 / 23.9 | +11.4 / +2.4 / +10.2 | +11.4 / +2.2 / +10.1 | +10.1 | 13.1 / 12.8 |
| unidirectional JONSWAP | 15.7 / 14.8 / 14.5 | +4.0 / +0.9 / +5.3 | +3.5 / +0.1* / +5.7 | +5.2 | 6.0 / 5.8 |
| directionally spread JONSWAP (cos^8) | 15.4 / 14.9 / 14.6 | +3.4 / +1.1 / +4.5 | +2.9 / +0.3* / +4.9 | +4.4 | 5.6 / 5.5 |
| Pierson-Moskowitz | 21.4 / 19.8 / 20.9 | +5.6 / +2.4 / +7.5 | +5.4 / +2.0 / +8.0 | +8.0 | 6.0 / 5.7 |
| long-period swell | 82.6 / 66.4 / 65.3 | +19.0 / +9.5 / +18.4 | +28.6 / +2.4 / +40.4 | +23.7 | 7.4 / 6.5 |
| short-period wind waves | 5.7 / 5.8 / 6.6 | +1.6 / +1.2 / -0.8 | +0.3 / +0.1* / -0.2* | -1.2 | 3.8 / 3.7 |
| mixed wind sea + swell | 184.0 / 50.6 / 146.4 | +1.3 / +1.1 / +3.1 | +2.0 / +1.0 / +4.0 | +11.9 | 59.7 / 57.9 |
| two systems, different frequencies | 66.9 / 40.0 / 52.3 | +3.8 / +0.2 / +5.4 | +4.3 / -1.4 / +6.5 | +16.1 | 9.5 / 8.9 |
| two crossing systems (90 deg) | 17.2 / 17.1 / 16.8 | +1.7 / +0.8 / +1.9 | +1.1 / +0.3 / +2.5 | +1.7 | 5.7 / 5.6 |
| changing wave direction (30 -> 120 deg) | 14.9 / 15.8 / 14.4 | +1.1 / +2.4 / +4.0 | +0.4 / +2.0 / +4.5 | +3.9 | 5.5 / 5.4 |
| abrupt change of wave period | 37.1 / 43.7 / 36.2 | +2.5 / +4.1 / +7.8 | +2.3 / +4.4 / +8.5 | +11.6 | 9.7 / 9.5 |
| wave onset and cessation | 11.6 / 11.5 / 11.2 | +2.6 / +1.0 / +3.2 | +65.4 / +64.5 / +71.4 | -3.5 | 9682 / 10216 |
| vessel heading changes (90 deg zigzag) | 15.2 / 14.9 / 14.5 | +1.4 / +0.7 / +1.9 | +1.0 / -0.2* / +2.3 | +4.5 | 4.9 / 4.8 |
| vessel manoeuvring | 78.4 / 79.7 / 14.5 | +0.9 / +0.5 / -6.8 | -0.5 / -1.1 / -6.3 | -12.7 | 59.5 / 59.3 |
| very low sea state (5 cm) | 3.9 / 4.0 / 4.4 | -0.3 / -0.1* / -0.0* | -25.3 / -26.3 / +7.7 | +3.3 | 3.5 / 3.5 |
| magnetic disturbances, calibration residuals | 15.7 / 15.4 / 14.7 | +3.9 / +0.5 / +4.3 | +3.3 / -0.4* / +4.7 | +4.2 | 6.3 / 6.3 |

| pooled over the 16 wave seas | x | y | z |
|---|---|---|---|
| CSM vs A41@C (shared modes alone) | +4.1 % [+3.9, +4.3] | +1.8 % [+1.5, +2.1] | +4.5 % [+4.4, +4.7] |
| CSM vs A41 | +9.3 % [+8.9, +9.6] | +5.6 % [+5.1, +6.1] | +14.1 % [+13.7, +14.5] |
| A41@C vs A41 (Evidence excitation alone) | +5.4 % [+5.0, +5.8] | +3.8 % [+3.4, +4.2] | +10.0 % [+9.6, +10.4] |
| SHM:est vs A41 (research filter, Legacy) | +6.0 % [+5.8, +6.2] | +2.6 % [+2.1, +3.0] | +6.0 % [+5.9, +6.2] |
| CSM:off vs A41@C (realization, no coupling) | -0.1 % [-0.2, -0.1] | -0.0 % [-0.1, +0.0] | -0.7 % [-0.8, -0.6] |

* **Where the confidence gives gain away.**  Mean applied x-z coherence,
  CSM / SHM:est: two systems at different frequencies 0.09 / 0.47, mixed sea
  (y-z) 0.06 / 0.38, heading zigzag 0.03 / 0.13.  The single shared mode
  cannot represent two wave systems, so the short-horizon monitor disagrees
  with the prior and the consistency test withdraws the coupling.  The
  yaw-rate hold withdraws it in turns, although a world-frame coupling stays
  valid through heading changes (coupled-prior study, section 7).
* **Where it removes a loss.**  Manoeuvring (speed and heading changes,
  Doppler-shifted encounter frequency): -6.8 % heave vs A41@C, against
  -12.7 % for SHM:est vs A41.  The remaining loss is the coupling's
  (CSM:off -1.2 %).  Wave onset: +3.2 % vs A41@C, against -3.5 % for SHM:est.
* **Evidence excitation** (A41@C vs A41) dominates calm water, onset and
  cessation (+64-70 %) and swell heave (+27 %).  It loses 25-26 %
  horizontally on the 5 cm sea (3.1 -> 3.9 mm), where it gains 8 % in heave.

### 7.2 Coupling-mismatch robustness (30 seeds)

External priors built from the true wave acceleration and then deliberately
falsified enter the filter, which must reject them from its own data.
Changes against A41; h = horizontal RMS.

| injected prior | sea | without confidence h / z (%) | NEES | with confidence h / z (%) | NEES |
|---|---|---|---|---|---|
| true coupling | directionally spread JONSWAP (cos^8) | +1.8 / +4.9 | 5.7 | +1.8 / +4.9 | 5.7 |
| true coupling | mixed wind sea + swell | +8.5 / +13.8 | 58.5 | +2.0 / +4.1 | 59.4 |
| true coupling | Pierson-Moskowitz | +4.4 / +8.5 | 6.1 | +4.1 / +8.1 | 6.1 |
| reversed (180 deg) | directionally spread JONSWAP (cos^8) | -8.1 / -27.7 | 6.9 | -0.5 / -0.7 | 5.9 |
| reversed (180 deg) | mixed wind sea + swell | -21.7 / -18.2 | 122.5 | -0.3 / +0.9 | 62.3 |
| reversed (180 deg) | Pierson-Moskowitz | -12.9 / -29.7 | 8.0 | -0.3 / +0.7 | 6.3 |
| turned 90 deg | directionally spread JONSWAP (cos^8) | -3.8 / -10.7 | 6.9 | -0.5 / -0.7 | 5.9 |
| turned 90 deg | mixed wind sea + swell | -25.4 / -21.0 | 133.6 | -0.3 / +0.9 | 62.3 |
| turned 90 deg | Pierson-Moskowitz | -4.4 / -10.5 | 7.7 | -0.3 / +0.7 | 6.3 |
| phase +45 deg | directionally spread JONSWAP (cos^8) | +3.1 / +0.3 | 6.1 | -0.5 / -0.7 | 5.9 |
| phase +45 deg | mixed wind sea + swell | -8.8 / +9.7 | 86.8 | -0.3 / +0.9 | 62.3 |
| phase +45 deg | Pierson-Moskowitz | +2.9 / +2.7 | 6.7 | -0.3 / +0.7 | 6.3 |
| doubled amplitude | directionally spread JONSWAP (cos^8) | -118.5 / -131.2 | 36.3 | -0.5 / -0.7 | 5.9 |
| doubled amplitude | mixed wind sea + swell | +0.1* / +37.5 | 109.4 | -0.3 / +0.9 | 62.3 |
| doubled amplitude | Pierson-Moskowitz | -91.4 / -89.6 | 32.6 | -0.3 / +0.7 | 6.3 |
| none (coupling off, CSM:off) | directionally spread JONSWAP (cos^8) | -- | -- | -0.5 / -0.7 | 5.9 |
| none (coupling off, CSM:off) | mixed wind sea + swell | -- | -- | -0.3 / +0.9 | 62.3 |
| none (coupling off, CSM:off) | Pierson-Moskowitz | -- | -- | -0.3 / +0.7 | 6.3 |
| estimated coupling | vessel heading changes (90 deg zigzag) | +1.8 / +5.1 | 4.8 | +0.4 / +2.4 | 4.9 |
| frozen at 400 s in the keel frame | vessel heading changes (90 deg zigzag) | -0.1* / -1.1 | 5.1 | +0.0* / +1.1 | 4.9 |

The estimated coupling (CSM, final rows of the same seeds) gives within 0.3
points of the true-coupling prior on all three seas (spread JONSWAP +1.7 /
+5.0 against +1.8 / +4.9).  With the confidence, every falsified prior gives
exactly the coupling-off result: the coupling was never applied.  Without it,
the reversed and turned priors lose 4-30 % and the doubled amplitude up to
131 %, with NEES up to 134.  The cost of that protection is the mixed sea: the
confidence also withdraws the true coupling there (+2.0 / +4.1 against
+8.5 / +13.8 %).

### 7.3 Spectra and coupling estimate

Energy-weighted relative error of the displacement cross-spectral matrix
against the truth's, 0.05-1.5 Hz, 10 seeds: JONSWAP A41 0.084, SHM:est 0.074,
CSM 0.073; crossing seas 0.082 / 0.075 / 0.074; mixed sea 0.758 / 0.620 /
0.741.  Every filter reproduces the truth's x-z phase at the spectral peak.
The estimated excitation phase is stable across seeds: x-z 124 +- 3 deg on
JONSWAP, 122 +- 3 deg on PM, 96 +- 7 deg on swell.

### 7.4 Consistency

Mean NEES of `[p, v]` (6 if consistent): CSM is at or below A41@C on every
sea except wave onset.  It is 4.8-6.5 on the irregular single-system seas,
3.5-3.7 on the very low sea and wind waves, and 0 in calm water (no error).
Every model is overconfident on the regular wave (CSM 12.8), the period
change (9.5), the two-peak seas (two systems 9.0, mixed 57.9: one CARMA mode
cannot represent two spectral peaks) and manoeuvring (59.3, where the
horizontal 95 % coverage is 0.81).  Wave onset is discussed in section 8;
its position-interval coverage stays at 0.96-0.97.

### 7.5 Embedded cost

| | A41 | CSM |
|---|---|---|
| host update, mean / p99 (us) | 3.2 / 7.0 | 6.6 / 32.9 |
| ESP32-S3 code / rodata (B) | 173,375 / 125,633 | 196,266 / 133,364 |
| QMEKF object (B) | 6,136 | 15,776 |
| largest stack frame / update-path stack (B) | 1,328 / 4,832 | 2,480 / 4,896 |

Xtensa GCC 14.2.0, -Os (`tests/kalman_carma_shared_modes/arduino/csm_size_probe`).
The confidence adds 0.2 us per update over SHM:est.

## 8. Limitations

* **Simulation only.**  No hardware or field evidence; the Monte Carlo hull,
  RAOs and sensor model are the CARMA3 benchmark's.
* **Evidence-mode covariance at wave onset.**  With no wave evidence the
  Evidence excitation is zero and the motion covariance collapses (sigma_z
  0.000 mm before onset).  When waves start, the filter carries a 10 mm heave
  error against that covariance: NEES about 10^4 for A41@C, CSM and CSM:off
  alike, while Legacy A41 stays at 5.7.  The accuracy is still better than
  Legacy's (10.5 against 100 mm over 300-400 s).  This belongs to the
  excitation mode the AtomS3R sketches deploy, not to the coupling; it is not
  changed here.
* **Conservative confidence.**  The consistency test and yaw-rate hold give
  away 6-11 points on two-system seas and about 3 on heading changes.
* **Manoeuvring.**  CSM still loses 6.8 % heave while manoeuvring.
* **Very low seas.**  Evidence excitation costs 25 % horizontally on a 5 cm
  sea (0.8 mm).
