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
  (`kalman_carma_shared_modes-math.tex`, the model;
  `kalman_carma_shared_modes-minimal.tex`, the minimal realization and the
  numerical study; `kalman_carma_shared_modes-charts.tex`, the replay chart
  gallery).

**Findings** (details in section 6; every number from the committed study data):

@@FINDINGS@@

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

Here `gamma^2` is the debiased multiple coherence of z on (x, y), and `b` is
scaled down if `b^H H^-1 b` would exceed `0.999 sigma_z^2`.  The excitation is
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
is kept.  The cost is one 9 x 9 LDLT and a 12 x 9 x 6 product, at most at
10 Hz and only while `c` falls.

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

@@VALIDATION@@

## 7. Embedded efficiency

@@EMBEDDED@@
