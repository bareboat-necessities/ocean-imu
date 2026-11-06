# Shared-modes QMEKF (`KalmanCarmaSharedModes`, experimental)

**Status: experimental, simulation only, not a default.**  The shipping
filter stays `SeaStateFusionFilter_OU_III`; the independent-axis CARMA3 QMEKF
(`src/kalman_carma/`, benchmark name `A41`) is unchanged and remains the
baseline.  Nothing here touches TFG or TFG-CARMA3.

This filter promotes Option C of the coupled-prior study
([carma-coupled-three-axis.md](carma-coupled-three-axis.md)) -- shared
directional stochastic modes -- into a standalone QMEKF.  It adds the
coupling confidence and fallback that the study showed to be necessary.

* Header: `src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h`
  (`CarmaCouplingConfidence.h`, `CarmaSharedModesMotion.h`).
* Tests and simulator: `tests/kalman_carma_shared_modes/`
  (`carma_shared_modes-test`, `kalman_carma_shared_modes-sim`).
* Validation: `tools/carma_shared_modes_validation.py` ->
  `reports/results/kalman_carma_shared_modes/`.
* Articles: `doc/kalman_carma_shared_modes/`
  (`kalman_carma_shared_modes-math.tex`, `kalman_carma_shared_modes-validation.tex`).

RESULTS_SUMMARY_PLACEHOLDER

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
| public API, diagnostics, direction estimate | `KalmanCarmaSharedModes.h` | new |

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

```cpp
#include "kalman_carma_shared_modes/KalmanCarmaSharedModes.h"
using ocean_imu::carma_shared_modes::KalmanCarmaSharedModes;

KalmanCarmaSharedModes::Config cfg;                // sensors in cfg.fusion (as SeaStateFusion_CARMA3)
cfg.fusion.sigma_a = ...; cfg.fusion.sigma_g = ...; cfg.fusion.sigma_m = ...;
KalmanCarmaSharedModes f(cfg);                     // or f.begin(cfg); f.reset() restarts
f.update(dt, gyr_body_ned, acc_body_ned, temp_c);  // every IMU sample
f.updateMag(mag_body_ned_uT);                      // at the magnetometer rate
f.isLive(); f.attitude(); f.position(); f.velocity(); f.waveAcceleration();
f.accelBias(); f.gyroBias(); f.positionVelocityCovariance();
f.seaState();        // frequency, operating frequency, kappa, committed accel std
f.waveDirection();   // direction of travel (deg from magnetic north), uncertainty, phase
f.coupling();        // mode, confidence c, state (warm-up/coupled/reduced/fallback),
                     // D^2, SNR, NIS ratio, applied coherence / phase
f.health();          // NIS, covariance finite / symmetric / min diagonal
```

`CouplingMode::Estimated` is the default and the only deployable mode.
`External` (with `setExternalCoupling`) is for research priors and for the
robustness tests.  `Off` gives independent axes.  `forceCouplingConfidence`
is a test hook.  Ground truth and RAO knowledge never enter the class.

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
oceanography-waves-lib v1.2.3 dataset with regression bars, and runs
carma3_bench's Monte Carlo engine with this filter registered next to A41 and
SHM:est.

RESULTS_PLACEHOLDER
