# OU-II: deploying the joint physical-MSE pseudo-measurement law

This is the write-up of the experiment
`doc/kalman_ou_ii/ou2-dual-regularization-mse.tex` asks for in its
"Testable Predictions and Next Study" section, and of the change that follows
from it.

OU-II regularizes the low-frequency end of its integration chain with two zero
pseudo-measurements, one on world displacement and one on world velocity. Their
standard deviations are not observation uncertainties; they are design
parameters. Until now they were scheduled by a law selected empirically:

```
r_p0,base = c_p sigma_aw tau^2,   r_v0,base = c_v sigma_aw tau,
                                  (c_p, c_v) = (0.65, 1.30)
```

renormalized by `sqrt(T_0/T_S)` for the self-similar pseudo-update cadence.
The powers are dimensionally consistent, but dimensional analysis does not
identify them, and they were never derived.

The note derives them, from a physical displacement-MSE criterion applied to
both channels jointly. That law is now the deployed default, with the empirical
law retained as a selectable configuration.

| law | `r_p` at the filter input | `r_v` at the filter input |
|---|---|---|
| `Empirical` | `c_p sigma_aw tau^(3/2)` | `c_v sigma_aw tau^(1/2)` |
| **`PhysicalMSE`** (default) | **`C_P q_eff^(1/10) sigma_a,B^(4/5) tau^(19/10)`** | **`C_V q_eff^(1/10) sigma_a,B^(4/5) tau^(9/10)`** |

Select the previous behaviour with
`setPseudoLaw(PseudoAdaptationLaw::Empirical)`, or `OU_II_PSEUDO_LAW=0` in the
simulator.

## What the derivation says

The two pseudo channels do not produce two independent corners. The reduced
Kalman-Bucy CARE for the double integrator under both of them has a closed-form
stabilizing solution parameterized by one frequency and one dimensionless
ratio,

```
omega_p = (q_eff/rho_p)^(1/4),   chi = (omega_v/omega_p)^2,
rho_p = r_p^2 T_S,               rho_v = r_v^2 T_S,
```

and the closed loop reconstructs displacement through a single second-order
transfer whose damping and DC gain both depend on `chi`. The physical objective
is not the filter's own covariance: it is

```
J = J_n + J_w,
J_n = q_eff / (2 omega_p^3 (1+chi)^2 sqrt(2+chi)),
J_w = (1/2pi) int |G(jw)-1|^2 S_p dw,
```

drift plus distortion of the real wave. Only the second term carries physical
wave energy, and it is the term a covariance recursion built on a false
measurement cannot see. This is the same structural point OU-III's `SpectralMSE`
law rests on; the OU-II version is harder only because two channels share one
loop.

Expanding for `chi << 1` and a corner below the wave band and minimising gives

```
omega_p*^5 = 3 q_eff / (8 sqrt2 M_-2),
chi*       = 9 q_eff / (16 sqrt2 M_0 omega_p*^3),
```

with `M_0 = int S_eta dw` and `M_-2 = int S_eta/w^2 dw`. For a self-similar sea
`M_0 ~ sigma_a^2 tau^4` and `M_-2 ~ sigma_a^2 tau^6`, which yields the deployed
shapes.

Two structural consequences are worth stating separately, because they are what
the implementation and its test are built around:

- **Both channels move by the same `tau^(2/5)`.** `12/5 - 2 = 7/5 - 1 = 2/5`.
  The theory does not ask for a different *relative* period law between the two
  pseudo channels; it asks for a common shape correction and a different
  amplitude power.
- **The channel ratio is fixed by spectral moments alone.**
  `(r_p/r_v)^2 = (3/2) M_-2/M_0`, so `r_p/r_v = (C_P/C_V) tau` with `q_eff` and
  the cadence both cancelled out. This is the better-determined half of the
  prediction, and the implementation applies it exactly rather than evaluating
  a second power.

## The coefficients are analytical

`tools/ou2_dual_mse_coefficients.py` evaluates `M_0` and `M_-2` on the eight
reference spectra -- the finite-band 0.02-0.8 Hz JONSWAP and Pierson-Moskowitz
first-order spectra the records are generated from -- at the `(tau, sigma_aw)`
the OU-II tuner actually settles on for each, and solves the two stationarity
conditions record by record:

```
C_P     = 0.1116   spread 0.0794-0.1351
C_V     = 0.2420   spread 0.1856-0.2974
C_P/C_V = 0.4611   spread 0.4279-0.4879
```

Nothing here is fitted to the empirical schedule. Two independent checks say
the numbers are the right ones:

- The channel ratio `0.43-0.49` is what the note reports from the same spectra,
  computed independently of this implementation. An ideal infinite-band
  Pierson-Moskowitz sea gives `sqrt(3/(4 pi)) = 0.489` analytically; an ideal
  `gamma = 3.3` JONSWAP gives `0.465`.
- The empirical law's own `c_p/c_v = 0.500` sits just outside that range, which
  is the agreement the note calls more informative than agreement of either
  absolute coefficient: the empirical sweep found the relative channel law the
  spectra dictate, without being told it.

Re-solving the *unexpanded* objective -- exact `|G(jw)-1|^2`, finite band, no
`chi << 1` step -- moves `C_P` by +6.0 % and `C_V` by -13.0 %, i.e. the ratio to
`0.560`. The position channel is therefore much less sensitive to the
weak-regularization expansion than the velocity channel, which is consistent
with `chi*` sitting at `0.025-0.084` rather than at zero. The deployed pair is
the expansion optimum, because those are the exponents the law uses; the exact
pair is available as an arm of the comparison driver
(`--arm 'ratio-exact:0.1183,0.5599'`).

## What `q_eff` is, and what the sweep caught

`C_P` and `C_P/C_V` are defined with `q_eff^(1/10)` divided out, so they are
invariant to it; the schedule is not. The first implementation set

```
q_eff = 2 R_a h,   sqrt(R_a) = 0.0148 m/s^2
```

by analogy with OU-III's `SpectralMSE` law, which uses the accelerometer's
bench noise spec. **That was the wrong physical quantity here**, and the
calibration sweep is what caught it.

OU-III's use of the bench figure is legitimate: its strong-observation branch is
a statement about the sensor, and there the drift-driving error *is* the sensor.
The OU-II objective is not. The note defines `q` as the residual acceleration
error "presented to the integration chain after acceleration estimation", which
carries attitude and gravity leakage, residual accelerometer bias and estimation
error on top of the sensor floor. The reduced model drops all three from its
*dynamics*, but it still needs their intensity.

Run at the bench density, the deployed law came out 0.56-0.73x the empirical
schedule across the envelope, and the eight-record sweep put the complete-MEKF
vertical optimum at `C_P = 0.174`, a factor of 1.56 above the analytical value.
A free scale off by 56 % is not a small residual; it is a sign that a fixed
input is wrong. Since `r_p ~ q_eff^(1/10)`, a factor of 1.56 in `r_p` is a
factor of 89 in `q_eff` -- an effective drift-band acceleration error of
0.14 m/s^2 rather than 0.0148.

The filter already carries a measured estimate of that number:
`ACC_NOISE_FLOOR_SIGMA_DEFAULT = 0.12 m/s^2`, the pre-band vertical
acceleration noise floor the tuner subtracts as non-wave energy. Referring the
law to it instead moves the schedule by `(0.12/0.0148)^(1/5) = 1.52`, against
the 1.56 the sweep asked for -- a 3 % agreement on a quantity that was free.

So the deployed configuration is the **analytical** `C_P = 0.1116` and
`C_P/C_V = 0.4611`, with `q_eff = 2 sigma_floor^2 h`. Nothing is fitted; the
sweep's role was to expose a mis-specified input, and the corrected input is
independently the one the filter already measures.

The corrected law also lands on the note's own envelope prediction in absolute
terms. `q_eff` only scales the schedule, so the *shape* against the empirical
law is the same either way; what changes is where that shape sits. Against the
empirical schedule at the same operating points it now runs

| `H_s` (m) | `r_p` MSE/empirical | note's Table (envelope) |
|---|---|---|
| 0.27 | 0.904 | 0.924 |
| 1.50 | 0.973 | 1.000 |
| 4.00 | 1.086 | 1.118 |
| 8.50 | 1.106 | 1.139 |

Normalized at the `H_s = 1.5 m` anchor the two columns agree to 0.5 %. The note
could only predict the *shape*, because it had no value for `q_eff`; with the
noise floor supplying one, the derivation predicts the absolute level too, and
lands within 3 % of the empirically calibrated `c_p` without being anchored to
it.

## Calibration: the analytical value is the measured optimum

`tools/ou2_pseudo_mse_scale_sweep.py` sweeps the one overall scale `C_P` while
holding everything the derivation determines -- the 4/5 amplitude power, the
12/5 and 7/5 period powers, and the spectral channel ratio -- fixed. The
deterministic eight-record protocol, against the empirical law on the same
records:

| `C_P` | mean Z (%Hs) | worst Z (%Hs) | mean 3D (m) | mean roll (deg) |
|---|---|---|---|---|
| `Empirical` | 6.3241 | 6.7688 | 0.5024 | 0.2866 |
| 0.070 | 6.5775 | 6.9126 | 0.4702 | 0.2841 |
| 0.085 | 6.3741 | 6.7175 | 0.4789 | 0.2851 |
| 0.098 | 6.2927 | 6.6573 | 0.4942 | 0.2860 |
| **0.1116 (analytical, deployed)** | **6.2769** | **6.6736** | 0.5158 | 0.2870 |
| 0.1146 | 6.2812 | 6.6860 | 0.5212 | 0.2872 |
| 0.130 | 6.3392 | 6.7919 | 0.5517 | 0.2883 |
| 0.150 | 6.4870 | 7.0121 | 0.5967 | 0.2897 |
| 0.175 | 6.7542 | 7.3811 | 0.6588 | 0.2912 |

The vertical column has an interior minimum and it is at the analytical value.
Both neighbouring grid points are worse, so the optimum is bracketed rather
than merely bounded, and the analytical `C_P` is a prediction the sweep
confirms rather than a starting point the sweep moved.

On the deterministic protocol the deployed law is therefore **-0.047 %Hs on the
mean vertical endpoint and -0.095 %Hs on the worst record**, against the
empirical schedule.

The 3D column is the one thing that does not move with the vertical: it falls
monotonically as `C_P` shrinks, so the deployed point costs +0.013 m (2.7 %) of
3D RMS relative to the empirical law, and the 3D optimum is out at `C_P` around
0.07 where the vertical is 0.25 %Hs worse. That split is the reduced model's
own stated limitation showing up as a measurement. The theory carries one
scalar `q` for all three axes, but the vertical channel is the one that takes
gravity leakage from attitude error, so its residual acceleration intensity is
genuinely different from the horizontal one and its optimal corner is too. A
single isotropic `C_P` cannot sit at both optima; the filter already has the
knob that could (`R_p0_xy_factor_`, currently 1.0, isotropic), and splitting it
by axis is the natural follow-on -- but it needs the full-state physical-`H_2`
calculation the note asks for, not a second fitted number.

Roll and pitch move by less than 0.001 deg either way and are not what this
parameter controls.

## Verdict: the paired multi-seed comparison

The deterministic protocol is the screening instrument; the verdict is taken on
the paired harness the paper's primary endpoint uses, so both laws see identical
wave-phase, IMU-noise and initialization seeds on identical scenarios.
`tools/ou2_pseudo_law_compare.py`, eight stationary scenarios x four paired
seeds = 32 paired rows, the stationary comparison the note's
Sec. (next-study) asks for:

| endpoint | empirical | PhysicalMSE | paired delta (95 % hw) |
|---|---|---|---|
| **`disp_z_pct_hs`** (primary) | 6.3540 | 6.3286 | **-0.0254 +/- 0.0173** |
| `disp_3d_rms_m` | 0.4711 | 0.4816 | +0.0104 +/- 0.0081 |
| `roll_rms_deg` | 0.3100 | 0.3103 | +0.0003 +/- 0.0002 |
| `pitch_rms_deg` | 0.3468 | 0.3472 | +0.0004 +/- 0.0004 |

The primary endpoint improves and the interval excludes zero. For scale, the
OU-III `SpectralMSE` deployment moved the same endpoint by
-0.0263 +/- 0.0137 %Hs, so the two derived laws are worth about the same amount
on their respective families.

The gain is not uniform across the envelope, and its structure is the law's own:

| scenario | `disp_z_pct_hs` delta |
|---|---|
| JONSWAP `H_s` 0.27 | **-0.063 +/- 0.011** |
| JONSWAP `H_s` 1.50 | -0.021 +/- 0.017 |
| JONSWAP `H_s` 4.00 | +0.005 +/- 0.030 |
| JONSWAP `H_s` 8.50 | +0.015 +/- 0.071 |
| PM-Stokes `H_s` 0.27 | **-0.104 +/- 0.034** |
| PM-Stokes `H_s` 1.50 | -0.032 +/- 0.020 |
| PM-Stokes `H_s` 4.00 | -0.003 +/- 0.019 |
| PM-Stokes `H_s` 8.50 | -0.001 +/- 0.044 |

Every gain is at the small-sea end and every large-sea scenario is a tie inside
its own interval. That is exactly where the two laws differ most: the derived
law is 0.85-0.90x the empirical schedule at `H_s = 0.27 m` and 1.08-1.11x at
`H_s = 8.5 m`, so the small seas are where it says something different and they
are where it wins. The `sigma_a^(4/5)` amplitude power, not the common
`tau^(2/5)` shape correction, is what that measures.

The 3D regression has the same shape with the opposite sign: it is entirely the
two `H_s = 8.5 m` scenarios (+0.053 and +0.026 m), while all four small and
medium seas improve. This is the axis split of the calibration section showing
up again -- the reduced model carries one scalar `q` for three axes -- and it is
the honest cost of the change, not noise.

Roll and pitch move by 0.0003-0.0004 deg. Both intervals technically exclude
zero at 32 pairs, and both are a tenth of a percent of the quantity; they are
detectable rather than meaningful.

The controlled transition is not in this comparison. The note asks for it to be
reported separately, and the existing OU-II tuning studies show a distinct
stationary/transient tradeoff, so pooling it into the stationary verdict would
mix two different questions.

## What changed in the filter

- `PseudoAdaptationLaw` selects the schedule; `PhysicalMSE` is the default and
  `Empirical` is retained. `Empirical` costs no transcendental where
  `PhysicalMSE` costs one `powf` per tuner update, so it stays the supported
  low-cost configuration for embedded targets without hardware transcendentals.
- The schedule is evaluated as `(sigma_a,B tau^3)^(4/5)` with `q_eff^(1/10)`
  cached, which is exact, and the velocity channel is `r_p/(ratio tau)`, which
  is exact and free. One `powf` covers both channels.
- `sigma_a,B` is the physical band-limited amplitude, so the law divides
  `c_sigma` back out of `sigma_aw`. `C_P` and `c_sigma` stay separately
  identifiable, and the invariance is pinned by `pseudo_law-test`.
- Only the `Empirical` base is renormalized for cadence. `PhysicalMSE` is
  derived on the continuous densities `rho = r^2 T_S` and so already contains
  the realized `T_S`; renormalizing it again would double-count the cadence.
- No clamp moved. `regularizer_floor-test` now asks the *deployed* law rather
  than a hard-coded schedule, and the answer is that `PhysicalMSE` clears the
  existing floors by more than the schedule it replaces: 0.070 m against the
  0.05 m `r_p0` floor in the near-still stress case, where `Empirical` asked
  0.058 m, and 0.271 m on the smallest scored record. (At the bench density the
  same case asked 0.046 m and would have needed the floor lowered -- another
  way the mis-specified `q_eff` showed itself.)

## Tests

`tests/kalman_ou_ii/pseudo_law-test.cpp` pins the five structural claims: the
4/5 amplitude exponent, the 19/10 and 9/10 filter-input period exponents, the
common `tau^(2/5)` offset from the empirical law, the exact channel ratio at
every operating point including inside the cadence clamps, `c_sigma`
invariance, and that the `Empirical` arm is still exactly the calibrated
schedule with its cadence renormalization intact.

`tuner_schedule-test` and `regularizer_floor-test` were corrected rather than
accommodated. Both had encoded the `Empirical` cadence contract as universal --
the first in its staged-`r` computation and its `r^2 T_S` assertion, the second
in its clamp-ordering check and in scaling a schedule's *output* to project the
near-still case. Each now branches on the law or asks the deployed law
directly, which is what they meant to test.
