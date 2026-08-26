# Why the vertical %H_s error is worse on the small seas

All three deployed families report a higher vertical RMS error *as a fraction of
the sea* on the smallest records than on the largest. Vertical RMS in percent of
`H_s` on the JONSWAP records, adaptive mode: OU-II and OU-III from the committed
ten-seed bundle (`reports/results/ou_validation/`, mean +/- sample sd), TFG from
its committed single-realization table:

| record | OU-II | OU-III | TFG |
| --- | ---: | ---: | ---: |
| H_s = 0.27 m | 6.86 +/- 0.23 | 4.96 +/- 0.14 | 4.67 |
| H_s = 1.50 m | 6.38 +/- 0.13 | 4.80 +/- 0.16 | 4.37 |
| H_s = 4.00 m | 6.39 +/- 0.33 | 4.84 +/- 0.19 | 4.26 |
| H_s = 8.50 m | 6.02 +/- 0.34 | 4.02 +/- 0.17 | 3.81 |

The gradient is the same sign in every family: about 14 % more relative error on
the smallest sea than on the largest for OU-II, 23 % for OU-III and 23 % for
TFG, over a 31x range of `H_s`. It is not monotone record by record -- OU-II and
OU-III both score the 4.00 m record marginally worse than the 1.50 m one -- but
the endpoints are separated by several times the seed spread.

The natural suspect is the `r_S` schedule, since it is the only sea-state-
dependent quantity between the tuner and the linear block, and a coefficient
fitted on the middle of the envelope would show up exactly this way.

It is not the `r_S` coefficient. This note is the measurement.
`tools/ou_low_sea_error_study.py` reproduces all three sections. The horizontal
side of the same question -- why OU-III raises 3D RMS while lowering vertical
RMS -- is [`ou-3d-error-attribution.md`](ou-3d-error-attribution.md).

## 1. C_J is already at its per-record optimum on every sea

`C_J` multiplies the whole SpectralMSE schedule, so sweeping it moves `r_S` by a
common factor on every record. Sweeping it one record at a time locates each
sea's own optimum. The deployed floor `MIN_R_S = 0.15` binds on the smallest
records -- the schedule asks for 0.16 to 0.19 there -- so the sweep is run with
the floor released to 0.001, otherwise the bottom of the grid measures the clamp
rather than the schedule.

Vertical RMS as a percentage of the reference RMS, trailing 900 s, default
seeds, `OU_III_R_S_MIN=0.001`:

| record | 0.0135 | 0.027 | **0.0538** | 0.108 | 0.215 | argmin | applied `r_S` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| JONSWAP H_s = 0.27 m | 30.403 | 21.301 | **17.877** | 20.736 | 29.570 | 0.0538 | 0.18 |
| JONSWAP H_s = 1.50 m | 25.732 | 18.800 | **16.786** | 20.367 | 29.123 | 0.0538 | 1.46 |
| JONSWAP H_s = 4.00 m | 23.914 | 17.716 | **15.741** | 18.089 | 24.573 | 0.0538 | 10.54 |
| JONSWAP H_s = 8.50 m | 23.734 | 17.406 | **14.686** | 15.970 | 23.064 | 0.0538 | 18.25 |

Three things follow.

- **Every record's argmin is the deployed 0.0538.** The optimum does not drift
  with sea state over a 16x sweep of the coefficient and a 31x range of `H_s`.
- **The minima are steep, not flat.** One grid step either way costs 8 to 20 %,
  so this is a resolved optimum on the small records too, not a plateau the
  sweep cannot see into.
- **The gradient survives per-record optimal tuning.** Reading down the bold
  column: 17.877, 16.786, 15.741, 14.686. Each entry is that record's *best
  achievable* score, and the small sea is still 22 % worse than the large one.

A single coefficient shifts the whole error curve up or down; it cannot tilt it.
Re-sweeping `C_J` therefore has nowhere to go, and could not remove this
gradient even if it did. Releasing the floor changes nothing either: the
smallest record scores 17.877 with the floor at 0.001 against 17.883 with it at
the deployed 0.15, so `MIN_R_S` is not costing anything on the deployed
envelope.

The same measurement was already available for OU-II and says the same thing.
From `reports/results/ou2_pseudo_mse_scale/raw.csv`, sweeping the pseudo-
measurement scale `C_P`, the per-record argmin moves only from 0.098 on the two
smallest seas to 0.1146 on the largest -- a 17 % drift across the envelope --
while the per-record *minimum* still falls monotonically, 6.657 to 5.927 %H_s.

## 2. What the gradient actually is: the amplitude axis

The deployed record set raises `H_s` and `T_z` together -- 0.27 m at
`T_z = 2.6 s` up to 8.5 m at `T_z = 8.3 s` -- so it cannot separate the two.
Rescaling one record's amplitude with its period and spectrum held fixed does.

JONSWAP H_s = 1.5 m, `T_z = 4.2 s`, amplitude scaled over 32x
(`scale_wave_motion`, which carries displacement, velocity, acceleration and the
attitude columns together so the record stays kinematically consistent):

| scale | H_s | z RMS (m) | z % ref RMS | roll (deg) | local slope |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.0625 | 0.094 | 0.0137 | 58.245 | 0.215 | |
| 0.125 | 0.189 | 0.0215 | 45.589 | 0.223 | -0.353 |
| 0.25 | 0.378 | 0.0300 | 31.772 | 0.221 | -0.521 |
| 0.5 | 0.755 | 0.0431 | 22.852 | 0.211 | -0.475 |
| 1 | 1.510 | 0.0636 | 16.834 | 0.190 | -0.441 |
| 2 | 3.020 | 0.0966 | 12.796 | 0.184 | -0.396 |

The slope is `d ln(relative error) / d ln(amplitude)`, and it is -0.44 over the
whole range. The absolute error falls as `H_s^0.56` -- it falls, but far more
slowly than the sea does, so the ratio rises.

That is the predicted shape, not a defect. The bias--variance optimum the
SpectralMSE law is derived from balances a drift term set by the *sensor*
against a distortion term set by the *sea*,

    J(w_R) = 3 q_eff / (2 w_R^3)  +  4 m_-4 w_R^4,

and at its own optimum leaves a residual relative error going as
`q_eff^(2/7) sigma_a^(-4/7) tau^(-2/7)`. With `tau` held fixed that is a slope
of -4/7 = -0.571 against amplitude. The measured -0.44 is shallower because a
scale-invariant in-band component dilutes it (section 3), but the sign and the
order are what the law predicts for a correctly tuned filter. `q_eff` does not
shrink when the sea does; a filter sitting on the MSE optimum still spends a
larger fraction of a small sea on the sensor floor.

Scale factors above 2 are not scored: `scale_wave_motion` scales the attitude
columns too, so at 4x the record carries a roll RMS error of 0.89 deg against
0.18 deg at 1x, and at 8x the estimator loses attitude altogether. Below 2x the
roll column stays flat at 0.18 to 0.22 deg, which is what makes the amplitude
axis clean over the range that is scored.

Note the size of this effect against the size of the question. Across the
deployed records `sigma_a` rises 4.1x, which at slope -0.44 alone would make the
smallest sea about 1.9x worse than the largest. The observed gradient is 1.22x.
The period axis, which grows with `H_s` over the same records, cancels most of
it. What is visible in the committed tables is the residue of two large and
nearly opposite effects, not a single uncorrected one.

## 3. How much of it is the sensor

Scoring the deployed records with and without the simulated IMU noise splits the
error into a sensor-driven part and a sensor-independent one (subtracted in
quadrature):

| record | with noise | no noise | noise part |
| --- | ---: | ---: | ---: |
| JONSWAP H_s = 0.27 m | 17.883 | 13.240 | 12.022 |
| JONSWAP H_s = 1.50 m | 16.786 | 11.615 | 12.119 |
| JONSWAP H_s = 4.00 m | 15.741 | 11.189 | 11.072 |
| JONSWAP H_s = 8.50 m | 14.686 | 11.407 | 9.250 |

The sensor-driven part carries the cleaner gradient, 12.02 down to 9.25, a
factor 1.30. The sensor-independent part is nearly flat and not monotone --
13.24, 11.62, 11.19, 11.41 -- with the smallest record the only outlier, and it
is the shortest-period record rather than merely the smallest. That residue is
in-band regularization bias and OU model mismatch, not drift.

So roughly two-thirds of the gradient is the sensor floor against a shrinking
signal, which is irreducible at fixed hardware, and the remainder is a
single-record model-mismatch effect on the 3 s sea rather than a trend.

## 4. One thing that is genuinely mistuned, off the deployed envelope

The amplitude sweep exposes a separate defect that the deployed record set
cannot see. `ACC_NOISE_FLOOR_SIGMA_DEFAULT = 0.12` m/s^2 is subtracted in
variance from the band acceleration before `sigma_a` is formed. The simulator's
white accelerometer noise is 0.0148 m/s^2, so on the deployed records the
subtraction is small -- on JONSWAP H_s = 0.27 m it removes about 9 % of the
band variance, 4.6 % of `sigma_a`. Down the amplitude sweep it stops being small.
Comparing `sigma_applied / c_sigma` against the value implied by linear scaling
from the unscaled record:

| H_s | `sigma_a` implied | `sigma_a` estimated | ratio |
| ---: | ---: | ---: | ---: |
| 3.000 | 1.870 | 1.872 | 1.00 |
| 1.500 | 0.935 | 0.935 | 1.00 |
| 0.750 | 0.468 | 0.431 | 0.92 |
| 0.375 | 0.234 | 0.115 | 0.49 |
| 0.188 | 0.117 | 0.038 | 0.33 |
| 0.094 | 0.058 | 0.001 | 0.02 |

At the bottom of the range the subtraction removes the entire wave signal and
`var_wave` lands on its 1e-6 floor. Since `r_S ~ sigma_a^(6/7)`, the schedule
collapses with it and the filter runs a far stronger anchor than the law asks
for.

This does not affect any committed record -- the deployed set pairs its small
`H_s` with short periods, which keeps `sigma_a` at 0.44 m/s^2 even on the 0.27 m
sea -- so it is not the cause of the gradient in section 1. It is a real
exposure for a *low long-period swell*, where a small amplitude and a long
period both push `sigma_a` down, and that sea is not in the record set. The
floor is reachable as `OU_ACC_NOISE_FLOOR_SIGMA` and should be re-derived
against the deployed sensor rather than left at a figure that is 8x the
simulator's own white noise.

## 5. Reproducing

```sh
make fetch-sim-data
make -C tests/kalman_ou_iii build

# section 1: per-record C_J optimum, r_S floor released
python3 tools/ou_low_sea_error_study.py cj

# section 2: amplitude alone, period and spectrum fixed
python3 tools/ou_low_sea_error_study.py amp --scales 0.0625,0.125,0.25,0.5,1,2

# section 3: with and without the simulated IMU noise
python3 tools/ou_low_sea_error_study.py noise
```
