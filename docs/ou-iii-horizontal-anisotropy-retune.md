# Retuning the horizontal integral-regularizer scale

**Outcome: `R_S_xy_factor` moves from 1.0 to 0.72 in OU-III.** 3D displacement
RMS falls 3.33% pooled over the eight scored records and three IMU seed
triplets, with the same sign in all 24 record x seed cells. Vertical RMS is
unchanged to the fourth digit. OU-II's analogous `R_p0_xy_factor` reproduces the
effect independently and its optimum is 0.65; that change is measured here but
**not** applied, for the reason in section 6.

This note came out of asking why every family reports a higher vertical error as
a fraction of `H_s` on the small seas
([`ou-low-sea-percent-error.md`](ou-low-sea-percent-error.md)). It does not
answer that question -- nothing here moves the vertical channel at all. It is
what the screen for that question turned up instead.

## 1. What the screen ruled out first

Six knobs were swept one at a time on the four JONSWAP records, scoring vertical
RMS as a percentage of reference RMS (geometric mean over records; the deployed
point scores 16.231):

| knob | deployed | best on the grid | best score | gain |
| --- | ---: | ---: | ---: | ---: |
| `OU_III_RS_MSE_COEFF` (`C_J`) | 0.0538 | 0.0538 | 16.231 | none, argmin on every record |
| `OU_ACC_NOISE_FLOOR_SIGMA` | 0.12 | -- | 16.231 | none, flat to 5 digits |
| `OU_III_TAU_COEFF` (`c_tau`) | 1.0 | 1.0 | 16.231 | none, argmin on every record |
| `OU_III_ACC_BIAS_INIT_STD` | 0.004 | 0.08 | 16.209 | 0.14% |
| `OU_III_SIGMA_COEFF` (`c_sigma`) | 0.9 | 0.4 | 16.151 | 0.50%, but 3D moves the other way |
| `OU_III_R_S_XY_FACTOR` | 1.0 | 0.72 | 16.233 | none vertically -- see below |

Two of these deserve a line each, because they were the two the question
predicted.

**The pre-band accelerometer noise floor does nothing.** `sigma_a` is formed by
subtracting `(0.12 * sqrt(band gain))^2` from the band acceleration variance,
and 0.12 is eight times the simulator's white accelerometer noise, so
over-subtraction on a small sea looked like the obvious suspect. Sweeping it
from 0.0148 to 0.18 moves neither `sigma_applied` nor the score on any deployed
record: the band gain is small enough that the subtraction is a few percent of
the band variance even on the 0.27 m sea. The over-subtraction that
[`ou-low-sea-percent-error.md`](ou-low-sea-percent-error.md) section 4 measures
is real, but it needs a *low long-period* sea to bite, and there is not one in
the record set.

**The accelerometer-bias prior is over-confident and it barely matters.** The
MEKF starts with `sigma_bacc0 = 0.004` m/s^2 against a true bias drawn uniform
in +/- 0.0785, and against its own process model's stationary 0.025. Widening it
to 0.08 is worth 0.14% of vertical RMS and 2% of the bias error itself. The
horizontal accelerometer bias stays at 70-78% of the true bias whatever the
prior says, because it is close to unobservable from a gravity measurement on a
platform that does not manoeuvre -- which is a structural limit, not a tuning
one.

`c_tau = 1.0` being the argmin on every record separately is worth noting on its
own: `tau = c_tau * T_z / 2` is exactly self-similar across a 31x range of sea
state, and the sweep is steep either side (0.9 and 1.15 cost 3.5% and 7.2%).

## 2. Where the headroom actually was

`tools/ou_axis_rs_optimum.py` measures, per axis, the two quantities the
per-axis MSE optimum depends on: the drift-noise intensity `q_i` from the
posterior acceleration-error spectrum, and the true-displacement moment
`m_{-4,i}` from the reference columns. On the four JONSWAP records, geometric
means over records:

| axis | `q/q_z` | `m_-4/m_-4z` | `r_S*/r_S*_z` |
| --- | ---: | ---: | ---: |
| x | 10.26 | 0.685 | 1.004 |
| y | 5.19 | 0.315 | 0.685 |

The horizontal channels carry 5 to 10 times the vertical drift-noise intensity,
and the exponents are wildly unequal -- `r_S* ~ q^(1/14) m^(3/7)` -- so what
sets the answer is the displacement moment, not the noise.

The x/y split here is an artifact of the records: all eight are generated at
+/-30 degrees, so world x carries three times the horizontal displacement world
y does. A deployed filter cannot assume that, and the filter has one scalar for
both horizontal axes anyway. The direction-independent form of the same
calculation takes `m_h/m_z = 1/2` per horizontal axis -- deep-water orbits make
total horizontal and vertical displacement variance equal -- and gives

    rho_xy = (q_h/q_z)^(1/14) * (1/2)^(3/7) = 7.30^(1/14) * 0.743 = 0.86,

against a deployed 1.0. That is the reduced-model reference, not the answer;
the complete-MEKF sweep below puts the optimum lower still, the same way the
applied `C_R` sits off its analytical reference.

## 3. Why this was not already measured

`R_S_xy_factor` has a history, and every previous sweep of it missed this
interval.

- It was 0.36, a small-sea optimum applied to every sea state. Moving the
  operating point to the wave band made 1.0 better than 0.36 on every
  stationary record, by 7 to 27% of 3D RMS in the two largest seas. **Nothing
  between 0.36 and 1.0 was scored at that time.**
- The two-knob study in
  [`ou-iii-anisotropy-consistency.md`](ou-iii-anisotropy-consistency.md) did
  reach below 1, at arm G = `(S_sigma, rho_xy) = (1.87, 0.53)`, and found it
  8.7% worse. But that is at `S_sigma = 1.87`, where the inflated horizontal
  acceleration prior already put the horizontal normalized corner 23% above the
  vertical one; further tightening over-anchored. Its arms at `S_sigma = 1` were
  `rho_xy = 1` and `rho_xy = 1.87` only.
- `S_factor` then moved to 1.0. That removed the implicit 23% horizontal
  high-pass the `rho_xy = 1` decision had been sitting on top of, and nothing
  re-opened the interval below 1 afterwards.

So the deployed value is not a measured optimum at the deployed operating point.
It is a measured optimum at an operating point that no longer ships.

## 4. The sweep

`tests/kalman_ou_iii/kalman_ou_iii-sim`, all eight scored records, three IMU
seed triplets (`W3D_SEED=1,7,99`), trailing 900 s, adaptive tuning. Each arm is
a paired per-(record, seed) comparison against `rho_xy = 1`, n = 24 cells.
Negative is better; `*` marks a delta with the same sign in every cell.

| metric | 0.8 | **0.72** | 0.6 | 0.5 |
| --- | ---: | ---: | ---: | ---: |
| `disp_x_rms_m` | -1.49 | **-1.09** | +1.18 | +5.30 |
| `disp_y_rms_m` | -6.05\* | **-7.81\*** | -9.24\* | -8.70\* |
| `disp_z_rms_m` | -0.00 | **-0.01** | -0.01 | -0.02 |
| `disp_3d_rms_m` | -2.90\* | **-3.33\*** | -2.72 | -0.51 |
| `roll_rms_deg` | -0.14 | **-0.19** | -0.25 | -0.24 |
| `pitch_rms_deg` | +0.27 | **+0.42** | +0.76 | +1.23 |
| `yaw_rms_deg` | +0.04 | **+0.06** | +0.09 | +0.14 |
| `accel_bias_3d_rms_mps2` | +0.03 | **+0.05** | +0.09 | +0.15 |
| `dir_travel_correct_pct` | -0.00 | **-0.01** | -0.01 | +0.00 |

An interior minimum in 3D RMS at 0.72, bracketed at -2.90 and -2.72 either side.

**0.72 is taken rather than the pooled minimum, because it is the last point
where both horizontal axes gain.** Below it x turns over -- +1.18% at 0.6, +5.30%
at 0.5 -- while y keeps improving, and the parameter goes back to trading the two
horizontal axes against each other. That trade is exactly what the old 0.36 was
doing and what the earlier sweeps correctly rejected it for; it also depends on
the records' fixed +/-30 degree heading, which a deployed filter cannot rely on.
At 0.72 the gain is on both axes and does not need the heading to hold.

Vertical is untouched at -0.01%, so this is not a vertical-versus-horizontal
trade. Attitude moves under half a percent in both directions, and the direction
estimator does not move at all.

## 5. Quality gates

Re-cut with `tools/ou_regauge_gates.py --family ou_iii` under the rule in
[`quality-gate-regauge.md`](quality-gate-regauge.md). Against the bars this
change inherited, only pitch fails; every other bar still holds and is re-cut to
the new worst plus half a percent.

| gate | was | now | binding record |
| --- | ---: | ---: | --- |
| Z %Hs JONSWAP | 4.489 | 4.489 | jonswap H0.27 |
| Z %Hs PM-Stokes | 4.462 | 4.462 | pmstokes H0.27 |
| yaw deg | 0.9004 | 0.9023 | jonswap H1.5 |
| roll deg | 0.3513 | 0.3493 | pmstokes H4.0 |
| pitch deg | 0.1975 | **0.2007** | pmstokes H4.0 |
| 3D % JONSWAP | 13.98 | **12.77** | jonswap H8.5 |
| 3D % PM-Stokes | 14.51 | **13.43** | pmstokes H8.5 |
| acc Z bias % | 4.301 | 4.3 | pmstokes H8.5 |
| bias 3D % | 78.92 | 79.0 | pmstokes H4.0 |
| gyro 3D bias % | 15.76 | 15.73 | pmstokes H0.27 |

The two 3D bars come down 8.6% and 7.4%, which is the change measured on the
gate the change is aimed at. Pitch goes up 1.6% on one record, against a pooled
multi-seed pitch delta of +0.42%; yaw's bar re-cuts +0.2% while pooled yaw moves
+0.06%, inside a record whose yaw spans 1.05 to 6.57 deg across five IMU seeds.

While re-cutting these it turned out `tools/ou_regauge_gates.py` carried a
`shipped` table two revisions stale -- it still listed the pre-`S_factor` bars,
so its printed "was" column and its "SHIPPED GATE NOW FAILS" annotations were
both wrong. Its own comment says a stale entry "silently reports the wrong `was`
without changing any `is`", which is what happened. It is refreshed here.

## 6. The certificates, and why OU-II is measured but not moved

Two proof tools asserted the literal `float R_S_xy_factor_ = 1.0f;` and used
isotropy in a load-bearing way:

- `ou3_p4_nonlinear_word_certificate.py` lower-bounds the smallest eigenvalue of
  every correction measurement covariance, taking `MIN_R_S^2` for the S=0 update.
  With an anisotropic `diag(rho_xy r_S, rho_xy r_S, r_S)^2` the correct bound is
  `(min(rho_xy,1) * MIN_R_S)^2`.
- `ou3_p5_first_s_gain_certificate.py` bounds the first-S-to-attitude gain by
  `sqrt(lambda_max(P_theta)) * rho * sqrt(D)/(D+r)`, which is decreasing in `r`,
  so the same per-axis minimum is the conservative choice.

Both now read `rho_xy` from the deployed member instead of asserting a literal,
so a future retune does not silently invalidate either bound, and both still
pass: P4 `PASS`, P5 `PASS` with the gain retained.

**OU-II is deliberately left at 1.0.** The same sweep protocol -- eight records,
three seed triplets, 24 paired cells -- gives the same picture on
`OU_II_R_P0_XY_FACTOR`, with the optimum at 0.65 by the same
both-axes-must-gain rule:

| metric | 0.8 | **0.65** | 0.55 |
| --- | ---: | ---: | ---: |
| `disp_x_rms_m` | -1.77 | **-0.86** | +1.35 |
| `disp_y_rms_m` | -6.69\* | **-9.90\*** | -10.63\* |
| `disp_z_rms_m` | -0.01 | **-0.01** | -0.02 |
| `disp_3d_rms_m` | -2.95\* | **-3.56\*** | -2.80\* |

That is an independent family reproducing the effect at nearly the same value,
which is the strongest evidence here that it is structural rather than a fit to
one wrapper. It is not applied because it has not been carried through the rest
of the checklist this note ran for OU-III: its own gate re-cut, its own contract
test, and its own paper statements. Doing that is the obvious next change and
nothing measured so far argues against it.

TFG is untouched and unmeasured. Its `R_S_x_factor` and `R_S_y_factor` are
separate per-axis knobs already sitting at 1.15, so it is not the same starting
point and deserves its own sweep rather than an assumed carry-over.

## 7. What is still stale

The evidence bundles in `reports/results/ou_validation` and
`reports/results/ou_robustness` were produced at `rho_xy = 1`, as were the
paper's results, intro and conclusion paragraphs about isotropic regularization
(`w3d-results.tex-part`, `w3d-intro.tex-part`, `w3d-conclusion-summary.tex-part`).
Only the two parameter statements are updated here -- the table row in
`w3d-fus-methods.tex-part` and the tuning-points caption. Regenerating the
bundles and rewriting the horizontal-drift analysis around a non-isotropic
regularizer is a larger job than this change, and is the same outstanding item
[`ou-iii-anisotropy-consistency.md`](ou-iii-anisotropy-consistency.md) already
records against the `S_factor` move.

## 8. Reproducing

```sh
make fetch-sim-data
make -C tests/kalman_ou_iii build

# the knob screen of section 1, one knob at a time
python3 tools/ou_low_sea_error_study.py screen

# the per-axis MSE optimum of section 2
python3 tools/ou_axis_rs_optimum.py --glob 'tests/kalman_ou_iii/w3d_jonswap_*_ou3.csv'

# the paired multi-seed sweep of section 4
python3 tools/ou_low_sea_error_study.py xy
python3 tools/ou_low_sea_error_study.py xy --family OU_II

# the gates of section 5
python3 tools/ou_regauge_gates.py --family ou_iii
```
