# Splitting the horizontal regularizer into independent x and y knobs

**Outcome.** All three families now express the horizontal integral/position
regularizer as two independent per-axis factors. The combined single knob is
gone from every wrapper, setter, environment variable and sweep tool; there is
no compatibility alias.

| family | member | before | after |
| --- | --- | ---: | ---: |
| OU-III | `R_S_x_factor_`, `R_S_y_factor_` | `R_S_xy_factor_ = 0.72` | 0.72, 0.72 |
| OU-II | `R_p0_x_factor_`, `R_p0_y_factor_` | `R_p0_xy_factor_ = 1.0` | **0.72, 0.72** |
| TFG | `R_S_x_factor_`, `R_S_y_factor_` | 1.15, 1.15 (already split) | 1.15, 1.15 |

Two of the three keep or take 0.72. **TFG does not, and the reason is
measured, not assumed** -- see section 4. The horizontal anisotropy is not a
constant the three families share.

## 1. Why two knobs rather than one

`tools/ou_axis_rs_optimum.py` reads the per-axis MSE optimum
`r_S*_i ~ q_i^(1/14) m_{-4,i}^(3/7)` off the simulator output. On the four
JONSWAP records it does not put the two horizontal axes in the same place:

| axis | `q/q_z` | `m_-4/m_-4z` | `r_S*/r_S*_z` |
| --- | ---: | ---: | ---: |
| x | 10.26 | 0.685 | 1.004 |
| y | 5.19 | 0.315 | 0.685 |

A single scalar cannot sit at both 1.004 and 0.685, and that is exactly what
stopped the OU-III retune at 0.72 rather than at the pooled minimum near 0.6:
below 0.72 the scalar starts buying y at x's expense. With two knobs that
trade is expressible instead of forced.

**The defaults are nonetheless equal, deliberately.** The x/y asymmetry above
is a property of the record set, not of the sea: every record is generated at
+/-30 degrees, so world x carries about three times the horizontal displacement
world y does. A filter that does not know its own heading relative to the
dominant sea has no basis for splitting them, and tuning to this record set's
particular heading would be fitting the fixture. The split is there so a
deployment that *does* know can use it, and so the next study of this parameter
can measure the two axes separately rather than through one scalar.

## 2. The API change

The combined knob is removed outright:

| family | removed | added |
| --- | --- | --- |
| OU-III | `setRSXYFactor()`, `OU_III_R_S_XY_FACTOR` | `setRSXFactor()`, `setRSYFactor()`, `OU_III_R_S_X_FACTOR`, `OU_III_R_S_Y_FACTOR` |
| OU-II | `setR_p0_XYFactor()`, `OU_R_P0_XY_FACTOR`, `OU_II_R_P0_XY_FACTOR` | `setR_p0_XFactor()`, `setR_p0_YFactor()`, `OU_{,II_}R_P0_{X,Y}_FACTOR` |
| TFG | -- | already `setRSXFactor()` / `setRSYFactor()` |

A sweep that means to move both horizontal axes now has to say so twice. That
is the point: the old single name made "move the horizontal scale" and "move
one axis" indistinguishable, and the two are different experiments.

OU-II's setter ceiling also moves from 1 to 4, matching OU-III's. A ceiling of
1 encodes the assumption that a horizontal anchor can only ever be tighter than
the vertical one, which is the assumption these knobs exist to test -- and which
TFG's measurement in section 4 refutes outright.

`tools/ou_low_sea_error_study.py xy` drives both axes together, because the
common horizontal scale is what a heading-agnostic deployment can actually set.
`tools/ou_tuning_sweep.py` and `tools/ou_anisotropy_ablation.py` set both.

## 3. OU-III is unchanged; OU-II moves to 0.72

**OU-III.** The split is bit-for-bit neutral. On JONSWAP H_s = 8.5 m the split
build reproduces the scalar build exactly:

```
disp_x_rms_m=0.730157793  disp_y_rms_m=0.606045604
disp_3d_rms_m=0.998257697 disp_z_pct_refrms=14.6870136
```

No gate moves and no certificate changes value; only the shape of the bound
changes (section 5).

**OU-II moves from 1.0 to 0.72**, which is the value the earlier sweep in
[`ou-iii-horizontal-anisotropy-retune.md`](ou-iii-horizontal-anisotropy-retune.md)
measured but did not apply. Paired over the eight scored records and three IMU
seed triplets, n = 24 cells, against `1.0`; `*` marks the same sign in every
cell:

| metric | OU-II at 0.72 |
| --- | ---: |
| `disp_x_rms_m` | -1.59 |
| `disp_y_rms_m` | -8.66\* |
| `disp_z_rms_m` | -0.01 |
| `disp_3d_rms_m` | **-3.50\*** |
| `roll_rms_deg` | -0.31 |
| `pitch_rms_deg` | +0.34 |
| `yaw_rms_deg` | +0.06 |
| `accel_bias_3d_rms_mps2` | +0.05 |

0.72 is not OU-II's own argmin -- 0.65 is, at -3.56% -- but the basin is flat
enough that the difference is six hundredths of a point, so the shared value
costs essentially nothing and keeps one constant across the two OU wrappers.

OU-II's gates, re-cut with `tools/ou_regauge_gates.py --family ou_ii`: every bar
holds or comes down, and nothing is loosened.

| gate | was | now | binding record |
| --- | ---: | ---: | --- |
| Z %Hs JONSWAP | 6.672 | 6.674 | jonswap H0.27 |
| Z %Hs PM-Stokes | 6.605 | 6.605 | pmstokes H0.27 |
| yaw deg | 1.048 | 1.041 | jonswap H4.0 |
| roll deg | 0.3896 | 0.3866 | jonswap H4.0 |
| pitch deg | 0.3296 | 0.3237 | jonswap H8.5 |
| 3D % JONSWAP | 17.18 | **15.54** | jonswap H8.5 |
| 3D % PM-Stokes | 17.94 | **16.62** | pmstokes H8.5 |
| acc Z bias % | 4.67 | 4.663 | jonswap H8.5 |
| bias 3D % | 84.2 | 83.97 | jonswap H4.0 |

## 4. TFG measures the opposite sign, and keeps 1.15

TFG's 1.15 had been carried from its historical operating point rather than
measured. It has now been measured on the same protocol -- eight records, three
IMU seed triplets, 24 paired cells -- and it is right. Paired % change against
1.15:

| metric | 1.4 | **1.15** | 0.9 | 0.72 |
| --- | ---: | ---: | ---: | ---: |
| `disp_x_rms_m` | -3.78 | — | +7.80\* | +17.45\* |
| `disp_y_rms_m` | +7.68\* | — | -6.06\* | -8.10 |
| `disp_z_rms_m` | -0.07 | — | +0.22\* | +0.55\* |
| `disp_3d_rms_m` | +1.03 | — | +2.23 | **+7.31\*** |
| `pitch_rms_deg` | +5.06 | — | +0.81 | +3.27 |
| `accel_bias_3d_rms_mps2` | +5.00 | — | +0.99 | +1.51 |

1.15 is an interior minimum in 3D RMS: 1.4 and 0.9 are both worse, and 0.72 is
worse by 7.31% with every one of the 24 cells agreeing. It also costs 0.55% of
*vertical* RMS unanimously, which neither OU family's move did -- the two OU
retunes left vertical at -0.01%.

So the horizontal anisotropy is a per-family quantity with a per-family sign:
the two OU wrappers want a horizontal anchor tighter than vertical, TFG wants
one looser. Whatever the mechanism is -- TFG carries a different linear block
and a different pseudo-measurement -- it is not shared, and the number should
not be carried between families without re-running this sweep.

Carrying 0.72 into TFG anyway would have meant loosening four of its seven
regression bars to absorb the regression: 3D JONSWAP 20.43 -> 22.22, 3D
PM-Stokes 20.15 -> 22.22, acc Z bias 4.532 -> 4.636, bias 3D 155.6 -> 162.3.
Re-cutting bars upward to accommodate a measured loss is not a retune, so TFG
keeps 1.15 and its bars are untouched. Its per-axis knobs are of course still
reachable; `TFG_R_S_X_FACTOR` / `TFG_R_S_Y_FACTOR` set them.

## 5. Certificates

`ou3_p4_nonlinear_word_certificate.py` lower-bounds the smallest eigenvalue of
every correction measurement covariance, and
`ou3_p5_first_s_gain_certificate.py` bounds the first-S-to-attitude gain by
`sqrt(D)/(D+r)`, which decreases in `r`. Both used to read one horizontal
factor. With `diag(rho_x r_S, rho_y r_S, r_S)^2` the correct quantity is
`min(rho_x, rho_y, 1) * MIN_R_S`, which is what both now compute, from the
deployed members rather than an asserted literal. Both still pass, with the
first-S gain retained.

## 6. Reproducing

```sh
make fetch-sim-data
make -C tests/kalman_ou_iii build && make -C tests/kalman_ou_ii build \
  && make -C tests/kalman_tfg build

# per-axis MSE optimum of section 1
python3 tools/ou_axis_rs_optimum.py --glob 'tests/kalman_ou_iii/w3d_jonswap_*_ou3.csv'

# the paired sweeps of sections 3 and 4
python3 tools/ou_low_sea_error_study.py xy --family OU_II
python3 tools/ou_low_sea_error_study.py xy --family OU_III

# TFG, whose knobs the study driver does not carry: set both axes directly
for r in 0.72 0.9 1.15 1.4; do
  TFG_R_S_X_FACTOR=$r TFG_R_S_Y_FACTOR=$r W3D_VALIDATION_WINDOW_SEC=900 \
    tests/kalman_tfg/kalman_tfg-sim --input <record>
done

# gates
python3 tools/ou_regauge_gates.py --family ou_ii
python3 tools/ou_regauge_gates.py --family tfg
```
