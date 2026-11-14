# Noise-free model-mismatch ablation

This study runs the standard eight stationary JONSWAP / PM-Stokes wave
records through OU-II, OU-III, and TFG with simulator-side sensor
corruption disabled (`--no-noise`).  The filters themselves keep their
normal deployed covariance assumptions, adaptation laws, pseudo-measurements,
startup, and regularization.

Scoring uses the trailing **900 s** of each 1200 s record.
Magnetometer updates remain enabled, but the magnetic measurements are ideal.
The reported floor therefore contains model/estimator mismatch, intentional
regularization bias, residual adaptation/startup effects, and numerical error;
it is not a claim of pure plant-model mismatch in isolation.

Source commit used for the replay: `f7cd0fa140a2f9647ef33d16bf40c25258442158`.

## Pooled RMS across the eight equal-duration records

Because every record contributes the same 900 s window at the same sample
rate, `sqrt(mean(record_RMS^2))` is the exact pooled RMS over their concatenation.

| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | 0.3675 | 0.2410 | 0.2839 | 0.5232 | 23.227 | 0.4934 | 0.3070 | 0.6068 | 0.098902 | 0.0000505 |
| OU-III | 0.3023 | 0.1890 | 0.1669 | 0.3936 | 13.652 | 0.0931 | 0.0478 | 0.3615 | 0.015831 | 0.0000459 |
| TFG | 0.3113 | 0.2263 | 0.1638 | 0.4183 | 13.400 | 0.2058 | 0.3443 | 0.7180 | 0.068407 | 0.0000353 |

## Per-record RMS

| Family | Sea | Hs [m] | X [m] | Y [m] | Z [m] | 3D [m] | Z / Hs [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | JONSWAP | 0.27 | 0.0169 | 0.0074 | 0.0706 | 0.0730 | 26.146 | 0.9705 | 0.6987 | 1.2242 | 0.204219 | 0.0000642 |
| OU-II | JONSWAP | 1.50 | 0.0903 | 0.0469 | 0.0951 | 0.1393 | 6.342 | 0.0428 | 0.0144 | 0.2589 | 0.003184 | 0.0000426 |
| OU-II | JONSWAP | 4.00 | 0.2717 | 0.1642 | 0.2337 | 0.3942 | 5.843 | 0.0468 | 0.0236 | 0.3969 | 0.001727 | 0.0000460 |
| OU-II | JONSWAP | 8.50 | 0.6854 | 0.4590 | 0.5092 | 0.9694 | 5.990 | 0.0562 | 0.0387 | 0.5284 | 0.001961 | 0.0000523 |
| OU-II | PM-Stokes | 0.27 | 0.0158 | 0.0072 | 0.0646 | 0.0668 | 23.910 | 0.9956 | 0.5112 | 0.6633 | 0.191087 | 0.0000587 |
| OU-II | PM-Stokes | 1.50 | 0.0864 | 0.0465 | 0.0932 | 0.1353 | 6.213 | 0.0424 | 0.0142 | 0.2000 | 0.003054 | 0.0000420 |
| OU-II | PM-Stokes | 4.00 | 0.2668 | 0.1635 | 0.2313 | 0.3892 | 5.784 | 0.0470 | 0.0231 | 0.5888 | 0.001855 | 0.0000448 |
| OU-II | PM-Stokes | 8.50 | 0.6705 | 0.4423 | 0.5007 | 0.9466 | 5.891 | 0.0559 | 0.0382 | 0.3419 | 0.002000 | 0.0000492 |
| OU-III | JONSWAP | 0.27 | 0.0255 | 0.0084 | 0.0206 | 0.0338 | 7.619 | 0.1232 | 0.0507 | 0.3320 | 0.021868 | 0.0000401 |
| OU-III | JONSWAP | 1.50 | 0.0976 | 0.0450 | 0.0666 | 0.1265 | 4.443 | 0.0743 | 0.0406 | 0.1743 | 0.012489 | 0.0000444 |
| OU-III | JONSWAP | 4.00 | 0.2268 | 0.1368 | 0.1383 | 0.2988 | 3.456 | 0.0580 | 0.0386 | 0.4102 | 0.008237 | 0.0000453 |
| OU-III | JONSWAP | 8.50 | 0.5564 | 0.3525 | 0.2950 | 0.7217 | 3.471 | 0.0598 | 0.0568 | 0.3957 | 0.009315 | 0.0000540 |
| OU-III | PM-Stokes | 0.27 | 0.0238 | 0.0081 | 0.0193 | 0.0318 | 7.164 | 0.1699 | 0.0550 | 0.3600 | 0.029857 | 0.0000401 |
| OU-III | PM-Stokes | 1.50 | 0.0898 | 0.0432 | 0.0641 | 0.1185 | 4.271 | 0.0791 | 0.0301 | 0.1577 | 0.012592 | 0.0000439 |
| OU-III | PM-Stokes | 4.00 | 0.2234 | 0.1396 | 0.1381 | 0.2974 | 3.453 | 0.0585 | 0.0385 | 0.6086 | 0.008331 | 0.0000465 |
| OU-III | PM-Stokes | 8.50 | 0.5489 | 0.3453 | 0.2971 | 0.7133 | 3.495 | 0.0566 | 0.0626 | 0.2347 | 0.009836 | 0.0000507 |
| TFG | JONSWAP | 0.27 | 0.0143 | 0.0060 | 0.0160 | 0.0223 | 5.923 | 0.0695 | 0.0325 | 0.2121 | 0.011852 | 0.0000215 |
| TFG | JONSWAP | 1.50 | 0.0859 | 0.0435 | 0.0578 | 0.1123 | 3.852 | 0.0493 | 0.0193 | 0.1182 | 0.005866 | 0.0000241 |
| TFG | JONSWAP | 4.00 | 0.2383 | 0.1375 | 0.1349 | 0.3064 | 3.373 | 0.0546 | 0.1161 | 0.2061 | 0.020543 | 0.0000300 |
| TFG | JONSWAP | 8.50 | 0.5705 | 0.4102 | 0.2897 | 0.7600 | 3.408 | 0.3180 | 0.5472 | 0.9308 | 0.108019 | 0.0000503 |
| TFG | PM-Stokes | 0.27 | 0.0144 | 0.0058 | 0.0162 | 0.0224 | 5.983 | 0.0784 | 0.0421 | 0.2372 | 0.014146 | 0.0000216 |
| TFG | PM-Stokes | 1.50 | 0.0844 | 0.0442 | 0.0579 | 0.1114 | 3.858 | 0.0416 | 0.0189 | 0.0555 | 0.003846 | 0.0000240 |
| TFG | PM-Stokes | 4.00 | 0.2362 | 0.1433 | 0.1356 | 0.3077 | 3.389 | 0.0436 | 0.1007 | 0.3522 | 0.017000 | 0.0000315 |
| TFG | PM-Stokes | 8.50 | 0.5679 | 0.4452 | 0.2949 | 0.7796 | 3.469 | 0.4668 | 0.7885 | 1.7243 | 0.157058 | 0.0000588 |

## Figures

- `ou_model_mismatch_floor.svg`: pooled floor per family across the
  displacement, attitude, and estimator-generated bias channels.
- `ou_model_mismatch_scaling.svg`: per-record residual against `Hs`, with a
  slope-one guide and the `Hs`-normalized vertical residual.

Both are mirrored byte-for-byte into `doc/kalman_ou_iii/` for the article.

## Interpretation boundary

`--no-noise` removes the simulator's stochastic and calibration-error
sensor terms before they reach the filters.  It does **not** remove wave
nonlinearity, attitude/translation coupling, OU/TFG prior mismatch, the
integral pseudo-measurements, finite adaptation bandwidth, startup residue,
or discretization.  Those effects are intentionally what this ablation
leaves visible.
