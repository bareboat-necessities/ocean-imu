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

Source commit used for the replay: `9de0543507509cabb478b751b83bf5c74711764d`.

## Pooled RMS across the eight equal-duration records

Because every record contributes the same 900 s window at the same sample
rate, `sqrt(mean(record_RMS^2))` is the exact pooled RMS over their concatenation.

| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | 0.3279 | 0.2049 | 0.2575 | 0.4645 | 21.189 | 0.1843 | 0.1796 | 0.5196 | 0.043224 | 0.0000174 |
| OU-III | 0.2106 | 0.1361 | 0.1444 | 0.2894 | 11.876 | 0.0951 | 0.0569 | 0.2467 | 0.017896 | 0.0000212 |
| TFG | 0.2342 | 0.1521 | 0.1456 | 0.3150 | 11.981 | 0.0618 | 0.1159 | 0.2120 | 0.021895 | 0.0000122 |

## Per-record RMS

| Family | Sea | Hs [m] | X [m] | Y [m] | Z [m] | 3D [m] | Z / Hs [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | JONSWAP | 0.27 | 0.0057 | 0.0028 | 0.0242 | 0.0251 | 8.978 | 0.2447 | 0.3612 | 1.0063 | 0.074677 | 0.0000079 |
| OU-II | JONSWAP | 1.50 | 0.0690 | 0.0311 | 0.0832 | 0.1125 | 5.543 | 0.0319 | 0.0249 | 0.1754 | 0.005915 | 0.0000091 |
| OU-II | JONSWAP | 4.00 | 0.2540 | 0.1416 | 0.2192 | 0.3641 | 5.480 | 0.0265 | 0.0658 | 0.1358 | 0.010800 | 0.0000173 |
| OU-II | JONSWAP | 8.50 | 0.6071 | 0.3905 | 0.4599 | 0.8559 | 5.410 | 0.2402 | 0.0593 | 0.2941 | 0.040101 | 0.0000270 |
| OU-II | PM-Stokes | 0.27 | 0.0058 | 0.0029 | 0.0209 | 0.0219 | 7.729 | 0.2370 | 0.3278 | 0.9228 | 0.069226 | 0.0000084 |
| OU-II | PM-Stokes | 1.50 | 0.0657 | 0.0312 | 0.0812 | 0.1090 | 5.412 | 0.0329 | 0.0312 | 0.1055 | 0.006479 | 0.0000090 |
| OU-II | PM-Stokes | 4.00 | 0.2425 | 0.1422 | 0.2161 | 0.3546 | 5.403 | 0.0467 | 0.0749 | 0.3456 | 0.013597 | 0.0000203 |
| OU-II | PM-Stokes | 8.50 | 0.5991 | 0.3757 | 0.4581 | 0.8426 | 5.390 | 0.3048 | 0.0709 | 0.1717 | 0.050903 | 0.0000259 |
| OU-III | JONSWAP | 0.27 | 0.0041 | 0.0019 | 0.0080 | 0.0092 | 2.959 | 0.1355 | 0.0766 | 0.3030 | 0.026583 | 0.0000031 |
| OU-III | JONSWAP | 1.50 | 0.0391 | 0.0177 | 0.0452 | 0.0624 | 3.016 | 0.0209 | 0.0384 | 0.1606 | 0.006440 | 0.0000082 |
| OU-III | JONSWAP | 4.00 | 0.1510 | 0.0897 | 0.1217 | 0.2137 | 3.042 | 0.0408 | 0.0626 | 0.1560 | 0.011658 | 0.0000180 |
| OU-III | JONSWAP | 8.50 | 0.3834 | 0.2512 | 0.2498 | 0.5220 | 2.939 | 0.1031 | 0.0376 | 0.2548 | 0.016477 | 0.0000331 |
| OU-III | PM-Stokes | 0.27 | 0.0043 | 0.0020 | 0.0080 | 0.0093 | 2.948 | 0.1413 | 0.0725 | 0.2922 | 0.027117 | 0.0000031 |
| OU-III | PM-Stokes | 1.50 | 0.0395 | 0.0191 | 0.0462 | 0.0637 | 3.079 | 0.0239 | 0.0426 | 0.0934 | 0.007080 | 0.0000094 |
| OU-III | PM-Stokes | 4.00 | 0.1510 | 0.0994 | 0.1259 | 0.2203 | 3.149 | 0.0333 | 0.0651 | 0.3940 | 0.010783 | 0.0000280 |
| OU-III | PM-Stokes | 8.50 | 0.3991 | 0.2578 | 0.2633 | 0.5432 | 3.097 | 0.1400 | 0.0437 | 0.1746 | 0.022466 | 0.0000348 |
| TFG | JONSWAP | 0.27 | 0.0031 | 0.0015 | 0.0084 | 0.0091 | 3.121 | 0.0466 | 0.1005 | 0.3088 | 0.018934 | 0.0000039 |
| TFG | JONSWAP | 1.50 | 0.0480 | 0.0221 | 0.0461 | 0.0701 | 3.075 | 0.0367 | 0.0651 | 0.0763 | 0.012144 | 0.0000057 |
| TFG | JONSWAP | 4.00 | 0.1776 | 0.1033 | 0.1222 | 0.2391 | 3.056 | 0.0284 | 0.1034 | 0.0953 | 0.017670 | 0.0000105 |
| TFG | JONSWAP | 8.50 | 0.4201 | 0.2780 | 0.2513 | 0.5630 | 2.957 | 0.0836 | 0.1446 | 0.1273 | 0.027684 | 0.0000204 |
| TFG | PM-Stokes | 0.27 | 0.0033 | 0.0016 | 0.0085 | 0.0093 | 3.148 | 0.0455 | 0.0968 | 0.2951 | 0.018280 | 0.0000039 |
| TFG | PM-Stokes | 1.50 | 0.0485 | 0.0238 | 0.0469 | 0.0715 | 3.126 | 0.0348 | 0.0791 | 0.0466 | 0.014105 | 0.0000062 |
| TFG | PM-Stokes | 4.00 | 0.1790 | 0.1146 | 0.1272 | 0.2477 | 3.181 | 0.0369 | 0.1146 | 0.2153 | 0.019916 | 0.0000121 |
| TFG | PM-Stokes | 8.50 | 0.4406 | 0.2882 | 0.2663 | 0.5900 | 3.133 | 0.1207 | 0.1810 | 0.3121 | 0.036342 | 0.0000205 |

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
