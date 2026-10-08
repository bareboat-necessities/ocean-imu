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

Source commit used for the replay: `de7f2861a7f3c157baac3fec26cfa7274e2bf6fc`.

## Pooled RMS across the eight equal-duration records

Because every record contributes the same 900 s window at the same sample
rate, `sqrt(mean(record_RMS^2))` is the exact pooled RMS over their concatenation.

| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | 0.3084 | 0.1920 | 0.2486 | 0.4402 | 20.456 | 0.0502 | 0.0324 | 0.2541 | 0.008441 | 0.0000219 |
| OU-III | 0.2402 | 0.1362 | 0.1446 | 0.3117 | 11.893 | 0.1006 | 0.0831 | 0.2650 | 0.021381 | 0.0000233 |
| TFG | 0.2579 | 0.1707 | 0.1452 | 0.3416 | 11.946 | 0.0610 | 0.1359 | 0.2420 | 0.024946 | 0.0000192 |

## Per-record RMS

| Family | Sea | Hs [m] | X [m] | Y [m] | Z [m] | 3D [m] | Z / Hs [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | JONSWAP | 0.27 | 0.0049 | 0.0024 | 0.0139 | 0.0149 | 5.139 | 0.0800 | 0.0394 | 0.1898 | 0.015186 | 0.0000041 |
| OU-II | JONSWAP | 1.50 | 0.0609 | 0.0275 | 0.0806 | 0.1047 | 5.372 | 0.0178 | 0.0175 | 0.2051 | 0.002735 | 0.0000105 |
| OU-II | JONSWAP | 4.00 | 0.2322 | 0.1307 | 0.2125 | 0.3409 | 5.313 | 0.0246 | 0.0262 | 0.2156 | 0.002982 | 0.0000216 |
| OU-II | JONSWAP | 8.50 | 0.5737 | 0.3673 | 0.4436 | 0.8129 | 5.218 | 0.0464 | 0.0379 | 0.2955 | 0.005762 | 0.0000321 |
| OU-II | PM-Stokes | 0.27 | 0.0049 | 0.0025 | 0.0135 | 0.0146 | 5.017 | 0.0771 | 0.0349 | 0.1737 | 0.014411 | 0.0000042 |
| OU-II | PM-Stokes | 1.50 | 0.0580 | 0.0276 | 0.0788 | 0.1017 | 5.254 | 0.0229 | 0.0214 | 0.1383 | 0.003584 | 0.0000108 |
| OU-II | PM-Stokes | 4.00 | 0.2223 | 0.1314 | 0.2101 | 0.3329 | 5.252 | 0.0278 | 0.0300 | 0.4636 | 0.003017 | 0.0000297 |
| OU-II | PM-Stokes | 8.50 | 0.5667 | 0.3522 | 0.4421 | 0.8004 | 5.201 | 0.0591 | 0.0431 | 0.2014 | 0.007756 | 0.0000346 |
| OU-III | JONSWAP | 0.27 | 0.0051 | 0.0020 | 0.0080 | 0.0097 | 2.957 | 0.1384 | 0.0892 | 0.3310 | 0.028113 | 0.0000036 |
| OU-III | JONSWAP | 1.50 | 0.0452 | 0.0175 | 0.0453 | 0.0663 | 3.019 | 0.0158 | 0.0833 | 0.0997 | 0.014073 | 0.0000108 |
| OU-III | JONSWAP | 4.00 | 0.1783 | 0.0938 | 0.1218 | 0.2354 | 3.044 | 0.0443 | 0.0829 | 0.1850 | 0.015081 | 0.0000226 |
| OU-III | JONSWAP | 8.50 | 0.4354 | 0.2497 | 0.2502 | 0.5608 | 2.944 | 0.1144 | 0.0703 | 0.1717 | 0.020995 | 0.0000338 |
| OU-III | PM-Stokes | 0.27 | 0.0054 | 0.0022 | 0.0080 | 0.0099 | 2.946 | 0.1440 | 0.0846 | 0.3188 | 0.028529 | 0.0000037 |
| OU-III | PM-Stokes | 1.50 | 0.0457 | 0.0189 | 0.0462 | 0.0677 | 3.082 | 0.0186 | 0.0893 | 0.1384 | 0.015069 | 0.0000120 |
| OU-III | PM-Stokes | 4.00 | 0.1784 | 0.1036 | 0.1260 | 0.2418 | 3.151 | 0.0350 | 0.0838 | 0.4260 | 0.014080 | 0.0000316 |
| OU-III | PM-Stokes | 8.50 | 0.4517 | 0.2568 | 0.2637 | 0.5827 | 3.103 | 0.1554 | 0.0795 | 0.2762 | 0.027591 | 0.0000374 |
| TFG | JONSWAP | 0.27 | 0.0033 | 0.0015 | 0.0084 | 0.0091 | 3.103 | 0.0634 | 0.1377 | 0.3997 | 0.025944 | 0.0000047 |
| TFG | JONSWAP | 1.50 | 0.0512 | 0.0244 | 0.0460 | 0.0730 | 3.066 | 0.0321 | 0.0850 | 0.0600 | 0.015041 | 0.0000096 |
| TFG | JONSWAP | 4.00 | 0.1926 | 0.1137 | 0.1221 | 0.2547 | 3.051 | 0.0401 | 0.1393 | 0.1678 | 0.024337 | 0.0000186 |
| TFG | JONSWAP | 8.50 | 0.4647 | 0.3137 | 0.2507 | 0.6142 | 2.950 | 0.0718 | 0.1406 | 0.1184 | 0.025806 | 0.0000310 |
| TFG | PM-Stokes | 0.27 | 0.0034 | 0.0017 | 0.0085 | 0.0093 | 3.132 | 0.0622 | 0.1326 | 0.3827 | 0.025054 | 0.0000047 |
| TFG | PM-Stokes | 1.50 | 0.0516 | 0.0262 | 0.0467 | 0.0744 | 3.113 | 0.0299 | 0.1013 | 0.0903 | 0.017584 | 0.0000103 |
| TFG | PM-Stokes | 4.00 | 0.1935 | 0.1256 | 0.1268 | 0.2632 | 3.170 | 0.0490 | 0.1530 | 0.1991 | 0.027080 | 0.0000207 |
| TFG | PM-Stokes | 8.50 | 0.4859 | 0.3237 | 0.2653 | 0.6413 | 3.121 | 0.1038 | 0.1767 | 0.2623 | 0.033922 | 0.0000311 |

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
