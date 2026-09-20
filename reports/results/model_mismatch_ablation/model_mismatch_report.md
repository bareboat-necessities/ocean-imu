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

Source commit used for the replay: `ffdd544dc4bc6ba18bc781c5fdef33eeafe2c290`.

## Pooled RMS across the eight equal-duration records

Because every record contributes the same 900 s window at the same sample
rate, `sqrt(mean(record_RMS^2))` is the exact pooled RMS over their concatenation.

| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | 0.3084 | 0.1920 | 0.2486 | 0.4402 | 20.456 | 0.0502 | 0.0324 | 0.2541 | 0.008441 | 0.0000219 |
| OU-III | 0.2402 | 0.1362 | 0.1446 | 0.3117 | 11.893 | 0.1006 | 0.0831 | 0.2650 | 0.021381 | 0.0000233 |
| TFG | 0.2342 | 0.1528 | 0.1453 | 0.3152 | 11.953 | 0.0644 | 0.1419 | 0.2516 | 0.026306 | 0.0000123 |

## Per-record RMS

| Family | Sea | Hs [m] | X [m] | Y [m] | Z [m] | 3D [m] | Z / Hs [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | JONSWAP | 0.27 | 0.0049 | 0.0024 | 0.0139 | 0.0149 | 5.139 | 0.0800 | 0.0394 | 0.1898 | 0.015186 | 0.0000041 |
| OU-II | JONSWAP | 1.50 | 0.0609 | 0.0275 | 0.0806 | 0.1047 | 5.372 | 0.0178 | 0.0175 | 0.2051 | 0.002735 | 0.0000105 |
| OU-II | JONSWAP | 4.00 | 0.2322 | 0.1307 | 0.2125 | 0.3409 | 5.314 | 0.0246 | 0.0262 | 0.2156 | 0.002982 | 0.0000216 |
| OU-II | JONSWAP | 8.50 | 0.5737 | 0.3673 | 0.4436 | 0.8129 | 5.218 | 0.0464 | 0.0379 | 0.2955 | 0.005762 | 0.0000321 |
| OU-II | PM-Stokes | 0.27 | 0.0049 | 0.0025 | 0.0135 | 0.0146 | 5.017 | 0.0771 | 0.0349 | 0.1737 | 0.014411 | 0.0000042 |
| OU-II | PM-Stokes | 1.50 | 0.0580 | 0.0276 | 0.0788 | 0.1017 | 5.254 | 0.0229 | 0.0214 | 0.1383 | 0.003583 | 0.0000108 |
| OU-II | PM-Stokes | 4.00 | 0.2223 | 0.1314 | 0.2101 | 0.3329 | 5.252 | 0.0278 | 0.0300 | 0.4636 | 0.003017 | 0.0000297 |
| OU-II | PM-Stokes | 8.50 | 0.5667 | 0.3522 | 0.4421 | 0.8004 | 5.201 | 0.0591 | 0.0431 | 0.2014 | 0.007755 | 0.0000346 |
| OU-III | JONSWAP | 0.27 | 0.0051 | 0.0020 | 0.0080 | 0.0097 | 2.957 | 0.1384 | 0.0892 | 0.3310 | 0.028113 | 0.0000036 |
| OU-III | JONSWAP | 1.50 | 0.0452 | 0.0175 | 0.0453 | 0.0663 | 3.019 | 0.0158 | 0.0834 | 0.0997 | 0.014074 | 0.0000108 |
| OU-III | JONSWAP | 4.00 | 0.1783 | 0.0938 | 0.1218 | 0.2354 | 3.044 | 0.0443 | 0.0829 | 0.1850 | 0.015081 | 0.0000226 |
| OU-III | JONSWAP | 8.50 | 0.4354 | 0.2497 | 0.2502 | 0.5608 | 2.944 | 0.1144 | 0.0703 | 0.1717 | 0.020995 | 0.0000338 |
| OU-III | PM-Stokes | 0.27 | 0.0054 | 0.0022 | 0.0080 | 0.0099 | 2.946 | 0.1440 | 0.0846 | 0.3188 | 0.028529 | 0.0000037 |
| OU-III | PM-Stokes | 1.50 | 0.0457 | 0.0189 | 0.0462 | 0.0677 | 3.082 | 0.0186 | 0.0893 | 0.1384 | 0.015069 | 0.0000120 |
| OU-III | PM-Stokes | 4.00 | 0.1784 | 0.1036 | 0.1260 | 0.2418 | 3.151 | 0.0350 | 0.0838 | 0.4260 | 0.014079 | 0.0000316 |
| OU-III | PM-Stokes | 8.50 | 0.4517 | 0.2568 | 0.2637 | 0.5827 | 3.103 | 0.1554 | 0.0795 | 0.2762 | 0.027591 | 0.0000374 |
| TFG | JONSWAP | 0.27 | 0.0031 | 0.0014 | 0.0084 | 0.0091 | 3.104 | 0.0624 | 0.1338 | 0.3920 | 0.025263 | 0.0000038 |
| TFG | JONSWAP | 1.50 | 0.0480 | 0.0222 | 0.0460 | 0.0701 | 3.065 | 0.0317 | 0.0808 | 0.0623 | 0.014374 | 0.0000058 |
| TFG | JONSWAP | 4.00 | 0.1780 | 0.1027 | 0.1221 | 0.2390 | 3.053 | 0.0423 | 0.1436 | 0.1664 | 0.025348 | 0.0000106 |
| TFG | JONSWAP | 8.50 | 0.4201 | 0.2800 | 0.2509 | 0.5638 | 2.952 | 0.0780 | 0.1629 | 0.1546 | 0.030205 | 0.0000205 |
| TFG | PM-Stokes | 0.27 | 0.0033 | 0.0016 | 0.0085 | 0.0092 | 3.132 | 0.0612 | 0.1290 | 0.3757 | 0.024426 | 0.0000039 |
| TFG | PM-Stokes | 1.50 | 0.0485 | 0.0239 | 0.0467 | 0.0714 | 3.112 | 0.0296 | 0.0972 | 0.0832 | 0.016913 | 0.0000063 |
| TFG | PM-Stokes | 4.00 | 0.1793 | 0.1138 | 0.1268 | 0.2474 | 3.171 | 0.0520 | 0.1558 | 0.1749 | 0.027907 | 0.0000121 |
| TFG | PM-Stokes | 8.50 | 0.4404 | 0.2897 | 0.2655 | 0.5902 | 3.123 | 0.1142 | 0.1980 | 0.3448 | 0.038403 | 0.0000206 |

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
