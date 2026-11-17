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

Source commit used for the replay: `71f410938ca1a8f393979f4933933868fb612461`.

## Pooled RMS across the eight equal-duration records

Because every record contributes the same 900 s window at the same sample
rate, `sqrt(mean(record_RMS^2))` is the exact pooled RMS over their concatenation.

| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | 0.2713 | 0.1889 | 0.2347 | 0.4054 | 18.693 | 0.0240 | 0.0204 | 0.1059 | 0.003830 | 0.0000246 |
| OU-III | 0.2053 | 0.1757 | 0.1362 | 0.3026 | 10.846 | 0.0612 | 0.0614 | 0.2694 | 0.014243 | 0.0000279 |
| TFG | 0.2278 | 0.2007 | 0.0776 | 0.3133 | 6.179 | 0.2432 | 0.1131 | 0.1943 | 0.045259 | 0.0000163 |

## Per-record RMS

| Family | Sea | Hs [m] | X [m] | Y [m] | Z [m] | 3D [m] | Z / Hs [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | JONSWAP | 0.27 | 0.0022 | 0.0013 | 0.0093 | 0.0096 | 3.438 | 0.0093 | 0.0111 | 0.0588 | 0.002471 | 0.0000112 |
| OU-II | JONSWAP | 1.50 | 0.0408 | 0.0224 | 0.0708 | 0.0847 | 4.718 | 0.0126 | 0.0099 | 0.0152 | 0.002504 | 0.0000079 |
| OU-II | JONSWAP | 4.00 | 0.1970 | 0.1051 | 0.2001 | 0.2999 | 5.003 | 0.0115 | 0.0223 | 0.0534 | 0.002720 | 0.0000180 |
| OU-II | JONSWAP | 8.50 | 0.5132 | 0.3629 | 0.4227 | 0.7575 | 4.973 | 0.0349 | 0.0297 | 0.1402 | 0.004403 | 0.0000350 |
| OU-II | PM-Stokes | 0.27 | 0.0023 | 0.0014 | 0.0094 | 0.0098 | 3.481 | 0.0120 | 0.0099 | 0.0586 | 0.002652 | 0.0000112 |
| OU-II | PM-Stokes | 1.50 | 0.0387 | 0.0243 | 0.0685 | 0.0824 | 4.570 | 0.0181 | 0.0117 | 0.0213 | 0.003503 | 0.0000084 |
| OU-II | PM-Stokes | 4.00 | 0.1894 | 0.1056 | 0.1971 | 0.2931 | 4.928 | 0.0117 | 0.0237 | 0.0639 | 0.003047 | 0.0000193 |
| OU-II | PM-Stokes | 8.50 | 0.4975 | 0.3612 | 0.4160 | 0.7423 | 4.894 | 0.0490 | 0.0305 | 0.2356 | 0.007041 | 0.0000504 |
| OU-III | JONSWAP | 0.27 | 0.0019 | 0.0011 | 0.0044 | 0.0049 | 1.629 | 0.0160 | 0.0196 | 0.0695 | 0.004317 | 0.0000106 |
| OU-III | JONSWAP | 1.50 | 0.0285 | 0.0219 | 0.0372 | 0.0517 | 2.483 | 0.0104 | 0.0197 | 0.1033 | 0.003567 | 0.0000112 |
| OU-III | JONSWAP | 4.00 | 0.1434 | 0.0875 | 0.1080 | 0.1997 | 2.699 | 0.0254 | 0.0766 | 0.1663 | 0.013423 | 0.0000203 |
| OU-III | JONSWAP | 8.50 | 0.3749 | 0.3061 | 0.2375 | 0.5391 | 2.794 | 0.0930 | 0.0856 | 0.4137 | 0.020592 | 0.0000372 |
| OU-III | PM-Stokes | 0.27 | 0.0021 | 0.0012 | 0.0047 | 0.0053 | 1.730 | 0.0199 | 0.0191 | 0.0707 | 0.004703 | 0.0000107 |
| OU-III | PM-Stokes | 1.50 | 0.0297 | 0.0276 | 0.0385 | 0.0559 | 2.563 | 0.0110 | 0.0206 | 0.1062 | 0.003651 | 0.0000122 |
| OU-III | PM-Stokes | 4.00 | 0.1483 | 0.1064 | 0.1133 | 0.2148 | 2.833 | 0.0219 | 0.0804 | 0.1676 | 0.013776 | 0.0000234 |
| OU-III | PM-Stokes | 8.50 | 0.3905 | 0.3646 | 0.2540 | 0.5916 | 2.988 | 0.1389 | 0.0947 | 0.5675 | 0.027609 | 0.0000580 |
| TFG | JONSWAP | 0.27 | 0.0015 | 0.0015 | 0.0016 | 0.0027 | 0.607 | 0.0351 | 0.0569 | 0.1382 | 0.011473 | 0.0000033 |
| TFG | JONSWAP | 1.50 | 0.0331 | 0.0244 | 0.0160 | 0.0441 | 1.063 | 0.0733 | 0.0482 | 0.1421 | 0.014870 | 0.0000081 |
| TFG | JONSWAP | 4.00 | 0.1551 | 0.1037 | 0.0542 | 0.1943 | 1.356 | 0.0657 | 0.1291 | 0.2565 | 0.024282 | 0.0000155 |
| TFG | JONSWAP | 8.50 | 0.4138 | 0.3553 | 0.1250 | 0.5596 | 1.470 | 0.3981 | 0.1286 | 0.1772 | 0.070337 | 0.0000235 |
| TFG | PM-Stokes | 0.27 | 0.0016 | 0.0016 | 0.0017 | 0.0028 | 0.630 | 0.0359 | 0.0558 | 0.1362 | 0.011386 | 0.0000032 |
| TFG | PM-Stokes | 1.50 | 0.0340 | 0.0288 | 0.0201 | 0.0489 | 1.337 | 0.0665 | 0.0596 | 0.1643 | 0.015099 | 0.0000098 |
| TFG | PM-Stokes | 4.00 | 0.1593 | 0.1225 | 0.0663 | 0.2116 | 1.657 | 0.0959 | 0.1347 | 0.2534 | 0.027662 | 0.0000179 |
| TFG | PM-Stokes | 8.50 | 0.4383 | 0.4107 | 0.1566 | 0.6208 | 1.843 | 0.5376 | 0.1967 | 0.2364 | 0.096822 | 0.0000286 |

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
