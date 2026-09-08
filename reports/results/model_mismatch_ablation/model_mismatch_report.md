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

Source commit used for the replay: `456430d4d2cb2449dd0ae3c7781f71c5eedf599e`.

## Pooled RMS across the eight equal-duration records

Because every record contributes the same 900 s window at the same sample
rate, `sqrt(mean(record_RMS^2))` is the exact pooled RMS over their concatenation.

| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | 0.3257 | 0.2035 | 0.2577 | 0.4625 | 21.203 | 0.2103 | 0.2438 | 0.6817 | 0.054493 | 0.0000180 |
| OU-III | 0.2473 | 0.1595 | 0.1444 | 0.3278 | 11.877 | 0.0985 | 0.0644 | 0.2557 | 0.018979 | 0.0000204 |
| TFG | 0.2342 | 0.1521 | 0.1456 | 0.3150 | 11.981 | 0.0618 | 0.1159 | 0.2119 | 0.021895 | 0.0000122 |

## Per-record RMS

| Family | Sea | Hs [m] | X [m] | Y [m] | Z [m] | 3D [m] | Z / Hs [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | JONSWAP | 0.27 | 0.0052 | 0.0025 | 0.0250 | 0.0256 | 9.245 | 0.3208 | 0.4997 | 1.3650 | 0.101710 | 0.0000111 |
| OU-II | JONSWAP | 1.50 | 0.0686 | 0.0310 | 0.0832 | 0.1122 | 5.548 | 0.0299 | 0.0256 | 0.1745 | 0.005729 | 0.0000093 |
| OU-II | JONSWAP | 4.00 | 0.2524 | 0.1407 | 0.2193 | 0.3627 | 5.483 | 0.0261 | 0.0666 | 0.1351 | 0.010922 | 0.0000176 |
| OU-II | JONSWAP | 8.50 | 0.6028 | 0.3878 | 0.4602 | 0.8518 | 5.415 | 0.2400 | 0.0594 | 0.2945 | 0.040085 | 0.0000273 |
| OU-II | PM-Stokes | 0.27 | 0.0051 | 0.0026 | 0.0213 | 0.0221 | 7.889 | 0.3097 | 0.4531 | 1.2492 | 0.094003 | 0.0000119 |
| OU-II | PM-Stokes | 1.50 | 0.0653 | 0.0311 | 0.0812 | 0.1088 | 5.415 | 0.0303 | 0.0325 | 0.1039 | 0.006350 | 0.0000091 |
| OU-II | PM-Stokes | 4.00 | 0.2411 | 0.1412 | 0.2162 | 0.3533 | 5.404 | 0.0455 | 0.0761 | 0.3431 | 0.013726 | 0.0000205 |
| OU-II | PM-Stokes | 8.50 | 0.5952 | 0.3730 | 0.4584 | 0.8388 | 5.393 | 0.3044 | 0.0709 | 0.1722 | 0.050861 | 0.0000261 |
| OU-III | JONSWAP | 0.27 | 0.0032 | 0.0014 | 0.0080 | 0.0087 | 2.950 | 0.1424 | 0.0925 | 0.3428 | 0.029052 | 0.0000033 |
| OU-III | JONSWAP | 1.50 | 0.0492 | 0.0223 | 0.0452 | 0.0705 | 3.016 | 0.0218 | 0.0404 | 0.1584 | 0.006726 | 0.0000082 |
| OU-III | JONSWAP | 4.00 | 0.1834 | 0.1082 | 0.1217 | 0.2453 | 3.042 | 0.0428 | 0.0661 | 0.1541 | 0.012139 | 0.0000178 |
| OU-III | JONSWAP | 8.50 | 0.4465 | 0.2929 | 0.2498 | 0.5896 | 2.939 | 0.1047 | 0.0412 | 0.2492 | 0.016524 | 0.0000309 |
| OU-III | PM-Stokes | 0.27 | 0.0033 | 0.0016 | 0.0079 | 0.0087 | 2.940 | 0.1489 | 0.0877 | 0.3307 | 0.029571 | 0.0000033 |
| OU-III | PM-Stokes | 1.50 | 0.0497 | 0.0240 | 0.0462 | 0.0720 | 3.080 | 0.0252 | 0.0453 | 0.0907 | 0.007562 | 0.0000092 |
| OU-III | PM-Stokes | 4.00 | 0.1835 | 0.1198 | 0.1260 | 0.2528 | 3.149 | 0.0359 | 0.0701 | 0.3867 | 0.011596 | 0.0000277 |
| OU-III | PM-Stokes | 8.50 | 0.4664 | 0.3010 | 0.2633 | 0.6144 | 3.097 | 0.1415 | 0.0477 | 0.1663 | 0.022491 | 0.0000336 |
| TFG | JONSWAP | 0.27 | 0.0031 | 0.0015 | 0.0084 | 0.0091 | 3.121 | 0.0466 | 0.1005 | 0.3088 | 0.018933 | 0.0000039 |
| TFG | JONSWAP | 1.50 | 0.0480 | 0.0221 | 0.0461 | 0.0701 | 3.075 | 0.0367 | 0.0650 | 0.0763 | 0.012140 | 0.0000057 |
| TFG | JONSWAP | 4.00 | 0.1776 | 0.1033 | 0.1222 | 0.2391 | 3.056 | 0.0284 | 0.1033 | 0.0954 | 0.017653 | 0.0000105 |
| TFG | JONSWAP | 8.50 | 0.4201 | 0.2780 | 0.2513 | 0.5630 | 2.957 | 0.0837 | 0.1448 | 0.1277 | 0.027733 | 0.0000204 |
| TFG | PM-Stokes | 0.27 | 0.0033 | 0.0016 | 0.0085 | 0.0093 | 3.148 | 0.0455 | 0.0968 | 0.2951 | 0.018284 | 0.0000039 |
| TFG | PM-Stokes | 1.50 | 0.0485 | 0.0238 | 0.0469 | 0.0715 | 3.126 | 0.0348 | 0.0790 | 0.0466 | 0.014104 | 0.0000062 |
| TFG | PM-Stokes | 4.00 | 0.1790 | 0.1146 | 0.1272 | 0.2477 | 3.181 | 0.0369 | 0.1147 | 0.2152 | 0.019925 | 0.0000121 |
| TFG | PM-Stokes | 8.50 | 0.4406 | 0.2882 | 0.2663 | 0.5900 | 3.133 | 0.1207 | 0.1808 | 0.3116 | 0.036309 | 0.0000206 |

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
