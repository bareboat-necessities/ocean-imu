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
| OU-II | 0.3650 | 0.2280 | 0.2802 | 0.5135 | 22.995 | 0.0296 | 0.0235 | 0.0561 | 0.004604 | 0.0000217 |
| OU-III | 0.2981 | 0.1744 | 0.1649 | 0.3827 | 13.533 | 0.0649 | 0.0374 | 0.1449 | 0.012030 | 0.0000231 |
| TFG | 0.3080 | 0.2220 | 0.0662 | 0.3854 | 5.434 | 0.1565 | 0.2678 | 0.6750 | 0.053117 | 0.0000241 |

## Per-record RMS

| Family | Sea | Hs [m] | X [m] | Y [m] | Z [m] | 3D [m] | Z / Hs [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | JONSWAP | 0.27 | 0.0127 | 0.0047 | 0.0179 | 0.0225 | 6.641 | 0.0434 | 0.0072 | 0.0480 | 0.007471 | 0.0000110 |
| OU-II | JONSWAP | 1.50 | 0.0881 | 0.0399 | 0.0925 | 0.1338 | 6.164 | 0.0087 | 0.0145 | 0.0880 | 0.001762 | 0.0000218 |
| OU-II | JONSWAP | 4.00 | 0.2702 | 0.1508 | 0.2322 | 0.3869 | 5.804 | 0.0250 | 0.0231 | 0.0490 | 0.003872 | 0.0000212 |
| OU-II | JONSWAP | 8.50 | 0.6804 | 0.4368 | 0.5055 | 0.9536 | 5.947 | 0.0249 | 0.0378 | 0.0441 | 0.001230 | 0.0000304 |
| OU-II | PM-Stokes | 0.27 | 0.0120 | 0.0046 | 0.0171 | 0.0214 | 6.336 | 0.0487 | 0.0077 | 0.0521 | 0.008380 | 0.0000112 |
| OU-II | PM-Stokes | 1.50 | 0.0848 | 0.0395 | 0.0913 | 0.1308 | 6.088 | 0.0088 | 0.0147 | 0.0322 | 0.001919 | 0.0000155 |
| OU-II | PM-Stokes | 4.00 | 0.2656 | 0.1487 | 0.2302 | 0.3816 | 5.755 | 0.0275 | 0.0227 | 0.0763 | 0.004327 | 0.0000235 |
| OU-II | PM-Stokes | 8.50 | 0.6664 | 0.4208 | 0.4980 | 0.9323 | 5.859 | 0.0241 | 0.0376 | 0.0347 | 0.001235 | 0.0000295 |
| OU-III | JONSWAP | 0.27 | 0.0259 | 0.0074 | 0.0209 | 0.0342 | 7.758 | 0.0891 | 0.0218 | 0.1036 | 0.015583 | 0.0000149 |
| OU-III | JONSWAP | 1.50 | 0.0991 | 0.0370 | 0.0674 | 0.1254 | 4.492 | 0.0319 | 0.0273 | 0.1240 | 0.006357 | 0.0000233 |
| OU-III | JONSWAP | 4.00 | 0.2281 | 0.1247 | 0.1381 | 0.2943 | 3.452 | 0.0686 | 0.0362 | 0.0885 | 0.012565 | 0.0000213 |
| OU-III | JONSWAP | 8.50 | 0.5436 | 0.3232 | 0.2894 | 0.6955 | 3.405 | 0.0186 | 0.0506 | 0.2003 | 0.006310 | 0.0000295 |
| OU-III | PM-Stokes | 0.27 | 0.0243 | 0.0071 | 0.0197 | 0.0321 | 7.310 | 0.1167 | 0.0232 | 0.1187 | 0.020275 | 0.0000145 |
| OU-III | PM-Stokes | 1.50 | 0.0915 | 0.0357 | 0.0652 | 0.1179 | 4.346 | 0.0233 | 0.0255 | 0.1209 | 0.005065 | 0.0000235 |
| OU-III | PM-Stokes | 4.00 | 0.2260 | 0.1215 | 0.1376 | 0.2912 | 3.441 | 0.0716 | 0.0403 | 0.1065 | 0.013402 | 0.0000231 |
| OU-III | PM-Stokes | 8.50 | 0.5413 | 0.3250 | 0.2934 | 0.6962 | 3.452 | 0.0188 | 0.0576 | 0.2323 | 0.007781 | 0.0000296 |
| TFG | JONSWAP | 0.27 | 0.0145 | 0.0047 | 0.0066 | 0.0166 | 2.447 | 0.0333 | 0.0121 | 0.0372 | 0.005964 | 0.0000071 |
| TFG | JONSWAP | 1.50 | 0.0855 | 0.0384 | 0.0168 | 0.0952 | 1.120 | 0.0226 | 0.0167 | 0.0197 | 0.004143 | 0.0000126 |
| TFG | JONSWAP | 4.00 | 0.2391 | 0.1432 | 0.0460 | 0.2824 | 1.150 | 0.0240 | 0.1303 | 0.2980 | 0.022301 | 0.0000184 |
| TFG | JONSWAP | 8.50 | 0.5612 | 0.4113 | 0.1097 | 0.7044 | 1.291 | 0.2616 | 0.4732 | 1.1881 | 0.092638 | 0.0000387 |
| TFG | PM-Stokes | 0.27 | 0.0145 | 0.0048 | 0.0079 | 0.0172 | 2.926 | 0.0379 | 0.0143 | 0.0446 | 0.006846 | 0.0000075 |
| TFG | PM-Stokes | 1.50 | 0.0841 | 0.0395 | 0.0200 | 0.0951 | 1.332 | 0.0127 | 0.0211 | 0.0286 | 0.003210 | 0.0000130 |
| TFG | PM-Stokes | 4.00 | 0.2387 | 0.1420 | 0.0539 | 0.2829 | 1.348 | 0.0378 | 0.1726 | 0.4011 | 0.030078 | 0.0000195 |
| TFG | PM-Stokes | 8.50 | 0.5611 | 0.4261 | 0.1312 | 0.7167 | 1.544 | 0.3498 | 0.5497 | 1.4070 | 0.111704 | 0.0000449 |

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
