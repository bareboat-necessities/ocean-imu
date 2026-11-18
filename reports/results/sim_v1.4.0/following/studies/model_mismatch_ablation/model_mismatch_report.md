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

Source commit used for the replay: `1ad97346644532cb44c5f32f29aaa3ec148140ff`.

## Pooled RMS across the eight equal-duration records

Because every record contributes the same 900 s window at the same sample
rate, `sqrt(mean(record_RMS^2))` is the exact pooled RMS over their concatenation.

| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | 0.3650 | 0.2280 | 0.2802 | 0.5135 | 22.995 | 0.0296 | 0.0235 | 0.0579 | 0.004603 | 0.0000219 |
| OU-III | 0.2981 | 0.1743 | 0.1649 | 0.3827 | 13.533 | 0.0648 | 0.0374 | 0.1452 | 0.012024 | 0.0000233 |
| TFG | 0.1344 | 0.1134 | 0.0661 | 0.1879 | 5.426 | 0.1537 | 0.2555 | 0.6467 | 0.051188 | 0.0000197 |

## Per-record RMS

| Family | Sea | Hs [m] | X [m] | Y [m] | Z [m] | 3D [m] | Z / Hs [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | JONSWAP | 0.27 | 0.0127 | 0.0047 | 0.0179 | 0.0225 | 6.641 | 0.0434 | 0.0072 | 0.0502 | 0.007466 | 0.0000115 |
| OU-II | JONSWAP | 1.50 | 0.0881 | 0.0399 | 0.0925 | 0.1338 | 6.164 | 0.0087 | 0.0145 | 0.0880 | 0.001762 | 0.0000218 |
| OU-II | JONSWAP | 4.00 | 0.2702 | 0.1508 | 0.2322 | 0.3869 | 5.804 | 0.0250 | 0.0231 | 0.0468 | 0.003872 | 0.0000211 |
| OU-II | JONSWAP | 8.50 | 0.6804 | 0.4368 | 0.5055 | 0.9536 | 5.947 | 0.0249 | 0.0378 | 0.0517 | 0.001229 | 0.0000307 |
| OU-II | PM-Stokes | 0.27 | 0.0120 | 0.0046 | 0.0171 | 0.0214 | 6.336 | 0.0487 | 0.0077 | 0.0542 | 0.008377 | 0.0000116 |
| OU-II | PM-Stokes | 1.50 | 0.0848 | 0.0395 | 0.0913 | 0.1308 | 6.088 | 0.0088 | 0.0147 | 0.0365 | 0.001915 | 0.0000159 |
| OU-II | PM-Stokes | 4.00 | 0.2656 | 0.1487 | 0.2302 | 0.3816 | 5.755 | 0.0275 | 0.0227 | 0.0761 | 0.004328 | 0.0000234 |
| OU-II | PM-Stokes | 8.50 | 0.6664 | 0.4208 | 0.4980 | 0.9323 | 5.859 | 0.0241 | 0.0376 | 0.0406 | 0.001237 | 0.0000298 |
| OU-III | JONSWAP | 0.27 | 0.0259 | 0.0074 | 0.0209 | 0.0342 | 7.758 | 0.0891 | 0.0218 | 0.1050 | 0.015579 | 0.0000152 |
| OU-III | JONSWAP | 1.50 | 0.0991 | 0.0370 | 0.0674 | 0.1254 | 4.492 | 0.0319 | 0.0273 | 0.1240 | 0.006357 | 0.0000233 |
| OU-III | JONSWAP | 4.00 | 0.2281 | 0.1247 | 0.1381 | 0.2943 | 3.452 | 0.0686 | 0.0362 | 0.0885 | 0.012560 | 0.0000213 |
| OU-III | JONSWAP | 8.50 | 0.5436 | 0.3232 | 0.2894 | 0.6955 | 3.405 | 0.0185 | 0.0504 | 0.2003 | 0.006277 | 0.0000300 |
| OU-III | PM-Stokes | 0.27 | 0.0243 | 0.0071 | 0.0197 | 0.0321 | 7.310 | 0.1167 | 0.0232 | 0.1199 | 0.020272 | 0.0000148 |
| OU-III | PM-Stokes | 1.50 | 0.0915 | 0.0357 | 0.0652 | 0.1179 | 4.346 | 0.0234 | 0.0255 | 0.1223 | 0.005062 | 0.0000240 |
| OU-III | PM-Stokes | 4.00 | 0.2260 | 0.1215 | 0.1376 | 0.2912 | 3.441 | 0.0716 | 0.0403 | 0.1064 | 0.013402 | 0.0000231 |
| OU-III | PM-Stokes | 8.50 | 0.5413 | 0.3250 | 0.2934 | 0.6962 | 3.452 | 0.0188 | 0.0574 | 0.2316 | 0.007753 | 0.0000299 |
| TFG | JONSWAP | 0.27 | 0.0066 | 0.0017 | 0.0066 | 0.0095 | 2.449 | 0.0323 | 0.0119 | 0.0365 | 0.005880 | 0.0000053 |
| TFG | JONSWAP | 1.50 | 0.0263 | 0.0126 | 0.0168 | 0.0336 | 1.120 | 0.0241 | 0.0170 | 0.0441 | 0.005113 | 0.0000069 |
| TFG | JONSWAP | 4.00 | 0.0897 | 0.0736 | 0.0459 | 0.1248 | 1.148 | 0.0154 | 0.1079 | 0.2489 | 0.018620 | 0.0000108 |
| TFG | JONSWAP | 8.50 | 0.2254 | 0.1908 | 0.1095 | 0.3150 | 1.289 | 0.2561 | 0.4501 | 1.1330 | 0.088924 | 0.0000327 |
| TFG | PM-Stokes | 0.27 | 0.0079 | 0.0021 | 0.0079 | 0.0114 | 2.928 | 0.0368 | 0.0137 | 0.0423 | 0.006697 | 0.0000059 |
| TFG | PM-Stokes | 1.50 | 0.0315 | 0.0158 | 0.0200 | 0.0405 | 1.332 | 0.0127 | 0.0143 | 0.0208 | 0.003145 | 0.0000079 |
| TFG | PM-Stokes | 4.00 | 0.1066 | 0.0751 | 0.0538 | 0.1411 | 1.345 | 0.0292 | 0.1472 | 0.3440 | 0.025764 | 0.0000129 |
| TFG | PM-Stokes | 8.50 | 0.2693 | 0.2346 | 0.1311 | 0.3804 | 1.542 | 0.3453 | 0.5345 | 1.3697 | 0.109215 | 0.0000398 |

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
