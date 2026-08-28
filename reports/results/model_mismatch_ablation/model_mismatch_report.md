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

Source commit used for the replay: `ca599fe5c35f1d7a462baec0fcacbfc3da70f762`.

## Pooled RMS across the eight equal-duration records

Because every record contributes the same 900 s window at the same sample
rate, `sqrt(mean(record_RMS^2))` is the exact pooled RMS over their concatenation.

| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | 0.3578 | 0.2478 | 0.2372 | 0.4957 | 20.088 | 0.1147 | 0.2277 | 0.5441 | 0.042484 | 0.0000112 |
| OU-III | 0.2772 | 0.2003 | 0.1392 | 0.3692 | 11.790 | 0.1151 | 0.0733 | 0.2609 | 0.021110 | 0.0000097 |
| TFG | 0.5388 | 0.2043 | 0.1494 | 0.5952 | 12.650 | 0.1179 | 0.1839 | 0.4806 | 0.031267 | 0.0002692 |

## Per-record RMS

| Family | Sea | Hs [m] | X [m] | Y [m] | Z [m] | 3D [m] | Z / Hs [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | JONSWAP | 0.27 | 0.0198 | 0.0133 | 0.0147 | 0.0280 | 5.435 | 0.0295 | 0.0510 | 0.0815 | 0.008644 | 0.0000039 |
| OU-II | JONSWAP | 1.50 | 0.1159 | 0.0695 | 0.0758 | 0.1550 | 5.055 | 0.0464 | 0.1706 | 0.3698 | 0.029530 | 0.0000056 |
| OU-II | JONSWAP | 4.00 | 0.3187 | 0.2112 | 0.1977 | 0.4304 | 4.942 | 0.1603 | 0.3869 | 0.6912 | 0.071358 | 0.0000120 |
| OU-II | JONSWAP | 8.50 | 0.6713 | 0.4744 | 0.4172 | 0.9218 | 4.908 | 0.1862 | 0.4647 | 1.1834 | 0.084724 | 0.0000225 |
| OU-II | PM-Stokes | 0.27 | 0.0187 | 0.0135 | 0.0146 | 0.0273 | 5.391 | 0.0228 | 0.0318 | 0.0216 | 0.002616 | 0.0000039 |
| OU-II | PM-Stokes | 1.50 | 0.1062 | 0.0667 | 0.0757 | 0.1465 | 5.046 | 0.0343 | 0.0494 | 0.2113 | 0.003250 | 0.0000060 |
| OU-II | PM-Stokes | 4.00 | 0.2857 | 0.1925 | 0.2029 | 0.3998 | 5.072 | 0.1228 | 0.0662 | 0.1706 | 0.020422 | 0.0000096 |
| OU-II | PM-Stokes | 8.50 | 0.6040 | 0.4184 | 0.4288 | 0.8507 | 5.044 | 0.1581 | 0.0978 | 0.5224 | 0.028069 | 0.0000128 |
| OU-III | JONSWAP | 0.27 | 0.0147 | 0.0098 | 0.0089 | 0.0198 | 3.307 | 0.0239 | 0.0311 | 0.0294 | 0.004745 | 0.0000032 |
| OU-III | JONSWAP | 1.50 | 0.0903 | 0.0548 | 0.0439 | 0.1144 | 2.924 | 0.0552 | 0.0672 | 0.1085 | 0.013019 | 0.0000050 |
| OU-III | JONSWAP | 4.00 | 0.2532 | 0.1709 | 0.1147 | 0.3263 | 2.868 | 0.0899 | 0.1048 | 0.2746 | 0.021845 | 0.0000097 |
| OU-III | JONSWAP | 8.50 | 0.5349 | 0.3832 | 0.2408 | 0.7007 | 2.833 | 0.2301 | 0.0969 | 0.3061 | 0.040577 | 0.0000173 |
| OU-III | PM-Stokes | 0.27 | 0.0140 | 0.0100 | 0.0092 | 0.0195 | 3.406 | 0.0211 | 0.0303 | 0.0236 | 0.002935 | 0.0000033 |
| OU-III | PM-Stokes | 1.50 | 0.0800 | 0.0518 | 0.0452 | 0.1055 | 3.016 | 0.0330 | 0.0475 | 0.2162 | 0.003472 | 0.0000059 |
| OU-III | PM-Stokes | 4.00 | 0.2167 | 0.1580 | 0.1209 | 0.2942 | 3.023 | 0.1240 | 0.0602 | 0.1473 | 0.020391 | 0.0000098 |
| OU-III | PM-Stokes | 8.50 | 0.4499 | 0.3377 | 0.2552 | 0.6177 | 3.002 | 0.1566 | 0.1016 | 0.5419 | 0.028511 | 0.0000135 |
| TFG | JONSWAP | 0.27 | 0.0341 | 0.0192 | 0.0095 | 0.0403 | 3.535 | 0.0275 | 0.0733 | 0.1580 | 0.008596 | 0.0000408 |
| TFG | JONSWAP | 1.50 | 0.1758 | 0.0679 | 0.0463 | 0.1941 | 3.087 | 0.0457 | 0.1202 | 0.2352 | 0.017449 | 0.0000838 |
| TFG | JONSWAP | 4.00 | 0.4950 | 0.2220 | 0.1209 | 0.5558 | 3.022 | 0.1303 | 0.2846 | 0.3763 | 0.048060 | 0.0001348 |
| TFG | JONSWAP | 8.50 | 0.8895 | 0.3509 | 0.2586 | 0.9906 | 3.042 | 0.2425 | 0.1537 | 0.4064 | 0.044365 | 0.0000318 |
| TFG | PM-Stokes | 0.27 | 0.0335 | 0.0202 | 0.0098 | 0.0403 | 3.615 | 0.0366 | 0.1182 | 0.2866 | 0.016740 | 0.0001964 |
| TFG | PM-Stokes | 1.50 | 0.1737 | 0.0683 | 0.0480 | 0.1927 | 3.201 | 0.0293 | 0.0758 | 0.2689 | 0.010552 | 0.0003319 |
| TFG | PM-Stokes | 4.00 | 0.4908 | 0.2323 | 0.1316 | 0.5587 | 3.289 | 0.1037 | 0.2002 | 0.5708 | 0.039416 | 0.0002683 |
| TFG | PM-Stokes | 8.50 | 0.9908 | 0.3119 | 0.2739 | 1.0743 | 3.223 | 0.1404 | 0.2936 | 0.9903 | 0.034922 | 0.0005755 |

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
