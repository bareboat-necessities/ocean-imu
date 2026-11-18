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

Source commit used for the replay: `99592c4078a5ba45c9a9c6edd455782ed440e495`.

## Pooled RMS across the eight equal-duration records

Because every record contributes the same 900 s window at the same sample
rate, `sqrt(mean(record_RMS^2))` is the exact pooled RMS over their concatenation.

| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | 0.2713 | 0.1889 | 0.2347 | 0.4054 | 18.693 | 0.0240 | 0.0204 | 0.0999 | 0.003828 | 0.0000239 |
| OU-III | 0.2054 | 0.1756 | 0.1362 | 0.3026 | 10.846 | 0.0612 | 0.0614 | 0.2646 | 0.014231 | 0.0000269 |
| TFG | 0.1448 | 0.1567 | 0.0775 | 0.2270 | 6.176 | 0.2420 | 0.1070 | 0.1737 | 0.044785 | 0.0000134 |

## Per-record RMS

| Family | Sea | Hs [m] | X [m] | Y [m] | Z [m] | 3D [m] | Z / Hs [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | JONSWAP | 0.27 | 0.0022 | 0.0013 | 0.0093 | 0.0096 | 3.438 | 0.0093 | 0.0111 | 0.0597 | 0.002472 | 0.0000113 |
| OU-II | JONSWAP | 1.50 | 0.0408 | 0.0224 | 0.0708 | 0.0847 | 4.718 | 0.0126 | 0.0099 | 0.0157 | 0.002507 | 0.0000079 |
| OU-II | JONSWAP | 4.00 | 0.1970 | 0.1051 | 0.2001 | 0.2999 | 5.003 | 0.0115 | 0.0223 | 0.0566 | 0.002725 | 0.0000183 |
| OU-II | JONSWAP | 8.50 | 0.5132 | 0.3628 | 0.4227 | 0.7574 | 4.973 | 0.0349 | 0.0297 | 0.1278 | 0.004396 | 0.0000336 |
| OU-II | PM-Stokes | 0.27 | 0.0023 | 0.0014 | 0.0094 | 0.0098 | 3.481 | 0.0120 | 0.0099 | 0.0595 | 0.002653 | 0.0000113 |
| OU-II | PM-Stokes | 1.50 | 0.0387 | 0.0243 | 0.0685 | 0.0824 | 4.570 | 0.0181 | 0.0117 | 0.0223 | 0.003505 | 0.0000085 |
| OU-II | PM-Stokes | 4.00 | 0.1894 | 0.1056 | 0.1971 | 0.2931 | 4.928 | 0.0117 | 0.0237 | 0.0669 | 0.003051 | 0.0000196 |
| OU-II | PM-Stokes | 8.50 | 0.4975 | 0.3611 | 0.4160 | 0.7423 | 4.894 | 0.0490 | 0.0305 | 0.2191 | 0.007031 | 0.0000482 |
| OU-III | JONSWAP | 0.27 | 0.0019 | 0.0011 | 0.0044 | 0.0049 | 1.629 | 0.0160 | 0.0196 | 0.0698 | 0.004317 | 0.0000106 |
| OU-III | JONSWAP | 1.50 | 0.0285 | 0.0219 | 0.0372 | 0.0517 | 2.483 | 0.0104 | 0.0197 | 0.1051 | 0.003568 | 0.0000117 |
| OU-III | JONSWAP | 4.00 | 0.1434 | 0.0875 | 0.1080 | 0.1997 | 2.699 | 0.0254 | 0.0766 | 0.1672 | 0.013422 | 0.0000205 |
| OU-III | JONSWAP | 8.50 | 0.3750 | 0.3060 | 0.2375 | 0.5391 | 2.794 | 0.0930 | 0.0855 | 0.4060 | 0.020567 | 0.0000355 |
| OU-III | PM-Stokes | 0.27 | 0.0021 | 0.0012 | 0.0047 | 0.0053 | 1.730 | 0.0199 | 0.0191 | 0.0710 | 0.004703 | 0.0000108 |
| OU-III | PM-Stokes | 1.50 | 0.0297 | 0.0276 | 0.0385 | 0.0559 | 2.563 | 0.0110 | 0.0206 | 0.1072 | 0.003651 | 0.0000125 |
| OU-III | PM-Stokes | 4.00 | 0.1483 | 0.1064 | 0.1133 | 0.2148 | 2.833 | 0.0219 | 0.0804 | 0.1686 | 0.013775 | 0.0000235 |
| OU-III | PM-Stokes | 8.50 | 0.3906 | 0.3645 | 0.2540 | 0.5915 | 2.988 | 0.1388 | 0.0945 | 0.5534 | 0.027577 | 0.0000550 |
| TFG | JONSWAP | 0.27 | 0.0008 | 0.0011 | 0.0016 | 0.0021 | 0.607 | 0.0345 | 0.0559 | 0.1357 | 0.011268 | 0.0000031 |
| TFG | JONSWAP | 1.50 | 0.0169 | 0.0166 | 0.0159 | 0.0286 | 1.063 | 0.0730 | 0.0397 | 0.1237 | 0.014125 | 0.0000061 |
| TFG | JONSWAP | 4.00 | 0.0838 | 0.0745 | 0.0542 | 0.1245 | 1.354 | 0.0670 | 0.1072 | 0.2090 | 0.021268 | 0.0000113 |
| TFG | JONSWAP | 8.50 | 0.2394 | 0.2568 | 0.1249 | 0.3726 | 1.469 | 0.3952 | 0.1313 | 0.1673 | 0.070325 | 0.0000176 |
| TFG | PM-Stokes | 0.27 | 0.0009 | 0.0011 | 0.0017 | 0.0022 | 0.630 | 0.0353 | 0.0547 | 0.1335 | 0.011169 | 0.0000031 |
| TFG | PM-Stokes | 1.50 | 0.0209 | 0.0214 | 0.0201 | 0.0360 | 1.337 | 0.0659 | 0.0512 | 0.1460 | 0.014160 | 0.0000083 |
| TFG | PM-Stokes | 4.00 | 0.1030 | 0.1001 | 0.0662 | 0.1581 | 1.655 | 0.0972 | 0.1122 | 0.2043 | 0.024846 | 0.0000154 |
| TFG | PM-Stokes | 8.50 | 0.3034 | 0.3378 | 0.1566 | 0.4803 | 1.842 | 0.5350 | 0.2000 | 0.2350 | 0.096843 | 0.0000252 |

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
