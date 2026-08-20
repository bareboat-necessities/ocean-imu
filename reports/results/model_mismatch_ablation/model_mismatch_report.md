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

Source commit used for the replay: `ef797b4d84bde50c8186e9a4e0bccf30f76a42b2`.

## Pooled RMS across the eight equal-duration records

Because every record contributes the same 900 s window at the same sample
rate, `sqrt(mean(record_RMS^2))` is the exact pooled RMS over their concatenation.

| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | 0.3328 | 0.2241 | 0.2399 | 0.4675 | 20.319 | 0.0900 | 0.2137 | 0.5560 | 0.038405 | 0.0000250 |
| OU-III | 0.2464 | 0.1736 | 0.1401 | 0.3323 | 11.867 | 0.0701 | 0.0647 | 0.2770 | 0.013235 | 0.0000243 |
| TFG | 0.5425 | 0.1993 | 0.1553 | 0.5984 | 13.152 | 0.0888 | 0.1339 | 0.3026 | 0.025373 | 0.0002315 |

## Per-record RMS

| Family | Sea | Hs [m] | X [m] | Y [m] | Z [m] | 3D [m] | Z / Hs [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s²] | Gyro bias 3D [rad/s] |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | JONSWAP | 0.27 | 0.0180 | 0.0121 | 0.0147 | 0.0262 | 5.453 | 0.0297 | 0.0551 | 0.0951 | 0.009612 | 0.0000032 |
| OU-II | JONSWAP | 1.50 | 0.1064 | 0.0636 | 0.0761 | 0.1455 | 5.075 | 0.0453 | 0.1742 | 0.3795 | 0.030200 | 0.0000045 |
| OU-II | JONSWAP | 4.00 | 0.2971 | 0.1949 | 0.2008 | 0.4081 | 5.019 | 0.1694 | 0.3798 | 0.6815 | 0.070801 | 0.0000110 |
| OU-II | JONSWAP | 8.50 | 0.6259 | 0.4317 | 0.4248 | 0.8709 | 4.997 | 0.0988 | 0.4173 | 1.1978 | 0.071982 | 0.0000139 |
| OU-II | PM-Stokes | 0.27 | 0.0169 | 0.0122 | 0.0146 | 0.0254 | 5.410 | 0.0217 | 0.0311 | 0.0236 | 0.002962 | 0.0000031 |
| OU-II | PM-Stokes | 1.50 | 0.0962 | 0.0601 | 0.0757 | 0.1364 | 5.049 | 0.0333 | 0.0474 | 0.2110 | 0.003273 | 0.0000049 |
| OU-II | PM-Stokes | 4.00 | 0.2601 | 0.1735 | 0.2035 | 0.3731 | 5.087 | 0.1229 | 0.0637 | 0.1695 | 0.020442 | 0.0000083 |
| OU-II | PM-Stokes | 8.50 | 0.5631 | 0.3734 | 0.4317 | 0.8018 | 5.079 | 0.0818 | 0.0788 | 0.5889 | 0.012779 | 0.0000674 |
| OU-III | JONSWAP | 0.27 | 0.0127 | 0.0085 | 0.0090 | 0.0177 | 3.329 | 0.0238 | 0.0319 | 0.0346 | 0.005139 | 0.0000025 |
| OU-III | JONSWAP | 1.50 | 0.0794 | 0.0481 | 0.0442 | 0.1028 | 2.949 | 0.0548 | 0.0669 | 0.1102 | 0.013123 | 0.0000037 |
| OU-III | JONSWAP | 4.00 | 0.2295 | 0.1521 | 0.1175 | 0.2993 | 2.936 | 0.1001 | 0.1036 | 0.2757 | 0.023109 | 0.0000106 |
| OU-III | JONSWAP | 8.50 | 0.4866 | 0.3390 | 0.2438 | 0.6412 | 2.869 | 0.0488 | 0.0684 | 0.3270 | 0.007841 | 0.0000095 |
| OU-III | PM-Stokes | 0.27 | 0.0118 | 0.0085 | 0.0092 | 0.0172 | 3.423 | 0.0204 | 0.0300 | 0.0258 | 0.003149 | 0.0000025 |
| OU-III | PM-Stokes | 1.50 | 0.0674 | 0.0438 | 0.0454 | 0.0924 | 3.029 | 0.0323 | 0.0459 | 0.2154 | 0.003453 | 0.0000047 |
| OU-III | PM-Stokes | 4.00 | 0.1832 | 0.1324 | 0.1213 | 0.2565 | 3.032 | 0.1232 | 0.0587 | 0.1553 | 0.020321 | 0.0000073 |
| OU-III | PM-Stokes | 8.50 | 0.3890 | 0.2849 | 0.2548 | 0.5453 | 2.997 | 0.0817 | 0.0782 | 0.5887 | 0.013131 | 0.0000665 |
| TFG | JONSWAP | 0.27 | 0.0308 | 0.0167 | 0.0084 | 0.0360 | 3.113 | 0.0344 | 0.0742 | 0.1228 | 0.012819 | 0.0001661 |
| TFG | JONSWAP | 1.50 | 0.1655 | 0.0617 | 0.0444 | 0.1821 | 2.961 | 0.0555 | 0.0642 | 0.1284 | 0.012764 | 0.0001318 |
| TFG | JONSWAP | 4.00 | 0.5015 | 0.2229 | 0.1261 | 0.5631 | 3.152 | 0.1458 | 0.2657 | 0.3635 | 0.046443 | 0.0001628 |
| TFG | JONSWAP | 8.50 | 0.9696 | 0.3539 | 0.2745 | 1.0681 | 3.230 | 0.1473 | 0.1033 | 0.2631 | 0.025818 | 0.0002879 |
| TFG | PM-Stokes | 0.27 | 0.0293 | 0.0173 | 0.0085 | 0.0351 | 3.155 | 0.0281 | 0.0855 | 0.2106 | 0.012513 | 0.0000262 |
| TFG | PM-Stokes | 1.50 | 0.1579 | 0.0618 | 0.0441 | 0.1752 | 2.940 | 0.0311 | 0.1176 | 0.4277 | 0.022933 | 0.0002983 |
| TFG | PM-Stokes | 4.00 | 0.4736 | 0.2188 | 0.1286 | 0.5373 | 3.214 | 0.0616 | 0.1312 | 0.4008 | 0.021116 | 0.0002280 |
| TFG | PM-Stokes | 8.50 | 0.9403 | 0.2946 | 0.2847 | 1.0256 | 3.350 | 0.1019 | 0.1188 | 0.3343 | 0.029512 | 0.0003641 |

## Interpretation boundary

`--no-noise` removes the simulator's stochastic and calibration-error
sensor terms before they reach the filters.  It does **not** remove wave
nonlinearity, attitude/translation coupling, OU/TFG prior mismatch, the
integral pseudo-measurements, finite adaptation bandwidth, startup residue,
or discretization.  Those effects are intentionally what this ablation
leaves visible.
