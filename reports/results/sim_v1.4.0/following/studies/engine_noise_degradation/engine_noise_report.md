# Engine-noise degradation study

How far the deployed estimators degrade when the IMU also records the
vibration of an inboard auxiliary diesel.  The eight stationary JONSWAP /
PM-Stokes records are replayed through OU-II, OU-III, and TFG with the
ordinary sensor noise models on and the engine vibration model added on
top; the filters keep their deployed covariances, adaptation, startup, and
regularization throughout.

The vessel modeled is a mid-size recreational cruising sailboat: a
naturally aspirated three-cylinder four-stroke diesel on flexible mounts,
a 2.6:1 reduction gear, and a three-blade fixed propeller.  The model is
a sensor-path model.  It adds crank orders, driveline shaft- and
blade-rate lines, a broadband structural floor, governor hunting, mount
transmissibility, the sensor's finite anti-alias bandwidth, accelerometer
vibration rectification, and gyroscope g-sensitivity.  It does not change
the vessel's rigid-body response to the sea.

Scoring uses the trailing **900 s** of each 1200 s record, and
every reported number pools the eight equal-duration records as
`sqrt(mean(record_RMS^2))`, which is the exact RMS over their concatenation.

Source commit used for the replay: `1ad97346644532cb44c5f32f29aaa3ec148140ff`.

## Engine speed

Vibration level fixed at 0.60 m/s² of hull broadband
RMS at 2400 rpm, sensor bandwidth 80 Hz.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3238 | 0.5762 | 0.1469 | 0.137 | 0.044 | 0.8287 | 1.00 |
| OU-II | 800 | 0.3751 | 1.5020 | 8.0979 | 1.9853 | 1.780 | 1.527 | 104.0786 | 14.06 |
| OU-II | 1200 | 0.3943 | 2.5533 | 12.7102 | 3.0092 | 2.773 | 2.477 | 104.0045 | 22.06 |
| OU-II | 1600 | 0.3810 | 4.9140 | 20.5497 | 3.7245 | 3.255 | 4.496 | 103.9287 | 35.67 |
| OU-II | 2000 | 0.4106 | 8.8925 | 36.4225 | 4.1999 | 3.028 | 8.403 | 103.7770 | 63.22 |
| OU-II | 2400 | 0.4543 | 12.3323 | 47.6041 | 3.7320 | 1.399 | 10.030 | 104.0896 | 82.62 |
| OU-II | 2800 | 0.5067 | 13.0866 | 49.0311 | 4.7224 | 2.961 | 11.189 | 103.8136 | 85.10 |
| OU-II | 3200 | 0.5657 | 18.3988 | 60.8936 | 8.3976 | 7.192 | 17.220 | 104.0100 | 105.69 |
| OU-III | off | 0.0000 | 0.2092 | 0.4853 | 0.1047 | 0.091 | 0.042 | 0.6607 | 1.00 |
| OU-III | 800 | 0.3751 | 2.4890 | 13.2823 | 2.0185 | 1.521 | 2.136 | 90.2977 | 27.37 |
| OU-III | 1200 | 0.3943 | 2.9588 | 18.9834 | 2.7663 | 2.109 | 2.085 | 90.8863 | 39.12 |
| OU-III | 1600 | 0.3810 | 4.4446 | 29.7795 | 2.8406 | 1.773 | 3.329 | 90.0627 | 61.36 |
| OU-III | 2000 | 0.4106 | 4.3331 | 32.1499 | 2.7757 | 1.451 | 2.551 | 89.9913 | 66.25 |
| OU-III | 2400 | 0.4543 | 100.5112 | 135.8805 | 15.9697 | 6.307 | 80.987 | 90.3402 | 279.98 |
| OU-III | 2800 | 0.5067 | 115.1893 | 145.3641 | 16.2713 | 11.824 | 97.622 | 90.1486 | 299.52 |
| OU-III | 3200 | 0.5657 | 75.3975 | 114.7377 | 15.6980 | 11.957 | 66.763 | 90.5881 | 236.42 |
| TFG | off | 0.0000 | 0.1469 | 0.2914 | 0.1225 | 0.113 | 0.032 | 0.6787 | 1.00 |
| TFG | 800 | 0.3751 | 0.7757 | 12.5431 | 7.6589 | 7.176 | 9.929 | 103.8611 | 43.04 |
| TFG | 1200 | 0.3943 | 2.6340 | 22.5522 | 6.6596 | 5.931 | 9.436 | 103.5253 | 77.39 |
| TFG | 1600 | 0.3810 | 1.8042 | 16.6716 | 9.6899 | 8.572 | 8.568 | 95.9615 | 57.21 |
| TFG | 2000 | 0.4106 | 4.3571 | 33.9589 | 9.0446 | 6.893 | 8.461 | 96.1295 | 116.53 |
| TFG | 2400 | 0.4543 | 12.4608 | 60.6691 | 8.8330 | 6.013 | 9.750 | 94.9436 | 208.18 |
| TFG | 2800 | 0.5067 | 19.1303 | 83.3627 | 12.3778 | 8.051 | 13.239 | 104.1174 | 286.05 |
| TFG | 3200 | 0.5657 | 23.6103 | 96.4026 | 13.4223 | 7.555 | 12.310 | 105.2118 | 330.80 |

## Vibration level

Engine speed fixed at 2400 rpm, sensor bandwidth 80 Hz.  The level is the hull broadband RMS
before the sensor's anti-alias filter, so a quiet, well-isolated
installation sits at the low end and a sensor near the engine bed at the
high end.

| Family | Hull level [m/s²] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3238 | 0.5762 | 0.1469 | 0.137 | 0.044 | 0.8287 | 1.00 |
| OU-II | 0.15 | 0.1136 | 0.3582 | 1.7356 | 1.0784 | 0.943 | 0.926 | 90.4366 | 3.01 |
| OU-II | 0.30 | 0.2272 | 1.2886 | 6.4580 | 1.5207 | 1.046 | 1.060 | 90.0850 | 11.21 |
| OU-II | 1.20 | 0.9087 | 263.4175 | 274.3070 | 22.9899 | 9.184 | 248.506 | 103.3161 | 476.10 |
| OU-II | 2.40 | 1.8173 | 545.7655 | 556.4335 | 28.8071 | 17.436 | 517.975 | 104.3908 | 965.77 |
| OU-III | off | 0.0000 | 0.2092 | 0.4853 | 0.1047 | 0.091 | 0.042 | 0.6607 | 1.00 |
| OU-III | 0.15 | 0.1136 | 0.3986 | 1.4401 | 1.6609 | 1.622 | 0.658 | 44.3625 | 2.97 |
| OU-III | 0.30 | 0.2272 | 0.9216 | 5.2749 | 1.0870 | 0.475 | 0.527 | 90.1098 | 10.87 |
| OU-III | 1.20 | 0.9087 | 167.9350 | 193.8082 | 19.9611 | 6.972 | 148.339 | 91.5097 | 399.34 |
| OU-III | 2.40 | 1.8173 | 312.7046 | 329.0288 | 29.9467 | 19.554 | 282.460 | 104.1998 | 677.97 |
| TFG | off | 0.0000 | 0.1469 | 0.2914 | 0.1225 | 0.113 | 0.032 | 0.6787 | 1.00 |
| TFG | 0.15 | 0.1136 | 0.1661 | 0.6540 | 3.0112 | 2.853 | 0.366 | 8.7262 | 2.24 |
| TFG | 0.30 | 0.2272 | 0.3890 | 2.1527 | 9.3714 | 9.140 | 1.384 | 35.6756 | 7.39 |
| TFG | 1.20 | 0.9087 | 23.3937 | 99.3654 | 26.1798 | 21.688 | 11.564 | 103.3844 | 340.97 |
| TFG | 2.40 | 1.8173 | 3.9578 | 7.0865 | 16.9943 | 15.906 | 4.638 | 103.8428 | 24.32 |

## Sensor anti-alias bandwidth

Engine speed fixed at 2400 rpm and the hull level at
0.60 m/s².  A narrower filter removes power *and*
stops the high crank orders from folding below Nyquist.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3238 | 0.5762 | 0.1469 | 0.137 | 0.044 | 0.8287 | 1.00 |
| OU-II | 20 | 0.1865 | 0.3630 | 2.0281 | 1.2410 | 1.127 | 1.202 | 96.1422 | 3.52 |
| OU-II | 40 | 0.2939 | 1.4827 | 7.0472 | 1.9625 | 1.552 | 1.289 | 91.3417 | 12.23 |
| OU-II | 160 | 0.6034 | 14.5922 | 54.3397 | 3.8841 | 0.917 | 10.409 | 100.0425 | 94.31 |
| OU-III | off | 0.0000 | 0.2092 | 0.4853 | 0.1047 | 0.091 | 0.042 | 0.6607 | 1.00 |
| OU-III | 20 | 0.1865 | 0.4203 | 2.3105 | 0.4963 | 0.207 | 0.849 | 89.8796 | 4.76 |
| OU-III | 40 | 0.2939 | 1.4068 | 8.6460 | 1.7522 | 1.266 | 0.839 | 90.1676 | 17.82 |
| OU-III | 160 | 0.6034 | 4.9869 | 34.6207 | 2.2902 | 0.217 | 3.311 | 89.9225 | 71.34 |
| TFG | off | 0.0000 | 0.1469 | 0.2914 | 0.1225 | 0.113 | 0.032 | 0.6787 | 1.00 |
| TFG | 20 | 0.1865 | 0.3377 | 2.0286 | 9.9071 | 9.584 | 1.279 | 35.9560 | 6.96 |
| TFG | 40 | 0.2939 | 3.1264 | 26.7795 | 9.0464 | 8.486 | 7.861 | 87.7109 | 91.89 |
| TFG | 160 | 0.6034 | 5.0758 | 39.1514 | 8.8581 | 7.769 | 9.247 | 86.2705 | 134.35 |

## Matched-power control

The same bandwidth sweep with the hull level rescaled so every cell
records the same vibration RMS as the nominal cruise cell.  The folded
line frequencies still differ from cell to cell; the recorded power no
longer does.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3238 | 0.5762 | 0.1469 | 0.137 | 0.044 | 0.8287 | 1.00 |
| OU-II | 20 | 0.4543 | 1.0950 | 5.8302 | 3.2135 | 3.010 | 1.092 | 104.0368 | 10.12 |
| OU-II | 40 | 0.4543 | 5.0173 | 24.3085 | 4.2965 | 3.435 | 5.384 | 104.0416 | 42.19 |
| OU-II | 160 | 0.4543 | 4.1626 | 21.8847 | 2.2540 | 0.501 | 3.640 | 90.7150 | 37.98 |
| OU-III | off | 0.0000 | 0.2092 | 0.4853 | 0.1047 | 0.091 | 0.042 | 0.6607 | 1.00 |
| OU-III | 20 | 0.4543 | 2.0329 | 14.3632 | 2.2181 | 1.348 | 1.033 | 90.7072 | 29.60 |
| OU-III | 40 | 0.4543 | 4.1578 | 29.7442 | 2.9601 | 1.862 | 2.198 | 90.3151 | 61.29 |
| OU-III | 160 | 0.4543 | 2.2848 | 13.4111 | 1.3827 | 0.217 | 1.291 | 89.8349 | 27.63 |
| TFG | off | 0.0000 | 0.1469 | 0.2914 | 0.1225 | 0.113 | 0.032 | 0.6787 | 1.00 |
| TFG | 20 | 0.4543 | 16.1332 | 75.1873 | 12.2545 | 8.281 | 9.195 | 103.9164 | 258.00 |
| TFG | 40 | 0.4543 | 10.4870 | 60.9158 | 9.9402 | 5.737 | 6.668 | 105.4900 | 209.03 |
| TFG | 160 | 0.4543 | 2.0284 | 18.5226 | 9.0651 | 8.761 | 7.180 | 69.5928 | 63.56 |

## Sensor attribution

The nominal cruise cell rerun with the model's gyroscope terms switched
off, so the accelerometer is the only perturbed sensor.  Compare against
the 2400 rpm row of the engine-speed table.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3238 | 0.5762 | 0.1469 | 0.137 | 0.044 | 0.8287 | 1.00 |
| OU-II | 2400 | 0.4543 | 12.6242 | 48.1380 | 3.7558 | 1.387 | 10.272 | 104.2115 | 83.55 |
| OU-III | off | 0.0000 | 0.2092 | 0.4853 | 0.1047 | 0.091 | 0.042 | 0.6607 | 1.00 |
| OU-III | 2400 | 0.4543 | 112.6954 | 147.3693 | 18.7298 | 11.402 | 95.713 | 89.3247 | 303.66 |
| TFG | off | 0.0000 | 0.1469 | 0.2914 | 0.1225 | 0.113 | 0.032 | 0.6787 | 1.00 |
| TFG | 2400 | 0.4543 | 8.5235 | 54.9839 | 8.8688 | 6.284 | 10.924 | 94.9431 | 188.67 |

## Figures

- `ou_engine_noise_speed.svg`: recorded vibration and the displacement,
  pitch, and yaw response against engine speed.
- `ou_engine_noise_mechanism.svg`: the level sweep, the bandwidth sweep with
  its matched-power control, every cell of the study collapsed onto
  recorded vibration RMS, and the rectified tilt offset that drives
  all of it.

Both are mirrored byte-for-byte into `doc/kalman_ou_iii/` for the article.

## Interpretation boundary

The engine model perturbs the accelerometer and gyroscope only.  The
magnetometer, the wave records, and the vessel's rigid-body motion are
unchanged, so this study bounds the sensor-path cost of motoring and not
the full difference between sailing and motoring.  A real passage under
power also changes the encounter spectrum, adds propeller-induced surge,
and runs the engine at a speed that itself varies with the sea; none of
that is modeled here.
