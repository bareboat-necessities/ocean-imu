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

Source commit used for the replay: `de7f2861a7f3c157baac3fec26cfa7274e2bf6fc`.

## Engine speed

Vibration level fixed at 0.60 m/s² of hull broadband
RMS at 2400 rpm, sensor bandwidth 80 Hz.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2890 | 0.5049 | 0.1794 | 0.171 | 0.041 | 0.5879 | 1.00 |
| OU-II | 800 | 0.3751 | 0.4142 | 2.2408 | 1.6530 | 1.597 | 0.800 | 104.0191 | 4.44 |
| OU-II | 1200 | 0.3943 | 0.6314 | 3.7180 | 2.4209 | 2.387 | 0.961 | 103.8107 | 7.36 |
| OU-II | 1600 | 0.3810 | 0.8681 | 5.6343 | 2.8326 | 2.789 | 1.055 | 103.8263 | 11.16 |
| OU-II | 2000 | 0.4106 | 2.0413 | 9.7845 | 3.2798 | 3.133 | 1.875 | 103.8013 | 19.38 |
| OU-II | 2400 | 0.4543 | 3.9001 | 13.2904 | 2.5900 | 2.310 | 3.411 | 103.8087 | 26.32 |
| OU-II | 2800 | 0.5067 | 4.8308 | 14.7896 | 1.3868 | 0.616 | 4.217 | 103.8790 | 29.29 |
| OU-II | 3200 | 0.5657 | 4.1819 | 15.4524 | 2.5478 | 2.183 | 3.784 | 104.0129 | 30.60 |
| OU-III | off | 0.0000 | 0.1792 | 0.4047 | 0.1293 | 0.117 | 0.035 | 0.5135 | 1.00 |
| OU-III | 800 | 0.3751 | 0.7708 | 3.9014 | 2.0738 | 2.017 | 0.841 | 99.2938 | 9.64 |
| OU-III | 1200 | 0.3943 | 1.3530 | 7.8580 | 3.1374 | 3.082 | 1.351 | 91.9877 | 19.42 |
| OU-III | 1600 | 0.3810 | 2.1631 | 12.7704 | 3.6403 | 3.564 | 2.083 | 90.2803 | 31.56 |
| OU-III | 2000 | 0.4106 | 4.3612 | 18.1897 | 3.4325 | 3.220 | 3.962 | 90.6552 | 44.95 |
| OU-III | 2400 | 0.4543 | 4.8645 | 19.7530 | 1.9236 | 1.390 | 4.394 | 90.3038 | 48.81 |
| OU-III | 2800 | 0.5067 | 4.8702 | 19.2358 | 1.8203 | 1.045 | 4.339 | 90.4152 | 47.53 |
| OU-III | 3200 | 0.5657 | 5.7406 | 20.3990 | 3.8852 | 3.396 | 5.068 | 91.1547 | 50.41 |
| TFG | off | 0.0000 | 0.1805 | 0.4083 | 0.0626 | 0.030 | 0.068 | 0.5860 | 1.00 |
| TFG | 800 | 0.3751 | 0.2742 | 0.8969 | 5.4966 | 5.132 | 0.434 | 57.4086 | 2.20 |
| TFG | 1200 | 0.3943 | 0.2627 | 0.9802 | 4.6315 | 3.247 | 0.416 | 61.0425 | 2.40 |
| TFG | 1600 | 0.3810 | 0.2682 | 1.2564 | 6.7973 | 5.567 | 0.501 | 69.3579 | 3.08 |
| TFG | 2000 | 0.4106 | 0.4821 | 10.4477 | 9.8238 | 7.340 | 9.776 | 97.5115 | 25.59 |
| TFG | 2400 | 0.4543 | 1.0014 | 11.9833 | 7.7176 | 6.417 | 7.759 | 104.3151 | 29.35 |
| TFG | 2800 | 0.5067 | 4.0216 | 35.9567 | 4.8360 | 2.492 | 8.177 | 104.2708 | 88.06 |
| TFG | 3200 | 0.5657 | 8.0561 | 63.0917 | 9.1070 | 5.382 | 13.426 | 104.2034 | 154.52 |

## Vibration level

Engine speed fixed at 2400 rpm, sensor bandwidth 80 Hz.  The level is the hull broadband RMS
before the sensor's anti-alias filter, so a quiet, well-isolated
installation sits at the low end and a sensor near the engine bed at the
high end.

| Family | Hull level [m/s²] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2890 | 0.5049 | 0.1794 | 0.171 | 0.041 | 0.5879 | 1.00 |
| OU-II | 0.15 | 0.1136 | 0.2851 | 0.7974 | 0.6132 | 0.514 | 0.281 | 65.1223 | 1.58 |
| OU-II | 0.30 | 0.2272 | 0.4405 | 2.4352 | 1.9778 | 1.863 | 0.625 | 104.3486 | 4.82 |
| OU-II | 1.20 | 0.9087 | 110.2032 | 130.2870 | 5.9467 | 3.832 | 101.909 | 103.8998 | 258.03 |
| OU-II | 2.40 | 1.8173 | 243.8788 | 252.3053 | 19.1508 | 15.786 | 236.636 | 103.7257 | 499.69 |
| OU-III | off | 0.0000 | 0.1792 | 0.4047 | 0.1293 | 0.117 | 0.035 | 0.5135 | 1.00 |
| OU-III | 0.15 | 0.1136 | 0.2609 | 0.9701 | 1.6485 | 1.602 | 0.331 | 49.1930 | 2.40 |
| OU-III | 0.30 | 0.2272 | 0.6418 | 3.5172 | 1.8882 | 1.801 | 0.631 | 90.0930 | 8.69 |
| OU-III | 1.20 | 0.9087 | 83.6152 | 108.1055 | 21.6007 | 18.738 | 79.731 | 104.3148 | 267.13 |
| OU-III | 2.40 | 1.8173 | 171.5063 | 180.8043 | 27.9033 | 25.557 | 165.808 | 103.8806 | 446.76 |
| TFG | off | 0.0000 | 0.1805 | 0.4083 | 0.0626 | 0.030 | 0.068 | 0.5860 | 1.00 |
| TFG | 0.15 | 0.1136 | 0.1947 | 0.5945 | 2.9880 | 2.773 | 0.335 | 9.6723 | 1.46 |
| TFG | 0.30 | 0.2272 | 0.2308 | 1.5029 | 8.6505 | 8.244 | 1.183 | 33.3644 | 3.68 |
| TFG | 1.20 | 0.9087 | 17.9067 | 93.1927 | 18.7018 | 14.473 | 13.469 | 103.3350 | 228.24 |
| TFG | 2.40 | 1.8173 | 6.0312 | 32.6952 | 15.7481 | 13.257 | 4.418 | 103.8784 | 80.07 |

## Sensor anti-alias bandwidth

Engine speed fixed at 2400 rpm and the hull level at
0.60 m/s².  A narrower filter removes power *and*
stops the high crank orders from folding below Nyquist.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2890 | 0.5049 | 0.1794 | 0.171 | 0.041 | 0.5879 | 1.00 |
| OU-II | 20 | 0.1865 | 0.2836 | 0.7478 | 1.4507 | 1.311 | 0.402 | 34.1251 | 1.48 |
| OU-II | 40 | 0.2939 | 0.3965 | 2.2380 | 0.9451 | 0.801 | 0.669 | 104.6272 | 4.43 |
| OU-II | 160 | 0.6034 | 2.4417 | 11.3922 | 1.6713 | 0.988 | 2.123 | 104.0588 | 22.56 |
| OU-III | off | 0.0000 | 0.1792 | 0.4047 | 0.1293 | 0.117 | 0.035 | 0.5135 | 1.00 |
| OU-III | 20 | 0.1865 | 0.3214 | 0.9884 | 1.5109 | 1.418 | 0.460 | 66.7206 | 2.44 |
| OU-III | 40 | 0.2939 | 0.6615 | 3.6447 | 0.8024 | 0.545 | 0.613 | 90.2023 | 9.01 |
| OU-III | 160 | 0.6034 | 5.0246 | 21.0116 | 1.8682 | 1.157 | 4.514 | 90.1981 | 51.92 |
| TFG | off | 0.0000 | 0.1805 | 0.4083 | 0.0626 | 0.030 | 0.068 | 0.5860 | 1.00 |
| TFG | 20 | 0.1865 | 0.2106 | 0.9191 | 7.7141 | 7.276 | 0.671 | 28.0731 | 2.25 |
| TFG | 40 | 0.2939 | 0.3229 | 2.9824 | 13.3501 | 13.070 | 2.544 | 78.4972 | 7.30 |
| TFG | 160 | 0.6034 | 3.3379 | 28.9240 | 8.3931 | 7.664 | 7.010 | 100.9113 | 70.84 |

## Matched-power control

The same bandwidth sweep with the hull level rescaled so every cell
records the same vibration RMS as the nominal cruise cell.  The folded
line frequencies still differ from cell to cell; the recorded power no
longer does.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2890 | 0.5049 | 0.1794 | 0.171 | 0.041 | 0.5879 | 1.00 |
| OU-II | 20 | 0.4543 | 0.2976 | 1.6880 | 1.5148 | 1.404 | 0.775 | 103.9127 | 3.34 |
| OU-II | 40 | 0.4543 | 1.2073 | 6.8226 | 1.0523 | 0.690 | 1.200 | 103.9424 | 13.51 |
| OU-II | 160 | 0.4543 | 0.9499 | 5.2654 | 1.4416 | 1.117 | 0.758 | 103.8939 | 10.43 |
| OU-III | off | 0.0000 | 0.1792 | 0.4047 | 0.1293 | 0.117 | 0.035 | 0.5135 | 1.00 |
| OU-III | 20 | 0.4543 | 0.3900 | 2.3237 | 0.9525 | 0.764 | 0.510 | 104.4087 | 5.74 |
| OU-III | 40 | 0.4543 | 2.3144 | 11.6680 | 1.4241 | 0.864 | 2.108 | 91.1989 | 28.83 |
| OU-III | 160 | 0.4543 | 1.9631 | 10.6535 | 1.7811 | 1.489 | 1.800 | 90.0535 | 26.32 |
| TFG | off | 0.0000 | 0.1805 | 0.4083 | 0.0626 | 0.030 | 0.068 | 0.5860 | 1.00 |
| TFG | 20 | 0.4543 | 0.4997 | 5.3162 | 9.3018 | 8.411 | 4.632 | 104.6862 | 13.02 |
| TFG | 40 | 0.4543 | 1.3872 | 15.2182 | 3.6434 | 2.992 | 7.421 | 103.9531 | 37.27 |
| TFG | 160 | 0.4543 | 0.4524 | 8.8625 | 10.1498 | 8.166 | 8.157 | 97.8811 | 21.71 |

## Sensor attribution

The nominal cruise cell rerun with the model's gyroscope terms switched
off, so the accelerometer is the only perturbed sensor.  Compare against
the 2400 rpm row of the engine-speed table.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2890 | 0.5049 | 0.1794 | 0.171 | 0.041 | 0.5879 | 1.00 |
| OU-II | 2400 | 0.4543 | 3.9010 | 13.2921 | 2.5894 | 2.309 | 3.412 | 103.8059 | 26.32 |
| OU-III | off | 0.0000 | 0.1792 | 0.4047 | 0.1293 | 0.117 | 0.035 | 0.5135 | 1.00 |
| OU-III | 2400 | 0.4543 | 4.8664 | 19.7557 | 1.9233 | 1.389 | 4.395 | 90.3050 | 48.82 |
| TFG | off | 0.0000 | 0.1805 | 0.4083 | 0.0626 | 0.030 | 0.068 | 0.5860 | 1.00 |
| TFG | 2400 | 0.4543 | 1.0029 | 11.9957 | 7.7151 | 6.415 | 7.762 | 104.3201 | 29.38 |

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
