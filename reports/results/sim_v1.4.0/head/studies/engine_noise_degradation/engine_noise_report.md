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

Source commit used for the replay: `37d4ddb7a2d8431f91d2807c526f6adc4e22b1f8`.

## Engine speed

Vibration level fixed at 0.60 m/s² of hull broadband
RMS at 2400 rpm, sensor bandwidth 80 Hz.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2733 | 0.4726 | 0.1657 | 0.158 | 0.034 | 0.8199 | 1.00 |
| OU-II | 800 | 0.3751 | 0.4035 | 1.8861 | 1.2337 | 1.188 | 1.250 | 87.0962 | 3.99 |
| OU-II | 1200 | 0.3943 | 0.5075 | 2.4966 | 1.6564 | 1.617 | 1.496 | 94.8596 | 5.28 |
| OU-II | 1600 | 0.3810 | 0.5990 | 2.9925 | 1.5692 | 1.523 | 0.951 | 106.2838 | 6.33 |
| OU-II | 2000 | 0.4106 | 0.9431 | 4.7709 | 1.7374 | 1.633 | 1.035 | 104.3051 | 10.10 |
| OU-II | 2400 | 0.4543 | 1.4308 | 6.9523 | 1.4810 | 1.229 | 1.502 | 104.1042 | 14.71 |
| OU-II | 2800 | 0.5067 | 1.8040 | 8.0242 | 1.3493 | 0.979 | 1.823 | 104.2513 | 16.98 |
| OU-II | 3200 | 0.5657 | 1.8266 | 8.1309 | 1.3283 | 1.131 | 1.847 | 104.0472 | 17.21 |
| OU-III | off | 0.0000 | 0.1687 | 0.3751 | 0.1713 | 0.161 | 0.020 | 0.8858 | 1.00 |
| OU-III | 800 | 0.3751 | 0.7204 | 2.6143 | 0.7056 | 0.545 | 0.954 | 105.8803 | 6.97 |
| OU-III | 1200 | 0.3943 | 0.8397 | 3.3121 | 1.3227 | 1.256 | 0.993 | 105.4609 | 8.83 |
| OU-III | 1600 | 0.3810 | 1.1792 | 5.6553 | 1.9930 | 1.940 | 1.224 | 94.0995 | 15.08 |
| OU-III | 2000 | 0.4106 | 1.3989 | 7.7427 | 1.5596 | 1.270 | 1.221 | 104.2103 | 20.64 |
| OU-III | 2400 | 0.4543 | 1.4775 | 8.5337 | 1.2031 | 0.609 | 1.252 | 98.9084 | 22.75 |
| OU-III | 2800 | 0.5067 | 1.7105 | 9.0675 | 1.0844 | 0.429 | 1.503 | 103.9999 | 24.17 |
| OU-III | 3200 | 0.5657 | 2.1622 | 10.1786 | 1.5082 | 1.229 | 2.085 | 103.8178 | 27.14 |
| TFG | off | 0.0000 | 0.1241 | 0.3597 | 0.1645 | 0.143 | 0.030 | 0.7567 | 1.00 |
| TFG | 800 | 0.3751 | 0.3524 | 3.0469 | 11.1738 | 10.855 | 1.942 | 68.9281 | 8.47 |
| TFG | 1200 | 0.3943 | 0.8116 | 9.9179 | 12.4289 | 12.068 | 5.737 | 99.0836 | 27.57 |
| TFG | 1600 | 0.3810 | 2.1911 | 21.9906 | 6.0169 | 4.180 | 7.054 | 109.9178 | 61.13 |
| TFG | 2000 | 0.4106 | 6.0518 | 50.3182 | 9.5265 | 6.659 | 7.390 | 105.3002 | 139.88 |
| TFG | 2400 | 0.4543 | 5.8983 | 51.2540 | 11.2648 | 8.401 | 8.458 | 106.3177 | 142.48 |
| TFG | 2800 | 0.5067 | 4.2354 | 32.6199 | 7.1318 | 6.110 | 6.175 | 104.0268 | 90.68 |
| TFG | 3200 | 0.5657 | 1.0069 | 13.4342 | 12.2114 | 11.438 | 7.800 | 108.1066 | 37.35 |

## Vibration level

Engine speed fixed at 2400 rpm, sensor bandwidth 80 Hz.  The level is the hull broadband RMS
before the sensor's anti-alias filter, so a quiet, well-isolated
installation sits at the low end and a sensor near the engine bed at the
high end.

| Family | Hull level [m/s²] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2733 | 0.4726 | 0.1657 | 0.158 | 0.034 | 0.8199 | 1.00 |
| OU-II | 0.15 | 0.1136 | 0.2720 | 0.6486 | 1.1267 | 1.039 | 0.332 | 16.0225 | 1.37 |
| OU-II | 0.30 | 0.2272 | 0.4473 | 2.2930 | 1.1369 | 1.098 | 1.154 | 97.7385 | 4.85 |
| OU-II | 1.20 | 0.9087 | 24.3077 | 41.7513 | 3.9588 | 2.988 | 19.823 | 104.0378 | 88.35 |
| OU-II | 2.40 | 1.8173 | 130.6178 | 139.0780 | 16.1623 | 13.078 | 126.697 | 104.2660 | 294.31 |
| OU-III | off | 0.0000 | 0.1687 | 0.3751 | 0.1713 | 0.161 | 0.020 | 0.8858 | 1.00 |
| OU-III | 0.15 | 0.1136 | 0.3199 | 0.7367 | 1.5634 | 1.539 | 0.327 | 23.8601 | 1.96 |
| OU-III | 0.30 | 0.2272 | 0.7770 | 3.3272 | 0.8371 | 0.718 | 0.890 | 90.3571 | 8.87 |
| OU-III | 1.20 | 0.9087 | 66.6914 | 89.6585 | 23.1457 | 19.259 | 63.107 | 103.3309 | 239.04 |
| OU-III | 2.40 | 1.8173 | 140.5605 | 150.2106 | 26.3342 | 22.016 | 133.235 | 102.8081 | 400.47 |
| TFG | off | 0.0000 | 0.1241 | 0.3597 | 0.1645 | 0.143 | 0.030 | 0.7567 | 1.00 |
| TFG | 0.15 | 0.1136 | 0.1442 | 0.8984 | 4.2288 | 3.970 | 0.420 | 12.7985 | 2.50 |
| TFG | 0.30 | 0.2272 | 0.3555 | 3.0610 | 11.2067 | 10.871 | 1.923 | 47.1563 | 8.51 |
| TFG | 1.20 | 0.9087 | 118.6651 | 204.0832 | 39.6248 | 34.232 | 72.582 | 103.7561 | 567.35 |
| TFG | 2.40 | 1.8173 | 125.6716 | 269.0608 | 27.4658 | 20.305 | 95.970 | 103.7708 | 747.98 |

## Sensor anti-alias bandwidth

Engine speed fixed at 2400 rpm and the hull level at
0.60 m/s².  A narrower filter removes power *and*
stops the high crank orders from folding below Nyquist.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2733 | 0.4726 | 0.1657 | 0.158 | 0.034 | 0.8199 | 1.00 |
| OU-II | 20 | 0.1865 | 0.2736 | 0.5066 | 0.3211 | 0.305 | 0.103 | 15.1472 | 1.07 |
| OU-II | 40 | 0.2939 | 0.4768 | 2.4911 | 0.7555 | 0.654 | 1.110 | 107.2246 | 5.27 |
| OU-II | 160 | 0.6034 | 0.9996 | 5.9641 | 1.0652 | 0.587 | 1.047 | 104.0505 | 12.62 |
| OU-III | off | 0.0000 | 0.1687 | 0.3751 | 0.1713 | 0.161 | 0.020 | 0.8858 | 1.00 |
| OU-III | 20 | 0.1865 | 0.1996 | 0.5416 | 1.2560 | 1.202 | 0.119 | 16.8196 | 1.44 |
| OU-III | 40 | 0.2939 | 0.8102 | 3.4042 | 0.5799 | 0.420 | 0.892 | 90.5178 | 9.08 |
| OU-III | 160 | 0.6034 | 1.4200 | 8.9598 | 1.1965 | 0.529 | 1.223 | 98.6102 | 23.89 |
| TFG | off | 0.0000 | 0.1241 | 0.3597 | 0.1645 | 0.143 | 0.030 | 0.7567 | 1.00 |
| TFG | 20 | 0.1865 | 0.1306 | 0.5333 | 3.0996 | 2.957 | 0.189 | 12.2284 | 1.48 |
| TFG | 40 | 0.2939 | 0.2816 | 2.6205 | 10.1288 | 9.789 | 1.640 | 44.6619 | 7.29 |
| TFG | 160 | 0.6034 | 5.5751 | 49.2063 | 10.2062 | 7.015 | 6.867 | 106.4940 | 136.79 |

## Matched-power control

The same bandwidth sweep with the hull level rescaled so every cell
records the same vibration RMS as the nominal cruise cell.  The folded
line frequencies still differ from cell to cell; the recorded power no
longer does.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2733 | 0.4726 | 0.1657 | 0.158 | 0.034 | 0.8199 | 1.00 |
| OU-II | 20 | 0.4543 | 0.3635 | 1.8694 | 1.2894 | 1.222 | 0.851 | 96.3430 | 3.96 |
| OU-II | 40 | 0.4543 | 0.8111 | 4.1072 | 1.1141 | 0.963 | 1.051 | 104.0683 | 8.69 |
| OU-II | 160 | 0.4543 | 0.6567 | 3.4743 | 0.8286 | 0.536 | 0.873 | 103.9331 | 7.35 |
| OU-III | off | 0.0000 | 0.1687 | 0.3751 | 0.1713 | 0.161 | 0.020 | 0.8858 | 1.00 |
| OU-III | 20 | 0.4543 | 0.6018 | 2.4160 | 0.6301 | 0.516 | 0.793 | 104.6442 | 6.44 |
| OU-III | 40 | 0.4543 | 1.3241 | 6.3331 | 0.9148 | 0.518 | 1.338 | 104.1810 | 16.88 |
| OU-III | 160 | 0.4543 | 1.1085 | 6.3158 | 0.9806 | 0.494 | 1.027 | 90.5525 | 16.84 |
| TFG | off | 0.0000 | 0.1241 | 0.3597 | 0.1645 | 0.143 | 0.030 | 0.7567 | 1.00 |
| TFG | 20 | 0.4543 | 0.2045 | 2.0646 | 8.4185 | 8.150 | 1.328 | 60.4624 | 5.74 |
| TFG | 40 | 0.4543 | 1.1876 | 14.8698 | 8.6629 | 7.426 | 6.106 | 110.4624 | 41.34 |
| TFG | 160 | 0.4543 | 1.7413 | 19.4994 | 6.6648 | 6.043 | 6.792 | 96.4814 | 54.21 |

## Sensor attribution

The nominal cruise cell rerun with the model's gyroscope terms switched
off, so the accelerometer is the only perturbed sensor.  Compare against
the 2400 rpm row of the engine-speed table.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2733 | 0.4726 | 0.1657 | 0.158 | 0.034 | 0.8199 | 1.00 |
| OU-II | 2400 | 0.4543 | 1.4430 | 6.9966 | 1.4783 | 1.223 | 1.514 | 104.0904 | 14.81 |
| OU-III | off | 0.0000 | 0.1687 | 0.3751 | 0.1713 | 0.161 | 0.020 | 0.8858 | 1.00 |
| OU-III | 2400 | 0.4543 | 1.4940 | 8.4589 | 1.1974 | 0.602 | 1.243 | 98.8594 | 22.55 |
| TFG | off | 0.0000 | 0.1241 | 0.3597 | 0.1645 | 0.143 | 0.030 | 0.7567 | 1.00 |
| TFG | 2400 | 0.4543 | 6.2587 | 51.3147 | 11.2561 | 8.442 | 7.081 | 106.2153 | 142.65 |

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
