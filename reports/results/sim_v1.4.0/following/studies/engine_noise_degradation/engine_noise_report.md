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

Source commit used for the replay: `71f410938ca1a8f393979f4933933868fb612461`.

## Engine speed

Vibration level fixed at 0.60 m/s² of hull broadband
RMS at 2400 rpm, sensor bandwidth 80 Hz.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3238 | 0.5762 | 0.1468 | 0.137 | 0.044 | 0.8278 | 1.00 |
| OU-II | 800 | 0.3751 | 1.5078 | 8.1210 | 1.9860 | 1.779 | 1.566 | 104.0233 | 14.09 |
| OU-II | 1200 | 0.3943 | 2.5681 | 12.7233 | 3.0123 | 2.776 | 2.467 | 104.0243 | 22.08 |
| OU-II | 1600 | 0.3810 | 4.9451 | 20.6071 | 3.7306 | 3.259 | 4.519 | 103.9159 | 35.77 |
| OU-II | 2000 | 0.4106 | 8.9140 | 36.5515 | 4.1985 | 3.021 | 8.289 | 103.8537 | 63.44 |
| OU-II | 2400 | 0.4543 | 12.3733 | 47.7629 | 3.7321 | 1.389 | 10.204 | 104.1173 | 82.90 |
| OU-II | 2800 | 0.5067 | 13.0675 | 49.0962 | 4.7225 | 2.959 | 11.217 | 103.8144 | 85.21 |
| OU-II | 3200 | 0.5657 | 18.1149 | 60.4244 | 8.3351 | 7.132 | 17.216 | 104.0392 | 104.87 |
| OU-III | off | 0.0000 | 0.2092 | 0.4853 | 0.1046 | 0.091 | 0.042 | 0.6600 | 1.00 |
| OU-III | 800 | 0.3751 | 2.5204 | 13.5467 | 2.0244 | 1.521 | 2.157 | 90.2562 | 27.91 |
| OU-III | 1200 | 0.3943 | 2.9318 | 18.9697 | 2.7560 | 2.084 | 2.060 | 90.8449 | 39.09 |
| OU-III | 1600 | 0.3810 | 4.1310 | 26.4915 | 2.9118 | 1.976 | 2.772 | 90.0991 | 54.59 |
| OU-III | 2000 | 0.4106 | 4.3083 | 32.0563 | 2.7580 | 1.423 | 2.537 | 90.0485 | 66.05 |
| OU-III | 2400 | 0.4543 | 73.8214 | 100.2985 | 13.2227 | 5.996 | 61.773 | 89.8144 | 206.67 |
| OU-III | 2800 | 0.5067 | 43.6604 | 78.7298 | 10.8901 | 8.631 | 35.533 | 90.1696 | 162.22 |
| OU-III | 3200 | 0.5657 | 111.2986 | 145.3359 | 16.8877 | 10.496 | 94.414 | 90.5346 | 299.47 |
| TFG | off | 0.0000 | 0.1470 | 0.4696 | 0.1346 | 0.124 | 0.043 | 0.6638 | 1.00 |
| TFG | 800 | 0.3751 | 0.8016 | 16.9447 | 7.6904 | 7.189 | 13.954 | 103.4041 | 36.09 |
| TFG | 1200 | 0.3943 | 2.6749 | 26.9495 | 6.6186 | 5.854 | 11.169 | 103.7063 | 57.39 |
| TFG | 1600 | 0.3810 | 2.1922 | 21.8564 | 9.7558 | 8.308 | 11.358 | 96.0921 | 46.55 |
| TFG | 2000 | 0.4106 | 4.6666 | 41.1194 | 9.1203 | 6.636 | 9.352 | 95.7954 | 87.57 |
| TFG | 2400 | 0.4543 | 17.4442 | 79.2994 | 9.2119 | 6.070 | 20.179 | 95.0374 | 168.88 |
| TFG | 2800 | 0.5067 | 23.3401 | 108.5013 | 12.4982 | 7.288 | 28.349 | 104.2912 | 231.07 |
| TFG | 3200 | 0.5657 | 27.2735 | 121.4831 | 13.4510 | 6.886 | 29.177 | 105.0608 | 258.72 |

## Vibration level

Engine speed fixed at 2400 rpm, sensor bandwidth 80 Hz.  The level is the hull broadband RMS
before the sensor's anti-alias filter, so a quiet, well-isolated
installation sits at the low end and a sensor near the engine bed at the
high end.

| Family | Hull level [m/s²] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3238 | 0.5762 | 0.1468 | 0.137 | 0.044 | 0.8278 | 1.00 |
| OU-II | 0.15 | 0.1136 | 0.3581 | 1.7353 | 1.0814 | 0.947 | 0.930 | 90.5313 | 3.01 |
| OU-II | 0.30 | 0.2272 | 1.2907 | 6.4693 | 1.5226 | 1.043 | 0.948 | 90.1518 | 11.23 |
| OU-II | 1.20 | 0.9087 | 287.8528 | 297.4451 | 22.7467 | 8.165 | 269.899 | 103.1756 | 516.25 |
| OU-II | 2.40 | 1.8173 | 551.9248 | 562.3188 | 29.4063 | 18.061 | 520.968 | 103.9067 | 975.98 |
| OU-III | off | 0.0000 | 0.2092 | 0.4853 | 0.1046 | 0.091 | 0.042 | 0.6600 | 1.00 |
| OU-III | 0.15 | 0.1136 | 0.4009 | 1.4591 | 1.6565 | 1.618 | 0.667 | 44.8240 | 3.01 |
| OU-III | 0.30 | 0.2272 | 0.9136 | 5.4088 | 1.0881 | 0.476 | 0.658 | 90.0911 | 11.15 |
| OU-III | 1.20 | 0.9087 | 176.1639 | 197.8592 | 20.6901 | 6.252 | 158.128 | 91.9801 | 407.69 |
| OU-III | 2.40 | 1.8173 | 331.1704 | 346.3468 | 30.3227 | 20.121 | 303.055 | 103.4460 | 713.65 |
| TFG | off | 0.0000 | 0.1470 | 0.4696 | 0.1346 | 0.124 | 0.043 | 0.6638 | 1.00 |
| TFG | 0.15 | 0.1136 | 0.1669 | 0.8842 | 3.0947 | 2.936 | 0.515 | 8.9289 | 1.88 |
| TFG | 0.30 | 0.2272 | 0.4079 | 2.7919 | 9.6090 | 9.384 | 1.959 | 36.8088 | 5.95 |
| TFG | 1.20 | 0.9087 | 24.7860 | 106.6748 | 25.2051 | 21.321 | 24.386 | 104.0998 | 227.18 |
| TFG | 2.40 | 1.8173 | 4.9268 | 8.7697 | 16.9125 | 15.773 | 5.863 | 103.8478 | 18.68 |

## Sensor anti-alias bandwidth

Engine speed fixed at 2400 rpm and the hull level at
0.60 m/s².  A narrower filter removes power *and*
stops the high crank orders from folding below Nyquist.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3238 | 0.5762 | 0.1468 | 0.137 | 0.044 | 0.8278 | 1.00 |
| OU-II | 20 | 0.1865 | 0.3405 | 1.9716 | 1.2075 | 1.092 | 1.160 | 94.4983 | 3.42 |
| OU-II | 40 | 0.2939 | 1.4867 | 7.0839 | 1.9675 | 1.555 | 1.287 | 91.3581 | 12.29 |
| OU-II | 160 | 0.6034 | 14.6282 | 54.4300 | 3.8879 | 0.917 | 10.443 | 100.0233 | 94.47 |
| OU-III | off | 0.0000 | 0.2092 | 0.4853 | 0.1046 | 0.091 | 0.042 | 0.6600 | 1.00 |
| OU-III | 20 | 0.1865 | 0.4214 | 2.3088 | 0.4995 | 0.211 | 0.843 | 89.5700 | 4.76 |
| OU-III | 40 | 0.2939 | 1.3715 | 8.5672 | 1.7591 | 1.271 | 0.861 | 90.1325 | 17.65 |
| OU-III | 160 | 0.6034 | 4.9569 | 34.5222 | 2.2890 | 0.227 | 3.293 | 89.9495 | 71.13 |
| TFG | off | 0.0000 | 0.1470 | 0.4696 | 0.1346 | 0.124 | 0.043 | 0.6638 | 1.00 |
| TFG | 20 | 0.1865 | 0.3451 | 2.6125 | 10.1333 | 9.814 | 1.807 | 36.8113 | 5.56 |
| TFG | 40 | 0.2939 | 3.6262 | 33.8834 | 9.1213 | 8.485 | 10.335 | 88.6852 | 72.16 |
| TFG | 160 | 0.6034 | 5.1800 | 47.3661 | 8.9267 | 7.792 | 9.968 | 86.8850 | 100.87 |

## Matched-power control

The same bandwidth sweep with the hull level rescaled so every cell
records the same vibration RMS as the nominal cruise cell.  The folded
line frequencies still differ from cell to cell; the recorded power no
longer does.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3238 | 0.5762 | 0.1468 | 0.137 | 0.044 | 0.8278 | 1.00 |
| OU-II | 20 | 0.4543 | 1.1000 | 5.8454 | 3.2164 | 3.012 | 1.146 | 104.1335 | 10.15 |
| OU-II | 40 | 0.4543 | 5.0233 | 24.4164 | 4.3007 | 3.438 | 5.424 | 104.0565 | 42.38 |
| OU-II | 160 | 0.4543 | 4.1815 | 22.0041 | 2.2583 | 0.502 | 3.846 | 90.7120 | 38.19 |
| OU-III | off | 0.0000 | 0.2092 | 0.4853 | 0.1046 | 0.091 | 0.042 | 0.6600 | 1.00 |
| OU-III | 20 | 0.4543 | 2.0733 | 14.3895 | 2.2143 | 1.341 | 0.995 | 90.7049 | 29.65 |
| OU-III | 40 | 0.4543 | 4.1356 | 29.7731 | 2.9796 | 1.883 | 2.149 | 90.3028 | 61.35 |
| OU-III | 160 | 0.4543 | 2.2959 | 13.5237 | 1.3878 | 0.213 | 1.289 | 89.8234 | 27.87 |
| TFG | off | 0.0000 | 0.1470 | 0.4696 | 0.1346 | 0.124 | 0.043 | 0.6638 | 1.00 |
| TFG | 20 | 0.4543 | 16.6880 | 90.0635 | 12.1239 | 7.920 | 10.743 | 103.8866 | 191.80 |
| TFG | 40 | 0.4543 | 12.6557 | 79.6719 | 10.3846 | 5.739 | 16.049 | 105.5473 | 169.67 |
| TFG | 160 | 0.4543 | 2.2847 | 23.4500 | 9.1427 | 8.813 | 9.629 | 70.5079 | 49.94 |

## Sensor attribution

The nominal cruise cell rerun with the model's gyroscope terms switched
off, so the accelerometer is the only perturbed sensor.  Compare against
the 2400 rpm row of the engine-speed table.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3238 | 0.5762 | 0.1468 | 0.137 | 0.044 | 0.8278 | 1.00 |
| OU-II | 2400 | 0.4543 | 12.6615 | 48.2957 | 3.7578 | 1.382 | 10.468 | 104.1990 | 83.82 |
| OU-III | off | 0.0000 | 0.2092 | 0.4853 | 0.1046 | 0.091 | 0.042 | 0.6600 | 1.00 |
| OU-III | 2400 | 0.4543 | 52.2161 | 84.8159 | 10.4849 | 7.508 | 40.556 | 89.6552 | 174.76 |
| TFG | off | 0.0000 | 0.1470 | 0.4696 | 0.1346 | 0.124 | 0.043 | 0.6638 | 1.00 |
| TFG | 2400 | 0.4543 | 17.6919 | 79.4720 | 9.5540 | 6.320 | 13.188 | 95.1478 | 169.25 |

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
