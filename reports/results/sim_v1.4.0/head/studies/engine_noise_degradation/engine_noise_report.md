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
| OU-II | off | 0.0000 | 0.2733 | 0.4725 | 0.1657 | 0.158 | 0.034 | 0.8197 | 1.00 |
| OU-II | 800 | 0.3751 | 0.3876 | 1.7572 | 1.2547 | 1.208 | 1.183 | 86.0945 | 3.72 |
| OU-II | 1200 | 0.3943 | 0.5085 | 2.4836 | 1.6795 | 1.638 | 1.527 | 91.8348 | 5.26 |
| OU-II | 1600 | 0.3810 | 0.6019 | 3.0278 | 1.5682 | 1.523 | 1.022 | 105.8549 | 6.41 |
| OU-II | 2000 | 0.4106 | 0.9417 | 4.7608 | 1.7373 | 1.633 | 1.027 | 104.3929 | 10.07 |
| OU-II | 2400 | 0.4543 | 1.4296 | 6.9332 | 1.4799 | 1.229 | 1.472 | 104.0526 | 14.67 |
| OU-II | 2800 | 0.5067 | 1.8017 | 8.0182 | 1.3486 | 0.980 | 1.838 | 104.2763 | 16.97 |
| OU-II | 3200 | 0.5657 | 1.8231 | 8.1252 | 1.3292 | 1.133 | 1.855 | 104.0913 | 17.19 |
| OU-III | off | 0.0000 | 0.1687 | 0.3751 | 0.1712 | 0.161 | 0.020 | 0.8846 | 1.00 |
| OU-III | 800 | 0.3751 | 0.7256 | 2.6123 | 0.7188 | 0.547 | 1.003 | 107.2921 | 6.96 |
| OU-III | 1200 | 0.3943 | 0.8319 | 3.2581 | 1.3140 | 1.246 | 1.003 | 105.3365 | 8.69 |
| OU-III | 1600 | 0.3810 | 1.1780 | 5.6720 | 1.9433 | 1.881 | 1.223 | 99.1139 | 15.12 |
| OU-III | 2000 | 0.4106 | 1.4035 | 7.7960 | 1.5595 | 1.270 | 1.224 | 104.1772 | 20.78 |
| OU-III | 2400 | 0.4543 | 1.4650 | 8.5721 | 1.2205 | 0.639 | 1.247 | 98.8217 | 22.85 |
| OU-III | 2800 | 0.5067 | 1.7274 | 9.1243 | 1.0848 | 0.432 | 1.487 | 104.1722 | 24.33 |
| OU-III | 3200 | 0.5657 | 2.1371 | 10.0299 | 1.4301 | 1.143 | 2.111 | 103.5571 | 26.74 |
| TFG | off | 0.0000 | 0.1241 | 0.2741 | 0.1567 | 0.138 | 0.024 | 0.7412 | 1.00 |
| TFG | 800 | 0.3751 | 0.3458 | 2.5318 | 10.9208 | 10.601 | 1.431 | 68.2243 | 9.24 |
| TFG | 1200 | 0.3943 | 0.7396 | 7.8206 | 12.4248 | 12.082 | 3.949 | 97.4515 | 28.54 |
| TFG | 1600 | 0.3810 | 1.9558 | 18.1266 | 5.9612 | 4.139 | 5.421 | 110.0544 | 66.14 |
| TFG | 2000 | 0.4106 | 5.4580 | 41.2652 | 10.1541 | 7.582 | 4.218 | 105.1403 | 150.57 |
| TFG | 2400 | 0.4543 | 5.3140 | 41.9327 | 11.6475 | 9.008 | 5.181 | 106.1694 | 153.01 |
| TFG | 2800 | 0.5067 | 4.0069 | 27.5393 | 7.0169 | 6.034 | 5.573 | 103.9674 | 100.49 |
| TFG | 3200 | 0.5657 | 0.8920 | 10.9334 | 12.0506 | 11.268 | 6.030 | 108.3031 | 39.90 |

## Vibration level

Engine speed fixed at 2400 rpm, sensor bandwidth 80 Hz.  The level is the hull broadband RMS
before the sensor's anti-alias filter, so a quiet, well-isolated
installation sits at the low end and a sensor near the engine bed at the
high end.

| Family | Hull level [m/s²] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2733 | 0.4725 | 0.1657 | 0.158 | 0.034 | 0.8197 | 1.00 |
| OU-II | 0.15 | 0.1136 | 0.2720 | 0.6479 | 1.1255 | 1.038 | 0.332 | 15.9922 | 1.37 |
| OU-II | 0.30 | 0.2272 | 0.4428 | 2.2418 | 1.1617 | 1.122 | 1.139 | 96.6929 | 4.74 |
| OU-II | 1.20 | 0.9087 | 51.9829 | 59.3785 | 8.9886 | 4.866 | 50.711 | 103.7124 | 125.66 |
| OU-II | 2.40 | 1.8173 | 141.2969 | 149.0672 | 16.7036 | 12.887 | 135.409 | 103.9246 | 315.46 |
| OU-III | off | 0.0000 | 0.1687 | 0.3751 | 0.1712 | 0.161 | 0.020 | 0.8846 | 1.00 |
| OU-III | 0.15 | 0.1136 | 0.3196 | 0.7346 | 1.5618 | 1.537 | 0.326 | 23.7581 | 1.96 |
| OU-III | 0.30 | 0.2272 | 0.7742 | 3.3225 | 0.8339 | 0.714 | 0.871 | 90.4449 | 8.86 |
| OU-III | 1.20 | 0.9087 | 85.3566 | 97.1728 | 27.4011 | 22.879 | 79.436 | 102.8810 | 259.07 |
| OU-III | 2.40 | 1.8173 | 147.1453 | 155.0035 | 33.8209 | 25.667 | 138.798 | 102.8208 | 413.25 |
| TFG | off | 0.0000 | 0.1241 | 0.2741 | 0.1567 | 0.138 | 0.024 | 0.7412 | 1.00 |
| TFG | 0.15 | 0.1136 | 0.1434 | 0.7655 | 4.0892 | 3.840 | 0.308 | 12.4452 | 2.79 |
| TFG | 0.30 | 0.2272 | 0.3442 | 2.5673 | 10.9655 | 10.626 | 1.417 | 45.8036 | 9.37 |
| TFG | 1.20 | 0.9087 | 68.1324 | 157.2660 | 36.7324 | 31.002 | 39.923 | 103.8425 | 573.85 |
| TFG | 2.40 | 1.8173 | 92.8881 | 147.1834 | 21.7886 | 15.351 | 47.207 | 104.0983 | 537.06 |

## Sensor anti-alias bandwidth

Engine speed fixed at 2400 rpm and the hull level at
0.60 m/s².  A narrower filter removes power *and*
stops the high crank orders from folding below Nyquist.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2733 | 0.4725 | 0.1657 | 0.158 | 0.034 | 0.8197 | 1.00 |
| OU-II | 20 | 0.1865 | 0.2736 | 0.5058 | 0.3187 | 0.303 | 0.103 | 15.1141 | 1.07 |
| OU-II | 40 | 0.2939 | 0.4980 | 2.5325 | 0.7738 | 0.664 | 1.230 | 108.1291 | 5.36 |
| OU-II | 160 | 0.6034 | 0.9992 | 5.9514 | 1.0700 | 0.595 | 1.064 | 104.0267 | 12.59 |
| OU-III | off | 0.0000 | 0.1687 | 0.3751 | 0.1712 | 0.161 | 0.020 | 0.8846 | 1.00 |
| OU-III | 20 | 0.1865 | 0.1994 | 0.5389 | 1.2501 | 1.194 | 0.118 | 16.7069 | 1.44 |
| OU-III | 40 | 0.2939 | 0.8074 | 3.3818 | 0.5762 | 0.417 | 0.915 | 90.3909 | 9.02 |
| OU-III | 160 | 0.6034 | 1.4150 | 8.9222 | 1.1972 | 0.529 | 1.223 | 98.4128 | 23.79 |
| TFG | off | 0.0000 | 0.1241 | 0.2741 | 0.1567 | 0.138 | 0.024 | 0.7412 | 1.00 |
| TFG | 20 | 0.1865 | 0.1304 | 0.4387 | 3.0447 | 2.906 | 0.138 | 12.1116 | 1.60 |
| TFG | 40 | 0.2939 | 0.2748 | 2.1969 | 9.9052 | 9.566 | 1.208 | 43.7425 | 8.02 |
| TFG | 160 | 0.6034 | 5.3404 | 41.3443 | 10.7137 | 7.780 | 5.898 | 106.5130 | 150.86 |

## Matched-power control

The same bandwidth sweep with the hull level rescaled so every cell
records the same vibration RMS as the nominal cruise cell.  The folded
line frequencies still differ from cell to cell; the recorded power no
longer does.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2733 | 0.4725 | 0.1657 | 0.158 | 0.034 | 0.8197 | 1.00 |
| OU-II | 20 | 0.4543 | 0.3685 | 1.8513 | 1.2947 | 1.225 | 0.872 | 95.6360 | 3.92 |
| OU-II | 40 | 0.4543 | 0.8094 | 4.0928 | 1.1165 | 0.967 | 1.032 | 104.0519 | 8.66 |
| OU-II | 160 | 0.4543 | 0.6553 | 3.4526 | 0.8250 | 0.535 | 0.839 | 103.8994 | 7.31 |
| OU-III | off | 0.0000 | 0.1687 | 0.3751 | 0.1712 | 0.161 | 0.020 | 0.8846 | 1.00 |
| OU-III | 20 | 0.4543 | 0.5960 | 2.3789 | 0.6331 | 0.517 | 0.785 | 104.7988 | 6.34 |
| OU-III | 40 | 0.4543 | 1.3258 | 6.4289 | 0.9094 | 0.513 | 1.330 | 104.3738 | 17.14 |
| OU-III | 160 | 0.4543 | 1.1085 | 6.2703 | 0.9804 | 0.498 | 1.022 | 90.5834 | 16.72 |
| TFG | off | 0.0000 | 0.1241 | 0.2741 | 0.1567 | 0.138 | 0.024 | 0.7412 | 1.00 |
| TFG | 20 | 0.4543 | 0.2017 | 1.7120 | 8.2200 | 7.953 | 0.977 | 60.0888 | 6.25 |
| TFG | 40 | 0.4543 | 1.1471 | 12.5479 | 8.7252 | 7.521 | 5.185 | 109.6338 | 45.79 |
| TFG | 160 | 0.4543 | 1.6822 | 16.5747 | 6.5949 | 5.994 | 5.856 | 96.4466 | 60.48 |

## Sensor attribution

The nominal cruise cell rerun with the model's gyroscope terms switched
off, so the accelerometer is the only perturbed sensor.  Compare against
the 2400 rpm row of the engine-speed table.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2733 | 0.4725 | 0.1657 | 0.158 | 0.034 | 0.8197 | 1.00 |
| OU-II | 2400 | 0.4543 | 1.4417 | 6.9771 | 1.4770 | 1.223 | 1.484 | 104.0374 | 14.77 |
| OU-III | off | 0.0000 | 0.1687 | 0.3751 | 0.1712 | 0.161 | 0.020 | 0.8846 | 1.00 |
| OU-III | 2400 | 0.4543 | 1.4840 | 8.6399 | 1.2226 | 0.637 | 1.273 | 98.9094 | 23.03 |
| TFG | off | 0.0000 | 0.1241 | 0.2741 | 0.1567 | 0.138 | 0.024 | 0.7412 | 1.00 |
| TFG | 2400 | 0.4543 | 5.4524 | 42.4452 | 11.6785 | 9.065 | 5.195 | 106.1425 | 154.88 |

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
