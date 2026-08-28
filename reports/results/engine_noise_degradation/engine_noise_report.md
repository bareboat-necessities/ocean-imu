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

Source commit used for the replay: `1e9b4e15098af637466a3b54373cb0787bca5b37`.

## Engine speed

Vibration level fixed at 0.60 m/s² of hull broadband
RMS at 2400 rpm, sensor bandwidth 80 Hz.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2794 | 0.6423 | 0.1965 | 0.162 | 0.178 | 0.6669 | 1.00 |
| OU-II | 800 | 0.3751 | 0.2969 | 0.8021 | 0.2509 | 0.203 | 0.272 | 5.1500 | 1.25 |
| OU-II | 1200 | 0.3943 | 0.7223 | 2.2745 | 2.3048 | 2.289 | 2.019 | 5.2214 | 3.54 |
| OU-II | 1600 | 0.3810 | 0.7007 | 2.3793 | 1.9178 | 1.887 | 2.012 | 10.4929 | 3.70 |
| OU-II | 2000 | 0.4106 | 0.4749 | 1.7384 | 1.4915 | 1.292 | 1.004 | 48.3277 | 2.71 |
| OU-II | 2400 | 0.4543 | 0.7400 | 3.3627 | 2.4362 | 2.303 | 2.862 | 58.5893 | 5.24 |
| OU-II | 2800 | 0.5067 | 1.2776 | 6.2175 | 2.7304 | 2.655 | 5.723 | 60.5089 | 9.68 |
| OU-II | 3200 | 0.5657 | 1.8877 | 10.4846 | 3.1089 | 3.075 | 9.992 | 58.9138 | 16.32 |
| OU-III | off | 0.0000 | 0.1791 | 0.5224 | 0.1582 | 0.144 | 0.140 | 0.6292 | 1.00 |
| OU-III | 800 | 0.3751 | 0.3275 | 0.8770 | 2.0752 | 2.014 | 0.511 | 12.6122 | 1.68 |
| OU-III | 1200 | 0.3943 | 0.3874 | 1.2680 | 1.7011 | 1.601 | 0.915 | 10.8525 | 2.43 |
| OU-III | 1600 | 0.3810 | 0.3825 | 1.4021 | 1.3930 | 1.211 | 0.905 | 17.8497 | 2.68 |
| OU-III | 2000 | 0.4106 | 0.5349 | 2.2379 | 2.7566 | 2.746 | 1.646 | 36.5838 | 4.28 |
| OU-III | 2400 | 0.4543 | 0.8279 | 4.2318 | 3.0844 | 3.075 | 3.875 | 45.3850 | 8.10 |
| OU-III | 2800 | 0.5067 | 1.1654 | 6.7458 | 3.2900 | 3.280 | 6.268 | 54.0595 | 12.91 |
| OU-III | 3200 | 0.5657 | 1.3978 | 8.6509 | 3.4297 | 3.283 | 8.282 | 66.5513 | 16.56 |
| TFG | off | 0.0000 | 0.1897 | 0.7446 | 0.3860 | 0.353 | 0.206 | 0.8042 | 1.00 |
| TFG | 800 | 0.3751 | 0.3738 | 1.4042 | 3.5404 | 3.471 | 0.980 | 9.6493 | 1.89 |
| TFG | 1200 | 0.3943 | 0.6071 | 2.2685 | 5.4900 | 5.343 | 1.814 | 13.1056 | 3.05 |
| TFG | 1600 | 0.3810 | 0.3974 | 1.7602 | 4.3086 | 4.190 | 1.080 | 11.3392 | 2.36 |
| TFG | 2000 | 0.4106 | 1.0273 | 5.1999 | 12.0642 | 11.947 | 4.470 | 38.9786 | 6.98 |
| TFG | 2400 | 0.4543 | 1.6100 | 14.0263 | 16.2053 | 16.127 | 12.903 | 74.3483 | 18.84 |
| TFG | 2800 | 0.5067 | 13.1957 | 107.1635 | 15.3802 | 9.643 | 104.235 | 108.7871 | 143.93 |
| TFG | 3200 | 0.5657 | 13.7924 | 108.2354 | 14.8353 | 5.215 | 105.011 | 121.9588 | 145.36 |

## Vibration level

Engine speed fixed at 2400 rpm, sensor bandwidth 80 Hz.  The level is the hull broadband RMS
before the sensor's anti-alias filter, so a quiet, well-isolated
installation sits at the low end and a sensor near the engine bed at the
high end.

| Family | Hull level [m/s²] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2794 | 0.6423 | 0.1965 | 0.162 | 0.178 | 0.6669 | 1.00 |
| OU-II | 0.15 | 0.1136 | 0.2732 | 0.6907 | 1.4496 | 1.336 | 0.289 | 5.7241 | 1.08 |
| OU-II | 0.30 | 0.2272 | 0.3419 | 1.0068 | 2.3116 | 2.248 | 0.638 | 13.7227 | 1.57 |
| OU-II | 1.20 | 0.9087 | 6.3098 | 38.0675 | 9.9458 | 9.536 | 15.695 | 103.7839 | 59.27 |
| OU-II | 2.40 | 1.8173 | 131.8057 | 147.3785 | 12.1401 | 7.830 | 125.315 | 103.8700 | 229.46 |
| OU-III | off | 0.0000 | 0.1791 | 0.5224 | 0.1582 | 0.144 | 0.140 | 0.6292 | 1.00 |
| OU-III | 0.15 | 0.1136 | 0.2723 | 0.7068 | 2.1120 | 1.957 | 0.403 | 7.3875 | 1.35 |
| OU-III | 0.30 | 0.2272 | 0.4391 | 1.3724 | 2.7577 | 2.744 | 1.153 | 15.2079 | 2.63 |
| OU-III | 1.20 | 0.9087 | 6.1144 | 37.5310 | 8.1179 | 7.705 | 10.637 | 103.8825 | 71.84 |
| OU-III | 2.40 | 1.8173 | 86.1529 | 121.1547 | 12.0612 | 5.368 | 79.214 | 103.9744 | 231.91 |
| TFG | off | 0.0000 | 0.1897 | 0.7446 | 0.3860 | 0.353 | 0.206 | 0.8042 | 1.00 |
| TFG | 0.15 | 0.1136 | 0.3266 | 1.1546 | 3.7174 | 3.560 | 0.843 | 9.6170 | 1.55 |
| TFG | 0.30 | 0.2272 | 0.8142 | 3.2714 | 10.0928 | 9.883 | 2.874 | 28.7637 | 4.39 |
| TFG | 1.20 | 0.9087 | 14.1465 | 110.2997 | 17.9828 | 10.002 | 106.557 | 118.2253 | 148.14 |
| TFG | 2.40 | 1.8173 | 28.7678 | 121.5837 | 17.3727 | 7.998 | 115.209 | 103.7911 | 163.29 |

## Sensor anti-alias bandwidth

Engine speed fixed at 2400 rpm and the hull level at
0.60 m/s².  A narrower filter removes power *and*
stops the high crank orders from folding below Nyquist.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2794 | 0.6423 | 0.1965 | 0.162 | 0.178 | 0.6669 | 1.00 |
| OU-II | 20 | 0.1865 | 0.3301 | 0.8685 | 2.2563 | 2.145 | 0.566 | 7.3300 | 1.35 |
| OU-II | 40 | 0.2939 | 0.5647 | 1.8946 | 2.7949 | 2.786 | 1.676 | 16.0341 | 2.95 |
| OU-II | 160 | 0.6034 | 1.5324 | 8.2069 | 3.2402 | 3.216 | 7.620 | 70.6418 | 12.78 |
| OU-III | off | 0.0000 | 0.1791 | 0.5224 | 0.1582 | 0.144 | 0.140 | 0.6292 | 1.00 |
| OU-III | 20 | 0.1865 | 0.3014 | 0.7913 | 2.3453 | 2.249 | 0.506 | 8.2050 | 1.51 |
| OU-III | 40 | 0.2939 | 0.5938 | 2.1993 | 2.9222 | 2.919 | 1.995 | 17.6960 | 4.21 |
| OU-III | 160 | 0.6034 | 1.1807 | 7.5695 | 2.7014 | 2.360 | 6.786 | 85.6991 | 14.49 |
| TFG | off | 0.0000 | 0.1897 | 0.7446 | 0.3860 | 0.353 | 0.206 | 0.8042 | 1.00 |
| TFG | 20 | 0.1865 | 0.4532 | 1.5741 | 6.4024 | 6.196 | 1.270 | 17.2386 | 2.11 |
| TFG | 40 | 0.2939 | 1.2438 | 6.2418 | 15.0069 | 14.868 | 5.643 | 52.8608 | 8.38 |
| TFG | 160 | 0.6034 | 22.7210 | 110.3906 | 16.1018 | 12.653 | 107.268 | 96.4302 | 148.26 |

## Matched-power control

The same bandwidth sweep with the hull level rescaled so every cell
records the same vibration RMS as the nominal cruise cell.  The folded
line frequencies still differ from cell to cell; the recorded power no
longer does.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2794 | 0.6423 | 0.1965 | 0.162 | 0.178 | 0.6669 | 1.00 |
| OU-II | 20 | 0.4543 | 0.6828 | 2.4724 | 2.9176 | 2.912 | 2.285 | 15.8407 | 3.85 |
| OU-II | 40 | 0.4543 | 0.9673 | 4.2143 | 2.9252 | 2.917 | 3.908 | 31.0985 | 6.56 |
| OU-II | 160 | 0.4543 | 1.0558 | 4.7903 | 2.9750 | 2.957 | 4.444 | 47.8148 | 7.46 |
| OU-III | off | 0.0000 | 0.1791 | 0.5224 | 0.1582 | 0.144 | 0.140 | 0.6292 | 1.00 |
| OU-III | 20 | 0.4543 | 0.5823 | 2.2978 | 2.9184 | 2.913 | 2.107 | 19.8624 | 4.40 |
| OU-III | 40 | 0.4543 | 0.9159 | 4.5457 | 3.1596 | 3.150 | 4.285 | 33.3650 | 8.70 |
| OU-III | 160 | 0.4543 | 0.7665 | 3.7766 | 3.0544 | 3.049 | 3.497 | 34.8317 | 7.23 |
| TFG | off | 0.0000 | 0.1897 | 0.7446 | 0.3860 | 0.353 | 0.206 | 0.8042 | 1.00 |
| TFG | 20 | 0.4543 | 1.3655 | 7.9144 | 15.6748 | 15.572 | 7.183 | 62.1873 | 10.63 |
| TFG | 40 | 0.4543 | 21.9552 | 108.9373 | 15.7816 | 11.840 | 106.232 | 102.1651 | 146.31 |
| TFG | 160 | 0.4543 | 1.5954 | 11.4030 | 16.3570 | 16.278 | 10.502 | 64.7276 | 15.31 |

## Sensor attribution

The nominal cruise cell rerun with the model's gyroscope terms switched
off, so the accelerometer is the only perturbed sensor.  Compare against
the 2400 rpm row of the engine-speed table.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2794 | 0.6423 | 0.1965 | 0.162 | 0.178 | 0.6669 | 1.00 |
| OU-II | 2400 | 0.4543 | 0.7402 | 3.3631 | 2.4357 | 2.303 | 2.862 | 58.5900 | 5.24 |
| OU-III | off | 0.0000 | 0.1791 | 0.5224 | 0.1582 | 0.144 | 0.140 | 0.6292 | 1.00 |
| OU-III | 2400 | 0.4543 | 0.8279 | 4.2316 | 3.0844 | 3.075 | 3.875 | 45.3841 | 8.10 |
| TFG | off | 0.0000 | 0.1897 | 0.7446 | 0.3860 | 0.353 | 0.206 | 0.8042 | 1.00 |
| TFG | 2400 | 0.4543 | 1.6102 | 14.0281 | 16.2045 | 16.126 | 12.905 | 74.3508 | 18.84 |

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
