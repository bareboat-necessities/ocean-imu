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

Source commit used for the replay: `f7cd0fa140a2f9647ef33d16bf40c25258442158`.

## Engine speed

Vibration level fixed at 0.60 m/s² of hull broadband
RMS at 2400 rpm, sensor bandwidth 80 Hz.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3249 | 0.5832 | 0.2350 | 0.218 | 0.046 | 0.7281 | 1.00 |
| OU-II | 800 | 0.3751 | 1.4131 | 8.4323 | 2.0574 | 1.924 | 1.481 | 103.9122 | 14.46 |
| OU-II | 1200 | 0.3943 | 2.4473 | 13.1846 | 3.3515 | 3.169 | 2.343 | 103.6841 | 22.61 |
| OU-II | 1600 | 0.3810 | 8.2836 | 21.4947 | 4.4083 | 4.076 | 6.904 | 91.9051 | 36.86 |
| OU-II | 2000 | 0.4106 | 30.5231 | 47.5301 | 5.2496 | 4.440 | 25.709 | 91.3362 | 81.50 |
| OU-II | 2400 | 0.4543 | 41.2519 | 62.7425 | 3.5968 | 1.379 | 35.799 | 92.0789 | 107.58 |
| OU-II | 2800 | 0.5067 | 40.4418 | 63.6832 | 6.2494 | 4.972 | 35.561 | 100.9689 | 109.19 |
| OU-II | 3200 | 0.5657 | 59.8499 | 90.3069 | 13.8027 | 12.752 | 55.837 | 104.2416 | 154.84 |
| OU-III | off | 0.0000 | 0.2106 | 0.4939 | 0.1436 | 0.133 | 0.045 | 0.6165 | 1.00 |
| OU-III | 800 | 0.3751 | 3.6329 | 18.1584 | 2.7038 | 2.498 | 3.254 | 90.1656 | 36.77 |
| OU-III | 1200 | 0.3943 | 7.9959 | 27.1992 | 4.3428 | 4.015 | 7.178 | 90.9724 | 55.07 |
| OU-III | 1600 | 0.3810 | 79.4084 | 118.9407 | 20.8805 | 17.402 | 73.811 | 89.9146 | 240.83 |
| OU-III | 2000 | 0.4106 | 14.8447 | 43.5922 | 4.3032 | 3.545 | 13.891 | 90.0211 | 88.26 |
| OU-III | 2400 | 0.4543 | 21.6181 | 54.6343 | 6.9697 | 5.744 | 17.437 | 89.9826 | 110.62 |
| OU-III | 2800 | 0.5067 | 74.0920 | 106.4156 | 17.3861 | 15.815 | 70.734 | 90.1934 | 215.47 |
| OU-III | 3200 | 0.5657 | 87.6732 | 114.2628 | 19.5095 | 17.651 | 83.441 | 90.5403 | 231.36 |
| TFG | off | 0.0000 | 0.2116 | 0.5007 | 0.2983 | 0.292 | 0.049 | 0.9720 | 1.00 |
| TFG | 800 | 0.3751 | 0.7318 | 14.5742 | 8.4240 | 8.109 | 12.841 | 103.1649 | 29.11 |
| TFG | 1200 | 0.3943 | 1.7436 | 18.3475 | 7.7004 | 7.406 | 11.233 | 103.2290 | 36.64 |
| TFG | 1600 | 0.3810 | 2.8961 | 25.5061 | 9.2044 | 8.575 | 9.789 | 104.1217 | 50.94 |
| TFG | 2000 | 0.4106 | 4.6472 | 39.2830 | 8.5633 | 7.659 | 9.753 | 101.2631 | 78.46 |
| TFG | 2400 | 0.4543 | 16.1362 | 80.9301 | 10.3295 | 6.012 | 10.797 | 105.6406 | 161.64 |
| TFG | 2800 | 0.5067 | 20.2902 | 90.6744 | 9.9044 | 4.726 | 14.935 | 105.3729 | 181.10 |
| TFG | 3200 | 0.5657 | 21.1029 | 99.3500 | 10.5627 | 5.252 | 15.065 | 103.9347 | 198.43 |

## Vibration level

Engine speed fixed at 2400 rpm, sensor bandwidth 80 Hz.  The level is the hull broadband RMS
before the sensor's anti-alias filter, so a quiet, well-isolated
installation sits at the low end and a sensor near the engine bed at the
high end.

| Family | Hull level [m/s²] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3249 | 0.5832 | 0.2350 | 0.218 | 0.046 | 0.7281 | 1.00 |
| OU-II | 0.15 | 0.1136 | 0.3591 | 1.7702 | 1.0947 | 0.956 | 0.864 | 89.7575 | 3.04 |
| OU-II | 0.30 | 0.2272 | 1.2844 | 6.1217 | 1.5876 | 1.234 | 0.909 | 90.2658 | 10.50 |
| OU-II | 1.20 | 0.9087 | 161.6778 | 174.7771 | 18.5684 | 16.242 | 158.233 | 104.0522 | 299.68 |
| OU-II | 2.40 | 1.8173 | 338.9756 | 349.9495 | 29.5111 | 26.238 | 328.708 | 103.1453 | 600.04 |
| OU-III | off | 0.0000 | 0.2106 | 0.4939 | 0.1436 | 0.133 | 0.045 | 0.6165 | 1.00 |
| OU-III | 0.15 | 0.1136 | 0.4559 | 1.8975 | 1.5551 | 1.426 | 0.847 | 61.6867 | 3.84 |
| OU-III | 0.30 | 0.2272 | 1.1440 | 5.8910 | 1.4733 | 1.197 | 0.885 | 89.8296 | 11.93 |
| OU-III | 1.20 | 0.9087 | 112.9868 | 147.7651 | 18.9279 | 16.366 | 107.520 | 91.4328 | 299.19 |
| OU-III | 2.40 | 1.8173 | 216.1141 | 226.7985 | 34.3755 | 31.288 | 208.678 | 103.6446 | 459.22 |
| TFG | off | 0.0000 | 0.2116 | 0.5007 | 0.2983 | 0.292 | 0.049 | 0.9720 | 1.00 |
| TFG | 0.15 | 0.1136 | 0.2364 | 0.9542 | 3.6630 | 3.471 | 0.544 | 11.2655 | 1.91 |
| TFG | 0.30 | 0.2272 | 0.6565 | 3.7021 | 11.1100 | 10.849 | 2.335 | 52.9883 | 7.39 |
| TFG | 1.20 | 0.9087 | 36.5177 | 138.7573 | 26.9945 | 20.020 | 21.089 | 103.6728 | 277.13 |
| TFG | 2.40 | 1.8173 | 3.0296 | 8.4919 | 15.1237 | 13.710 | 6.491 | 103.8590 | 16.96 |

## Sensor anti-alias bandwidth

Engine speed fixed at 2400 rpm and the hull level at
0.60 m/s².  A narrower filter removes power *and*
stops the high crank orders from folding below Nyquist.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3249 | 0.5832 | 0.2350 | 0.218 | 0.046 | 0.7281 | 1.00 |
| OU-II | 20 | 0.1865 | 0.3280 | 1.8219 | 1.2154 | 1.080 | 0.982 | 92.3833 | 3.12 |
| OU-II | 40 | 0.2939 | 1.4625 | 6.9600 | 1.9367 | 1.527 | 1.199 | 91.0462 | 11.93 |
| OU-II | 160 | 0.6034 | 31.8381 | 63.3265 | 3.9085 | 1.263 | 26.095 | 91.5702 | 108.58 |
| OU-III | off | 0.0000 | 0.2106 | 0.4939 | 0.1436 | 0.133 | 0.045 | 0.6165 | 1.00 |
| OU-III | 20 | 0.1865 | 0.3465 | 2.1157 | 0.5788 | 0.331 | 0.711 | 89.9469 | 4.28 |
| OU-III | 40 | 0.2939 | 1.9605 | 10.3407 | 1.9868 | 1.512 | 1.425 | 90.1402 | 20.94 |
| OU-III | 160 | 0.6034 | 12.0757 | 40.0926 | 2.1918 | 0.604 | 10.701 | 89.9079 | 81.18 |
| TFG | off | 0.0000 | 0.2116 | 0.5007 | 0.2983 | 0.292 | 0.049 | 0.9720 | 1.00 |
| TFG | 20 | 0.1865 | 0.4426 | 2.8270 | 11.9855 | 11.593 | 1.938 | 48.4872 | 5.65 |
| TFG | 40 | 0.2939 | 3.3419 | 30.2359 | 9.0731 | 7.786 | 9.734 | 102.8084 | 60.39 |
| TFG | 160 | 0.6034 | 5.2447 | 45.4693 | 9.3131 | 7.729 | 9.661 | 98.7113 | 90.81 |

## Matched-power control

The same bandwidth sweep with the hull level rescaled so every cell
records the same vibration RMS as the nominal cruise cell.  The folded
line frequencies still differ from cell to cell; the recorded power no
longer does.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3249 | 0.5832 | 0.2350 | 0.218 | 0.046 | 0.7281 | 1.00 |
| OU-II | 20 | 0.4543 | 1.0559 | 5.9358 | 3.5942 | 3.345 | 0.992 | 104.0415 | 10.18 |
| OU-II | 40 | 0.4543 | 12.8924 | 28.6520 | 5.7864 | 4.921 | 10.500 | 104.6201 | 49.13 |
| OU-II | 160 | 0.4543 | 6.4457 | 22.5162 | 2.1501 | 0.401 | 4.827 | 90.4515 | 38.61 |
| OU-III | off | 0.0000 | 0.2106 | 0.4939 | 0.1436 | 0.133 | 0.045 | 0.6165 | 1.00 |
| OU-III | 20 | 0.4543 | 2.7521 | 16.1472 | 2.8585 | 2.059 | 1.775 | 90.8014 | 32.69 |
| OU-III | 40 | 0.4543 | 10.1602 | 38.6211 | 6.0314 | 5.427 | 9.142 | 90.2882 | 78.20 |
| OU-III | 160 | 0.4543 | 3.4812 | 16.9449 | 1.5338 | 0.865 | 2.923 | 89.7057 | 34.31 |
| TFG | off | 0.0000 | 0.2116 | 0.5007 | 0.2983 | 0.292 | 0.049 | 0.9720 | 1.00 |
| TFG | 20 | 0.4543 | 18.0333 | 89.0081 | 11.6143 | 8.061 | 12.832 | 103.9845 | 177.77 |
| TFG | 40 | 0.4543 | 24.3674 | 96.2669 | 12.1513 | 8.002 | 11.702 | 104.1728 | 192.27 |
| TFG | 160 | 0.4543 | 2.3015 | 21.6690 | 9.8727 | 9.479 | 9.015 | 75.0771 | 43.28 |

## Sensor attribution

The nominal cruise cell rerun with the model's gyroscope terms switched
off, so the accelerometer is the only perturbed sensor.  Compare against
the 2400 rpm row of the engine-speed table.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.3249 | 0.5832 | 0.2350 | 0.218 | 0.046 | 0.7281 | 1.00 |
| OU-II | 2400 | 0.4543 | 41.1428 | 62.4963 | 3.5672 | 1.395 | 35.660 | 92.1947 | 107.16 |
| OU-III | off | 0.0000 | 0.2106 | 0.4939 | 0.1436 | 0.133 | 0.045 | 0.6165 | 1.00 |
| OU-III | 2400 | 0.4543 | 44.9216 | 69.0966 | 10.8713 | 9.493 | 40.955 | 90.0023 | 139.91 |
| TFG | off | 0.0000 | 0.2116 | 0.5007 | 0.2983 | 0.292 | 0.049 | 0.9720 | 1.00 |
| TFG | 2400 | 0.4543 | 5.3951 | 49.3361 | 8.1301 | 5.136 | 9.408 | 105.6162 | 98.54 |

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
