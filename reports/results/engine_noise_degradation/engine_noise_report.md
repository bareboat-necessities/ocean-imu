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

Source commit used for the replay: `6b20d608bca82f6487fd6d4c1d8344f4af829eef`.

## Engine speed

Vibration level fixed at 0.60 m/s² of hull broadband
RMS at 2400 rpm, sensor bandwidth 80 Hz.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2927 | 0.5152 | 0.1864 | 0.174 | 0.062 | 0.5949 | 1.00 |
| OU-II | 800 | 0.3751 | 0.7047 | 2.3448 | 2.4832 | 2.417 | 1.914 | 67.3041 | 4.55 |
| OU-II | 1200 | 0.3943 | 1.0776 | 3.9688 | 2.7652 | 2.724 | 3.495 | 67.3517 | 7.70 |
| OU-II | 1600 | 0.3810 | 1.4482 | 7.8353 | 4.6296 | 4.554 | 3.491 | 100.8458 | 15.21 |
| OU-II | 2000 | 0.4106 | 3.2175 | 14.9927 | 5.9885 | 5.894 | 4.702 | 104.1636 | 29.10 |
| OU-II | 2400 | 0.4543 | 5.9056 | 18.6099 | 4.8961 | 4.754 | 6.557 | 103.7679 | 36.12 |
| OU-II | 2800 | 0.5067 | 6.8846 | 19.0747 | 2.5891 | 2.014 | 7.043 | 103.9298 | 37.03 |
| OU-II | 3200 | 0.5657 | 6.9245 | 21.6500 | 5.2234 | 4.632 | 7.205 | 103.5643 | 42.03 |
| OU-III | off | 0.0000 | 0.1791 | 0.4052 | 0.1559 | 0.146 | 0.043 | 0.5699 | 1.00 |
| OU-III | 800 | 0.3751 | 0.7899 | 4.4617 | 2.2537 | 2.208 | 0.877 | 104.1290 | 11.01 |
| OU-III | 1200 | 0.3943 | 1.4033 | 9.0256 | 3.1859 | 3.128 | 1.402 | 103.8734 | 22.27 |
| OU-III | 1600 | 0.3810 | 2.6458 | 15.7608 | 3.9928 | 3.896 | 2.540 | 103.9108 | 38.89 |
| OU-III | 2000 | 0.4106 | 5.0602 | 20.4924 | 3.6519 | 3.428 | 4.551 | 103.8704 | 50.57 |
| OU-III | 2400 | 0.4543 | 5.7349 | 22.1951 | 1.9885 | 1.424 | 5.159 | 103.9038 | 54.77 |
| OU-III | 2800 | 0.5067 | 6.6243 | 23.1718 | 2.3765 | 1.608 | 5.839 | 103.8927 | 57.18 |
| OU-III | 3200 | 0.5657 | 9.3926 | 26.3465 | 5.3037 | 4.739 | 8.216 | 103.9416 | 65.01 |
| TFG | off | 0.0000 | 0.1805 | 0.4164 | 0.0833 | 0.064 | 0.089 | 0.5550 | 1.00 |
| TFG | 800 | 0.3751 | 0.2660 | 0.7343 | 4.5830 | 3.898 | 0.340 | 49.4267 | 1.76 |
| TFG | 1200 | 0.3943 | 0.2549 | 0.8543 | 4.3199 | 1.871 | 0.483 | 56.1184 | 2.05 |
| TFG | 1600 | 0.3810 | 0.2682 | 0.9229 | 5.4284 | 4.113 | 0.360 | 62.8588 | 2.22 |
| TFG | 2000 | 0.4106 | 0.6224 | 8.0528 | 11.9899 | 8.852 | 7.356 | 109.0520 | 19.34 |
| TFG | 2400 | 0.4543 | 0.3614 | 9.5741 | 7.5455 | 6.332 | 8.820 | 105.2866 | 22.99 |
| TFG | 2800 | 0.5067 | 2.2182 | 21.0053 | 3.5607 | 2.407 | 6.742 | 103.9191 | 50.44 |
| TFG | 3200 | 0.5657 | 5.3597 | 49.3688 | 6.7570 | 3.774 | 10.838 | 104.0951 | 118.56 |

## Vibration level

Engine speed fixed at 2400 rpm, sensor bandwidth 80 Hz.  The level is the hull broadband RMS
before the sensor's anti-alias filter, so a quiet, well-isolated
installation sits at the low end and a sensor near the engine bed at the
high end.

| Family | Hull level [m/s²] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2927 | 0.5152 | 0.1864 | 0.174 | 0.062 | 0.5949 | 1.00 |
| OU-II | 0.15 | 0.1136 | 0.2908 | 0.6537 | 2.0461 | 1.969 | 0.269 | 17.3246 | 1.27 |
| OU-II | 0.30 | 0.2272 | 0.4695 | 1.6336 | 2.5279 | 2.418 | 1.193 | 81.6288 | 3.17 |
| OU-II | 1.20 | 0.9087 | 120.3010 | 127.7845 | 11.8559 | 8.604 | 113.121 | 103.9740 | 248.05 |
| OU-II | 2.40 | 1.8173 | 391.8066 | 400.8692 | 32.4415 | 23.077 | 360.889 | 103.7307 | 778.14 |
| OU-III | off | 0.0000 | 0.1791 | 0.4052 | 0.1559 | 0.146 | 0.043 | 0.5699 | 1.00 |
| OU-III | 0.15 | 0.1136 | 0.2611 | 0.9833 | 1.9270 | 1.876 | 0.355 | 53.5014 | 2.43 |
| OU-III | 0.30 | 0.2272 | 0.6515 | 3.8651 | 1.9927 | 1.904 | 0.656 | 103.9585 | 9.54 |
| OU-III | 1.20 | 0.9087 | 117.1975 | 125.2986 | 29.2488 | 26.237 | 114.923 | 103.4224 | 309.19 |
| OU-III | 2.40 | 1.8173 | 248.6947 | 257.5318 | 36.7700 | 30.063 | 235.821 | 103.3773 | 635.49 |
| TFG | off | 0.0000 | 0.1805 | 0.4164 | 0.0833 | 0.064 | 0.089 | 0.5550 | 1.00 |
| TFG | 0.15 | 0.1136 | 0.1922 | 0.5574 | 2.6000 | 2.397 | 0.295 | 8.6578 | 1.34 |
| TFG | 0.30 | 0.2272 | 0.2244 | 1.2253 | 7.8180 | 7.415 | 0.918 | 30.3742 | 2.94 |
| TFG | 1.20 | 0.9087 | 14.6410 | 85.8477 | 15.4099 | 10.864 | 16.934 | 103.7003 | 206.16 |
| TFG | 2.40 | 1.8173 | 7.2184 | 43.8402 | 15.2515 | 11.868 | 4.158 | 103.8954 | 105.28 |

## Sensor anti-alias bandwidth

Engine speed fixed at 2400 rpm and the hull level at
0.60 m/s².  A narrower filter removes power *and*
stops the high crank orders from folding below Nyquist.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2927 | 0.5152 | 0.1864 | 0.174 | 0.062 | 0.5949 | 1.00 |
| OU-II | 20 | 0.1865 | 0.3733 | 0.8821 | 2.2187 | 2.105 | 0.608 | 55.3575 | 1.71 |
| OU-II | 40 | 0.2939 | 0.7795 | 3.0683 | 2.4796 | 2.293 | 2.737 | 88.8462 | 5.96 |
| OU-II | 160 | 0.6034 | 3.3983 | 11.9778 | 2.9199 | 2.594 | 4.490 | 103.9050 | 23.25 |
| OU-III | off | 0.0000 | 0.1791 | 0.4052 | 0.1559 | 0.146 | 0.043 | 0.5699 | 1.00 |
| OU-III | 20 | 0.1865 | 0.3212 | 1.0202 | 1.4967 | 1.404 | 0.511 | 84.2525 | 2.52 |
| OU-III | 40 | 0.2939 | 0.6702 | 3.9858 | 0.8099 | 0.529 | 0.644 | 103.9926 | 9.84 |
| OU-III | 160 | 0.6034 | 5.7565 | 23.5255 | 1.9339 | 1.245 | 5.200 | 103.8312 | 58.05 |
| TFG | off | 0.0000 | 0.1805 | 0.4164 | 0.0833 | 0.064 | 0.089 | 0.5550 | 1.00 |
| TFG | 20 | 0.1865 | 0.2017 | 0.7202 | 6.5507 | 6.141 | 0.476 | 23.7288 | 1.73 |
| TFG | 40 | 0.2939 | 0.2783 | 2.2009 | 12.5369 | 12.201 | 1.807 | 68.1824 | 5.29 |
| TFG | 160 | 0.6034 | 1.8899 | 17.7333 | 6.7937 | 5.604 | 6.979 | 106.2763 | 42.59 |

## Matched-power control

The same bandwidth sweep with the hull level rescaled so every cell
records the same vibration RMS as the nominal cruise cell.  The folded
line frequencies still differ from cell to cell; the recorded power no
longer does.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2927 | 0.5152 | 0.1864 | 0.174 | 0.062 | 0.5949 | 1.00 |
| OU-II | 20 | 0.4543 | 0.9713 | 3.9293 | 3.3301 | 3.157 | 3.513 | 88.5674 | 7.63 |
| OU-II | 40 | 0.4543 | 1.1626 | 7.1195 | 2.8410 | 2.568 | 4.232 | 103.8204 | 13.82 |
| OU-II | 160 | 0.4543 | 1.0228 | 6.2316 | 3.2134 | 3.094 | 3.743 | 104.0871 | 12.10 |
| OU-III | off | 0.0000 | 0.1791 | 0.4052 | 0.1559 | 0.146 | 0.043 | 0.5699 | 1.00 |
| OU-III | 20 | 0.4543 | 0.3940 | 2.4782 | 1.0406 | 0.848 | 0.550 | 103.7975 | 6.12 |
| OU-III | 40 | 0.4543 | 2.5267 | 13.2158 | 1.6971 | 1.190 | 2.285 | 103.9666 | 32.61 |
| OU-III | 160 | 0.4543 | 2.1107 | 11.9719 | 1.8559 | 1.575 | 1.927 | 103.8127 | 29.54 |
| TFG | off | 0.0000 | 0.1805 | 0.4164 | 0.0833 | 0.064 | 0.089 | 0.5550 | 1.00 |
| TFG | 20 | 0.4543 | 0.3432 | 2.9426 | 11.3879 | 10.874 | 2.482 | 91.1222 | 7.07 |
| TFG | 40 | 0.4543 | 0.3214 | 9.7014 | 3.2649 | 2.960 | 9.040 | 103.4354 | 23.30 |
| TFG | 160 | 0.4543 | 0.3374 | 9.5738 | 8.9964 | 8.369 | 9.030 | 97.7787 | 22.99 |

## Sensor attribution

The nominal cruise cell rerun with the model's gyroscope terms switched
off, so the accelerometer is the only perturbed sensor.  Compare against
the 2400 rpm row of the engine-speed table.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2927 | 0.5152 | 0.1864 | 0.174 | 0.062 | 0.5949 | 1.00 |
| OU-II | 2400 | 0.4543 | 5.9073 | 18.6126 | 4.8953 | 4.753 | 6.558 | 103.7727 | 36.13 |
| OU-III | off | 0.0000 | 0.1791 | 0.4052 | 0.1559 | 0.146 | 0.043 | 0.5699 | 1.00 |
| OU-III | 2400 | 0.4543 | 5.7373 | 22.1987 | 1.9880 | 1.424 | 5.160 | 103.8989 | 54.78 |
| TFG | off | 0.0000 | 0.1805 | 0.4164 | 0.0833 | 0.064 | 0.089 | 0.5550 | 1.00 |
| TFG | 2400 | 0.4543 | 0.3614 | 9.5743 | 7.5429 | 6.331 | 8.820 | 105.3151 | 22.99 |

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
