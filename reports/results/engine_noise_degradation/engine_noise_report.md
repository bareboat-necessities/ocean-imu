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

Source commit used for the replay: `1cb65a24379b74bf59443701c9ecb3be03f0574d`.

## Engine speed

Vibration level fixed at 0.60 m/s² of hull broadband
RMS at 2400 rpm, sensor bandwidth 80 Hz.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2922 | 0.5115 | 0.1888 | 0.173 | 0.061 | 0.5994 | 1.00 |
| OU-II | 800 | 0.3751 | 0.6977 | 2.2947 | 2.4782 | 2.412 | 1.877 | 67.1905 | 4.49 |
| OU-II | 1200 | 0.3943 | 1.0668 | 3.8905 | 2.7647 | 2.724 | 3.431 | 67.3281 | 7.61 |
| OU-II | 1600 | 0.3810 | 1.4189 | 7.6736 | 4.6093 | 4.533 | 3.460 | 100.2663 | 15.00 |
| OU-II | 2000 | 0.4106 | 3.1594 | 14.6213 | 5.9740 | 5.880 | 4.620 | 104.1480 | 28.59 |
| OU-II | 2400 | 0.4543 | 5.7983 | 18.1659 | 4.8985 | 4.758 | 6.440 | 103.7712 | 35.52 |
| OU-II | 2800 | 0.5067 | 6.7552 | 18.6366 | 2.5678 | 2.003 | 6.911 | 103.9569 | 36.44 |
| OU-II | 3200 | 0.5657 | 6.7272 | 21.1169 | 5.0940 | 4.503 | 7.028 | 103.7110 | 41.28 |
| OU-III | off | 0.0000 | 0.1790 | 0.4046 | 0.1661 | 0.157 | 0.047 | 0.5712 | 1.00 |
| OU-III | 800 | 0.3751 | 0.8348 | 5.2136 | 2.7570 | 2.699 | 1.118 | 104.0807 | 12.88 |
| OU-III | 1200 | 0.3943 | 1.5409 | 10.8662 | 3.8958 | 3.822 | 1.647 | 103.7307 | 26.85 |
| OU-III | 1600 | 0.3810 | 2.9830 | 19.7699 | 4.9310 | 4.837 | 2.935 | 103.9573 | 48.86 |
| OU-III | 2000 | 0.4106 | 6.3478 | 25.6342 | 3.8575 | 3.627 | 5.616 | 103.8580 | 63.35 |
| OU-III | 2400 | 0.4543 | 7.8390 | 27.2576 | 2.0921 | 1.166 | 6.775 | 103.8768 | 67.36 |
| OU-III | 2800 | 0.5067 | 8.4342 | 27.6871 | 3.3934 | 2.566 | 7.174 | 103.8182 | 68.42 |
| OU-III | 3200 | 0.5657 | 11.3565 | 31.3423 | 6.9313 | 6.295 | 9.745 | 103.8913 | 77.46 |
| TFG | off | 0.0000 | 0.1805 | 0.4164 | 0.0833 | 0.064 | 0.089 | 0.5549 | 1.00 |
| TFG | 800 | 0.3751 | 0.2660 | 0.7343 | 4.5834 | 3.898 | 0.340 | 49.4251 | 1.76 |
| TFG | 1200 | 0.3943 | 0.2549 | 0.8543 | 4.3199 | 1.871 | 0.483 | 56.1186 | 2.05 |
| TFG | 1600 | 0.3810 | 0.2682 | 0.9229 | 5.4286 | 4.113 | 0.360 | 62.8584 | 2.22 |
| TFG | 2000 | 0.4106 | 0.6225 | 8.0555 | 12.0031 | 8.867 | 7.359 | 109.0165 | 19.35 |
| TFG | 2400 | 0.4543 | 0.3614 | 9.5740 | 7.5454 | 6.332 | 8.820 | 105.2877 | 22.99 |
| TFG | 2800 | 0.5067 | 2.2182 | 21.0043 | 3.5607 | 2.407 | 6.742 | 103.9179 | 50.44 |
| TFG | 3200 | 0.5657 | 5.3616 | 49.3773 | 6.7576 | 3.775 | 10.868 | 104.0925 | 118.58 |

## Vibration level

Engine speed fixed at 2400 rpm, sensor bandwidth 80 Hz.  The level is the hull broadband RMS
before the sensor's anti-alias filter, so a quiet, well-isolated
installation sits at the low end and a sensor near the engine bed at the
high end.

| Family | Hull level [m/s²] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2922 | 0.5115 | 0.1888 | 0.173 | 0.061 | 0.5994 | 1.00 |
| OU-II | 0.15 | 0.1136 | 0.2896 | 0.6466 | 2.0452 | 1.968 | 0.264 | 17.3293 | 1.26 |
| OU-II | 0.30 | 0.2272 | 0.4653 | 1.6042 | 2.5248 | 2.415 | 1.169 | 81.6023 | 3.14 |
| OU-II | 1.20 | 0.9087 | 113.7468 | 123.9065 | 9.7480 | 5.368 | 106.401 | 103.9185 | 242.24 |
| OU-II | 2.40 | 1.8173 | 320.0658 | 332.0988 | 32.2941 | 22.999 | 301.034 | 103.4830 | 649.27 |
| OU-III | off | 0.0000 | 0.1790 | 0.4046 | 0.1661 | 0.157 | 0.047 | 0.5712 | 1.00 |
| OU-III | 0.15 | 0.1136 | 0.2634 | 1.0340 | 1.9523 | 1.908 | 0.448 | 53.6154 | 2.56 |
| OU-III | 0.30 | 0.2272 | 0.6929 | 4.5247 | 2.2195 | 2.129 | 0.925 | 103.9165 | 11.18 |
| OU-III | 1.20 | 0.9087 | 104.2448 | 112.3932 | 31.0379 | 27.927 | 101.934 | 103.5991 | 277.76 |
| OU-III | 2.40 | 1.8173 | 252.3653 | 264.5747 | 35.5550 | 28.631 | 237.799 | 103.4869 | 653.85 |
| TFG | off | 0.0000 | 0.1805 | 0.4164 | 0.0833 | 0.064 | 0.089 | 0.5549 | 1.00 |
| TFG | 0.15 | 0.1136 | 0.1922 | 0.5574 | 2.5999 | 2.397 | 0.295 | 8.6576 | 1.34 |
| TFG | 0.30 | 0.2272 | 0.2244 | 1.2253 | 7.8180 | 7.415 | 0.918 | 30.3745 | 2.94 |
| TFG | 1.20 | 0.9087 | 14.6427 | 85.8522 | 15.4097 | 10.863 | 16.968 | 103.6998 | 206.18 |
| TFG | 2.40 | 1.8173 | 7.2179 | 43.8339 | 15.2511 | 11.868 | 4.143 | 103.8949 | 105.27 |

## Sensor anti-alias bandwidth

Engine speed fixed at 2400 rpm and the hull level at
0.60 m/s².  A narrower filter removes power *and*
stops the high crank orders from folding below Nyquist.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2922 | 0.5115 | 0.1888 | 0.173 | 0.061 | 0.5994 | 1.00 |
| OU-II | 20 | 0.1865 | 0.3706 | 0.8694 | 2.2175 | 2.104 | 0.596 | 55.3711 | 1.70 |
| OU-II | 40 | 0.2939 | 0.7713 | 3.0096 | 2.4749 | 2.289 | 2.683 | 88.8562 | 5.88 |
| OU-II | 160 | 0.6034 | 3.3521 | 11.7207 | 2.9293 | 2.607 | 4.445 | 103.9218 | 22.91 |
| OU-III | off | 0.0000 | 0.1790 | 0.4046 | 0.1661 | 0.157 | 0.047 | 0.5712 | 1.00 |
| OU-III | 20 | 0.1865 | 0.3248 | 1.1086 | 1.5612 | 1.462 | 0.647 | 84.4592 | 2.74 |
| OU-III | 40 | 0.2939 | 0.7133 | 4.7066 | 0.7982 | 0.496 | 0.916 | 103.9999 | 11.63 |
| OU-III | 160 | 0.6034 | 7.7579 | 27.7936 | 2.2568 | 1.224 | 6.739 | 103.7262 | 68.69 |
| TFG | off | 0.0000 | 0.1805 | 0.4164 | 0.0833 | 0.064 | 0.089 | 0.5549 | 1.00 |
| TFG | 20 | 0.1865 | 0.2017 | 0.7201 | 6.5506 | 6.141 | 0.476 | 23.7285 | 1.73 |
| TFG | 40 | 0.2939 | 0.2783 | 2.2009 | 12.5367 | 12.201 | 1.807 | 68.1819 | 5.29 |
| TFG | 160 | 0.6034 | 1.8901 | 17.7335 | 6.7936 | 5.604 | 6.980 | 106.2768 | 42.59 |

## Matched-power control

The same bandwidth sweep with the hull level rescaled so every cell
records the same vibration RMS as the nominal cruise cell.  The folded
line frequencies still differ from cell to cell; the recorded power no
longer does.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2922 | 0.5115 | 0.1888 | 0.173 | 0.061 | 0.5994 | 1.00 |
| OU-II | 20 | 0.4543 | 0.9608 | 3.8528 | 3.3082 | 3.136 | 3.444 | 88.5228 | 7.53 |
| OU-II | 40 | 0.4543 | 1.1445 | 6.9516 | 2.8028 | 2.534 | 4.151 | 103.7763 | 13.59 |
| OU-II | 160 | 0.4543 | 1.0092 | 6.1170 | 3.2038 | 3.085 | 3.710 | 104.0668 | 11.96 |
| OU-III | off | 0.0000 | 0.1790 | 0.4046 | 0.1661 | 0.157 | 0.047 | 0.5712 | 1.00 |
| OU-III | 20 | 0.4543 | 0.4113 | 2.8371 | 1.1381 | 0.899 | 0.872 | 103.8218 | 7.01 |
| OU-III | 40 | 0.4543 | 2.9320 | 16.5942 | 2.2646 | 1.666 | 2.718 | 103.8956 | 41.01 |
| OU-III | 160 | 0.4543 | 2.4262 | 14.9307 | 2.1125 | 1.801 | 2.284 | 103.7347 | 36.90 |
| TFG | off | 0.0000 | 0.1805 | 0.4164 | 0.0833 | 0.064 | 0.089 | 0.5549 | 1.00 |
| TFG | 20 | 0.4543 | 0.3432 | 2.9427 | 11.3878 | 10.874 | 2.482 | 91.1242 | 7.07 |
| TFG | 40 | 0.4543 | 0.3213 | 9.7010 | 3.2648 | 2.960 | 9.040 | 103.4329 | 23.30 |
| TFG | 160 | 0.4543 | 0.3374 | 9.5737 | 8.9963 | 8.369 | 9.030 | 97.7829 | 22.99 |

## Sensor attribution

The nominal cruise cell rerun with the model's gyroscope terms switched
off, so the accelerometer is the only perturbed sensor.  Compare against
the 2400 rpm row of the engine-speed table.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2922 | 0.5115 | 0.1888 | 0.173 | 0.061 | 0.5994 | 1.00 |
| OU-II | 2400 | 0.4543 | 5.8000 | 18.1685 | 4.8977 | 4.757 | 6.441 | 103.7714 | 35.52 |
| OU-III | off | 0.0000 | 0.1790 | 0.4046 | 0.1661 | 0.157 | 0.047 | 0.5712 | 1.00 |
| OU-III | 2400 | 0.4543 | 7.8428 | 27.2629 | 2.0919 | 1.165 | 6.778 | 103.8804 | 67.38 |
| TFG | off | 0.0000 | 0.1805 | 0.4164 | 0.0833 | 0.064 | 0.089 | 0.5549 | 1.00 |
| TFG | 2400 | 0.4543 | 0.3613 | 9.5742 | 7.5430 | 6.331 | 8.819 | 105.3145 | 22.99 |

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
