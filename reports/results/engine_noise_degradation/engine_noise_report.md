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

Source commit used for the replay: `9de0543507509cabb478b751b83bf5c74711764d`.

## Engine speed

Vibration level fixed at 0.60 m/s² of hull broadband
RMS at 2400 rpm, sensor bandwidth 80 Hz.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2927 | 0.5152 | 0.1870 | 0.174 | 0.062 | 0.5969 | 1.00 |
| OU-II | 800 | 0.3751 | 0.6883 | 2.1993 | 1.2009 | 1.077 | 1.909 | 43.7937 | 4.27 |
| OU-II | 1200 | 0.3943 | 1.0712 | 3.9223 | 2.0868 | 2.018 | 3.494 | 68.1664 | 7.61 |
| OU-II | 1600 | 0.3810 | 1.4410 | 7.7977 | 3.8550 | 3.766 | 3.489 | 86.6455 | 15.14 |
| OU-II | 2000 | 0.4106 | 3.2056 | 14.9466 | 5.8938 | 5.799 | 4.695 | 104.2585 | 29.01 |
| OU-II | 2400 | 0.4543 | 5.9023 | 18.5935 | 4.8459 | 4.711 | 6.555 | 91.9107 | 36.09 |
| OU-II | 2800 | 0.5067 | 6.8804 | 19.0529 | 2.1714 | 1.675 | 7.040 | 98.0105 | 36.98 |
| OU-II | 3200 | 0.5657 | 6.9159 | 21.6179 | 4.1167 | 3.587 | 7.199 | 97.4059 | 41.96 |
| OU-III | off | 0.0000 | 0.1791 | 0.4052 | 0.1561 | 0.146 | 0.043 | 0.5763 | 1.00 |
| OU-III | 800 | 0.3751 | 0.7816 | 4.3703 | 2.1118 | 2.055 | 0.870 | 99.7571 | 10.79 |
| OU-III | 1200 | 0.3943 | 1.4013 | 9.0080 | 3.2147 | 3.160 | 1.400 | 92.1063 | 22.23 |
| OU-III | 1600 | 0.3810 | 2.6442 | 15.7470 | 3.8791 | 3.781 | 2.538 | 90.2429 | 38.86 |
| OU-III | 2000 | 0.4106 | 5.0586 | 20.4740 | 3.6018 | 3.374 | 4.550 | 90.6529 | 50.53 |
| OU-III | 2400 | 0.4543 | 5.7341 | 22.1857 | 1.9721 | 1.411 | 5.158 | 90.2918 | 54.75 |
| OU-III | 2800 | 0.5067 | 6.6233 | 23.1593 | 2.3363 | 1.570 | 5.838 | 90.5730 | 57.16 |
| OU-III | 3200 | 0.5657 | 9.3909 | 26.3297 | 5.2249 | 4.665 | 8.214 | 91.2403 | 64.98 |
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
| OU-II | off | 0.0000 | 0.2927 | 0.5152 | 0.1870 | 0.174 | 0.062 | 0.5969 | 1.00 |
| OU-II | 0.15 | 0.1136 | 0.2906 | 0.6536 | 1.5444 | 1.444 | 0.269 | 11.5038 | 1.27 |
| OU-II | 0.30 | 0.2272 | 0.4695 | 1.6335 | 2.1423 | 2.027 | 1.193 | 63.0706 | 3.17 |
| OU-II | 1.20 | 0.9087 | 118.6358 | 126.6301 | 11.5919 | 7.778 | 112.012 | 103.8674 | 245.81 |
| OU-II | 2.40 | 1.8173 | 363.3440 | 374.5676 | 30.2517 | 21.394 | 340.469 | 103.3630 | 727.09 |
| OU-III | off | 0.0000 | 0.1791 | 0.4052 | 0.1561 | 0.146 | 0.043 | 0.5763 | 1.00 |
| OU-III | 0.15 | 0.1136 | 0.2609 | 0.9830 | 1.6535 | 1.607 | 0.354 | 48.8667 | 2.43 |
| OU-III | 0.30 | 0.2272 | 0.6512 | 3.8638 | 1.9140 | 1.828 | 0.656 | 90.1523 | 9.54 |
| OU-III | 1.20 | 0.9087 | 117.1914 | 125.2212 | 28.2149 | 25.126 | 114.877 | 103.6796 | 309.05 |
| OU-III | 2.40 | 1.8173 | 232.7603 | 240.0385 | 33.3475 | 26.975 | 223.191 | 103.3755 | 592.42 |
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
| OU-II | off | 0.0000 | 0.2927 | 0.5152 | 0.1870 | 0.174 | 0.062 | 0.5969 | 1.00 |
| OU-II | 20 | 0.1865 | 0.3726 | 0.8802 | 2.4599 | 2.394 | 0.607 | 11.7906 | 1.71 |
| OU-II | 40 | 0.2939 | 0.7794 | 3.0663 | 2.3712 | 2.187 | 2.736 | 72.3988 | 5.95 |
| OU-II | 160 | 0.6034 | 3.3964 | 11.9673 | 2.9059 | 2.576 | 4.489 | 97.8470 | 23.23 |
| OU-III | off | 0.0000 | 0.1791 | 0.4052 | 0.1561 | 0.146 | 0.043 | 0.5763 | 1.00 |
| OU-III | 20 | 0.1865 | 0.3212 | 1.0195 | 1.5171 | 1.423 | 0.511 | 66.6800 | 2.52 |
| OU-III | 40 | 0.2939 | 0.6695 | 3.9811 | 0.8055 | 0.546 | 0.644 | 90.1870 | 9.83 |
| OU-III | 160 | 0.6034 | 5.7562 | 23.5226 | 1.8617 | 1.129 | 5.200 | 90.1459 | 58.05 |
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
| OU-II | off | 0.0000 | 0.2927 | 0.5152 | 0.1870 | 0.174 | 0.062 | 0.5969 | 1.00 |
| OU-II | 20 | 0.4543 | 0.9682 | 3.9085 | 2.4069 | 2.281 | 3.512 | 73.6650 | 7.59 |
| OU-II | 40 | 0.4543 | 1.1519 | 7.0830 | 1.7965 | 1.599 | 4.230 | 90.7109 | 13.75 |
| OU-II | 160 | 0.4543 | 1.0223 | 6.2293 | 3.2727 | 3.157 | 3.743 | 90.2116 | 12.09 |
| OU-III | off | 0.0000 | 0.1791 | 0.4052 | 0.1561 | 0.146 | 0.043 | 0.5763 | 1.00 |
| OU-III | 20 | 0.4543 | 0.3905 | 2.4492 | 0.9586 | 0.774 | 0.548 | 104.4269 | 6.04 |
| OU-III | 40 | 0.4543 | 2.5251 | 13.1994 | 1.5234 | 0.983 | 2.284 | 91.2836 | 32.58 |
| OU-III | 160 | 0.4543 | 2.1105 | 11.9702 | 1.7832 | 1.488 | 1.927 | 90.0437 | 29.54 |
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
| OU-II | off | 0.0000 | 0.2927 | 0.5152 | 0.1870 | 0.174 | 0.062 | 0.5969 | 1.00 |
| OU-II | 2400 | 0.4543 | 5.9041 | 18.5962 | 4.8450 | 4.711 | 6.556 | 91.8860 | 36.10 |
| OU-III | off | 0.0000 | 0.1791 | 0.4052 | 0.1561 | 0.146 | 0.043 | 0.5763 | 1.00 |
| OU-III | 2400 | 0.4543 | 5.7364 | 22.1892 | 1.9718 | 1.410 | 5.160 | 90.2934 | 54.76 |
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
