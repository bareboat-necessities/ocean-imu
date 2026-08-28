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

Source commit used for the replay: `a86d6cf7d88aca32e9b8e1298f698fe9eb6739cc`.

## Engine speed

Vibration level fixed at 0.60 m/s² of hull broadband
RMS at 2400 rpm, sensor bandwidth 80 Hz.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2794 | 0.6423 | 0.1965 | 0.162 | 0.178 | 0.6669 | 1.00 |
| OU-II | 800 | 0.3751 | 0.2969 | 0.8021 | 0.2509 | 0.203 | 0.272 | 5.1500 | 1.25 |
| OU-II | 1200 | 0.3943 | 0.7224 | 2.2746 | 2.3048 | 2.289 | 2.019 | 5.2214 | 3.54 |
| OU-II | 1600 | 0.3810 | 0.7007 | 2.3793 | 1.9178 | 1.887 | 2.012 | 10.4929 | 3.70 |
| OU-II | 2000 | 0.4106 | 0.4756 | 1.7395 | 1.4937 | 1.294 | 1.005 | 48.2908 | 2.71 |
| OU-II | 2400 | 0.4543 | 0.7404 | 3.3669 | 2.4362 | 2.304 | 2.859 | 58.5476 | 5.24 |
| OU-II | 2800 | 0.5067 | 1.2783 | 6.2225 | 2.7312 | 2.656 | 5.729 | 60.4988 | 9.69 |
| OU-II | 3200 | 0.5657 | 1.8904 | 10.5055 | 3.1172 | 3.084 | 10.014 | 58.9536 | 16.36 |
| OU-III | off | 0.0000 | 0.1791 | 0.5224 | 0.1582 | 0.144 | 0.140 | 0.6292 | 1.00 |
| OU-III | 800 | 0.3751 | 0.3275 | 0.8770 | 2.0752 | 2.014 | 0.511 | 12.6122 | 1.68 |
| OU-III | 1200 | 0.3943 | 0.3873 | 1.2680 | 1.7011 | 1.601 | 0.915 | 10.8525 | 2.43 |
| OU-III | 1600 | 0.3810 | 0.3825 | 1.4021 | 1.3930 | 1.211 | 0.905 | 17.8497 | 2.68 |
| OU-III | 2000 | 0.4106 | 0.5333 | 2.2381 | 2.7568 | 2.746 | 1.646 | 36.5903 | 4.28 |
| OU-III | 2400 | 0.4543 | 0.8284 | 4.2341 | 3.0843 | 3.075 | 3.874 | 45.3828 | 8.10 |
| OU-III | 2800 | 0.5067 | 1.1656 | 6.7475 | 3.2899 | 3.280 | 6.270 | 54.1010 | 12.92 |
| OU-III | 3200 | 0.5657 | 1.4006 | 8.6774 | 3.4296 | 3.286 | 8.306 | 67.0417 | 16.61 |
| TFG | off | 0.0000 | 0.1897 | 0.7446 | 0.3860 | 0.353 | 0.206 | 0.8042 | 1.00 |
| TFG | 800 | 0.3751 | 0.3738 | 1.4042 | 3.5404 | 3.471 | 0.980 | 9.6493 | 1.89 |
| TFG | 1200 | 0.3943 | 0.6071 | 2.2685 | 5.4900 | 5.343 | 1.814 | 13.1057 | 3.05 |
| TFG | 1600 | 0.3810 | 0.3974 | 1.7602 | 4.3086 | 4.190 | 1.080 | 11.3392 | 2.36 |
| TFG | 2000 | 0.4106 | 1.0255 | 5.1995 | 12.0647 | 11.947 | 4.470 | 38.9851 | 6.98 |
| TFG | 2400 | 0.4543 | 2.9348 | 36.5099 | 16.1051 | 15.056 | 35.528 | 77.6807 | 49.03 |
| TFG | 2800 | 0.5067 | 13.2731 | 107.5711 | 15.3844 | 9.653 | 104.580 | 108.7273 | 144.47 |
| TFG | 3200 | 0.5657 | 14.1098 | 108.5562 | 14.8532 | 5.205 | 105.271 | 122.0850 | 145.80 |

## Vibration level

Engine speed fixed at 2400 rpm, sensor bandwidth 80 Hz.  The level is the hull broadband RMS
before the sensor's anti-alias filter, so a quiet, well-isolated
installation sits at the low end and a sensor near the engine bed at the
high end.

| Family | Hull level [m/s²] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2794 | 0.6423 | 0.1965 | 0.162 | 0.178 | 0.6669 | 1.00 |
| OU-II | 0.15 | 0.1136 | 0.2731 | 0.6906 | 1.4464 | 1.333 | 0.289 | 5.7188 | 1.08 |
| OU-II | 0.30 | 0.2272 | 0.3421 | 1.0068 | 2.3086 | 2.246 | 0.637 | 13.7228 | 1.57 |
| OU-II | 1.20 | 0.9087 | 6.3050 | 38.0500 | 9.9520 | 9.542 | 15.631 | 103.7336 | 59.24 |
| OU-II | 2.40 | 1.8173 | 144.7942 | 165.9562 | 17.0174 | 12.453 | 137.890 | 104.2149 | 258.38 |
| OU-III | off | 0.0000 | 0.1791 | 0.5224 | 0.1582 | 0.144 | 0.140 | 0.6292 | 1.00 |
| OU-III | 0.15 | 0.1136 | 0.2721 | 0.7073 | 2.1095 | 1.955 | 0.403 | 7.3832 | 1.35 |
| OU-III | 0.30 | 0.2272 | 0.4392 | 1.3722 | 2.7589 | 2.745 | 1.153 | 15.2141 | 2.63 |
| OU-III | 1.20 | 0.9087 | 6.1227 | 37.5749 | 8.1383 | 7.722 | 10.650 | 103.8804 | 71.92 |
| OU-III | 2.40 | 1.8173 | 85.1608 | 124.3681 | 12.1334 | 5.188 | 80.611 | 103.9427 | 238.06 |
| TFG | off | 0.0000 | 0.1897 | 0.7446 | 0.3860 | 0.353 | 0.206 | 0.8042 | 1.00 |
| TFG | 0.15 | 0.1136 | 0.3264 | 1.1538 | 3.7156 | 3.559 | 0.842 | 9.6125 | 1.55 |
| TFG | 0.30 | 0.2272 | 0.8142 | 3.2680 | 10.0955 | 9.888 | 2.871 | 28.7707 | 4.39 |
| TFG | 1.20 | 0.9087 | 14.0772 | 109.8450 | 17.9422 | 10.080 | 106.216 | 118.2798 | 147.53 |
| TFG | 2.40 | 1.8173 | 30.7412 | 121.7265 | 17.3391 | 7.936 | 115.340 | 103.7506 | 163.48 |

## Sensor anti-alias bandwidth

Engine speed fixed at 2400 rpm and the hull level at
0.60 m/s².  A narrower filter removes power *and*
stops the high crank orders from folding below Nyquist.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2794 | 0.6423 | 0.1965 | 0.162 | 0.178 | 0.6669 | 1.00 |
| OU-II | 20 | 0.1865 | 0.3302 | 0.8687 | 2.2560 | 2.145 | 0.566 | 7.3297 | 1.35 |
| OU-II | 40 | 0.2939 | 0.5647 | 1.8942 | 2.7954 | 2.787 | 1.675 | 16.0374 | 2.95 |
| OU-II | 160 | 0.6034 | 1.5378 | 8.2359 | 3.2360 | 3.210 | 7.629 | 70.6166 | 12.82 |
| OU-III | off | 0.0000 | 0.1791 | 0.5224 | 0.1582 | 0.144 | 0.140 | 0.6292 | 1.00 |
| OU-III | 20 | 0.1865 | 0.3013 | 0.7913 | 2.3452 | 2.249 | 0.506 | 8.2050 | 1.51 |
| OU-III | 40 | 0.2939 | 0.5938 | 2.1990 | 2.9224 | 2.919 | 1.995 | 17.6981 | 4.21 |
| OU-III | 160 | 0.6034 | 1.1851 | 7.5916 | 2.6912 | 2.350 | 6.805 | 85.9540 | 14.53 |
| TFG | off | 0.0000 | 0.1897 | 0.7446 | 0.3860 | 0.353 | 0.206 | 0.8042 | 1.00 |
| TFG | 20 | 0.1865 | 0.4532 | 1.5739 | 6.4024 | 6.196 | 1.270 | 17.2385 | 2.11 |
| TFG | 40 | 0.2939 | 1.2437 | 6.2387 | 15.0080 | 14.870 | 5.641 | 52.8639 | 8.38 |
| TFG | 160 | 0.6034 | 18.9119 | 109.5622 | 16.1126 | 12.685 | 106.651 | 96.4791 | 147.15 |

## Matched-power control

The same bandwidth sweep with the hull level rescaled so every cell
records the same vibration RMS as the nominal cruise cell.  The folded
line frequencies still differ from cell to cell; the recorded power no
longer does.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2794 | 0.6423 | 0.1965 | 0.162 | 0.178 | 0.6669 | 1.00 |
| OU-II | 20 | 0.4543 | 0.6827 | 2.4721 | 2.9178 | 2.913 | 2.285 | 15.8411 | 3.85 |
| OU-II | 40 | 0.4543 | 0.9674 | 4.2134 | 2.9254 | 2.917 | 3.906 | 31.1009 | 6.56 |
| OU-II | 160 | 0.4543 | 1.0584 | 4.8045 | 2.9817 | 2.965 | 4.446 | 46.7712 | 7.48 |
| OU-III | off | 0.0000 | 0.1791 | 0.5224 | 0.1582 | 0.144 | 0.140 | 0.6292 | 1.00 |
| OU-III | 20 | 0.4543 | 0.5823 | 2.2978 | 2.9186 | 2.913 | 2.107 | 19.8623 | 4.40 |
| OU-III | 40 | 0.4543 | 0.9159 | 4.5446 | 3.1597 | 3.150 | 4.283 | 33.3620 | 8.70 |
| OU-III | 160 | 0.4543 | 0.7678 | 3.7852 | 3.0535 | 3.049 | 3.501 | 34.8566 | 7.25 |
| TFG | off | 0.0000 | 0.1897 | 0.7446 | 0.3860 | 0.353 | 0.206 | 0.8042 | 1.00 |
| TFG | 20 | 0.4543 | 1.3654 | 7.9128 | 15.6754 | 15.573 | 7.182 | 62.1918 | 10.63 |
| TFG | 40 | 0.4543 | 21.9375 | 108.8916 | 15.7806 | 11.840 | 106.197 | 102.1950 | 146.25 |
| TFG | 160 | 0.4543 | 1.5942 | 11.3727 | 16.3510 | 16.275 | 10.478 | 64.6342 | 15.27 |

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
