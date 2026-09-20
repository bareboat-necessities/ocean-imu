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

Source commit used for the replay: `02adc1b0d12b34ff91b1c30952aef2f53d68db4f`.

## Engine speed

Vibration level fixed at 0.60 m/s² of hull broadband
RMS at 2400 rpm, sensor bandwidth 80 Hz.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2890 | 0.5049 | 0.1794 | 0.171 | 0.041 | 0.5879 | 1.00 |
| OU-II | 800 | 0.3751 | 0.4142 | 2.2408 | 1.6530 | 1.597 | 0.800 | 104.0191 | 4.44 |
| OU-II | 1200 | 0.3943 | 0.6314 | 3.7180 | 2.4209 | 2.387 | 0.961 | 103.8108 | 7.36 |
| OU-II | 1600 | 0.3810 | 0.8681 | 5.6343 | 2.8326 | 2.789 | 1.055 | 103.8264 | 11.16 |
| OU-II | 2000 | 0.4106 | 2.0413 | 9.7845 | 3.2798 | 3.133 | 1.875 | 103.8013 | 19.38 |
| OU-II | 2400 | 0.4543 | 3.9001 | 13.2904 | 2.5900 | 2.310 | 3.411 | 103.8087 | 26.32 |
| OU-II | 2800 | 0.5067 | 4.8308 | 14.7896 | 1.3868 | 0.616 | 4.217 | 103.8788 | 29.29 |
| OU-II | 3200 | 0.5657 | 4.1819 | 15.4524 | 2.5478 | 2.183 | 3.784 | 104.0129 | 30.60 |
| OU-III | off | 0.0000 | 0.1792 | 0.4047 | 0.1293 | 0.117 | 0.035 | 0.5136 | 1.00 |
| OU-III | 800 | 0.3751 | 0.7708 | 3.9014 | 2.0738 | 2.017 | 0.841 | 99.2937 | 9.64 |
| OU-III | 1200 | 0.3943 | 1.3530 | 7.8580 | 3.1374 | 3.082 | 1.351 | 91.9878 | 19.42 |
| OU-III | 1600 | 0.3810 | 2.1632 | 12.7705 | 3.6403 | 3.564 | 2.083 | 90.2810 | 31.56 |
| OU-III | 2000 | 0.4106 | 4.3613 | 18.1898 | 3.4325 | 3.220 | 3.962 | 90.6552 | 44.95 |
| OU-III | 2400 | 0.4543 | 4.8645 | 19.7530 | 1.9236 | 1.390 | 4.394 | 90.3037 | 48.81 |
| OU-III | 2800 | 0.5067 | 4.8701 | 19.2356 | 1.8201 | 1.045 | 4.339 | 90.4153 | 47.53 |
| OU-III | 3200 | 0.5657 | 5.7406 | 20.3989 | 3.8852 | 3.396 | 5.068 | 91.1547 | 50.40 |
| TFG | off | 0.0000 | 0.1804 | 0.4169 | 0.0693 | 0.043 | 0.090 | 0.5825 | 1.00 |
| TFG | 800 | 0.3751 | 0.2738 | 0.8626 | 5.6856 | 5.493 | 0.435 | 55.8642 | 2.07 |
| TFG | 1200 | 0.3943 | 0.2620 | 0.8726 | 4.7760 | 4.116 | 0.371 | 61.2541 | 2.09 |
| TFG | 1600 | 0.3810 | 0.2676 | 1.3649 | 7.0948 | 6.043 | 0.612 | 73.0643 | 3.27 |
| TFG | 2000 | 0.4106 | 0.4360 | 9.6337 | 10.2204 | 9.026 | 9.053 | 104.6522 | 23.11 |
| TFG | 2400 | 0.4543 | 0.6141 | 9.5405 | 7.4100 | 6.088 | 7.917 | 104.8840 | 22.89 |
| TFG | 2800 | 0.5067 | 3.5952 | 33.0472 | 4.2440 | 2.130 | 8.235 | 104.3124 | 79.27 |
| TFG | 3200 | 0.5657 | 8.9438 | 66.1066 | 9.0656 | 5.245 | 10.372 | 104.1843 | 158.57 |

## Vibration level

Engine speed fixed at 2400 rpm, sensor bandwidth 80 Hz.  The level is the hull broadband RMS
before the sensor's anti-alias filter, so a quiet, well-isolated
installation sits at the low end and a sensor near the engine bed at the
high end.

| Family | Hull level [m/s²] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2890 | 0.5049 | 0.1794 | 0.171 | 0.041 | 0.5879 | 1.00 |
| OU-II | 0.15 | 0.1136 | 0.2851 | 0.7974 | 0.6133 | 0.514 | 0.281 | 65.1224 | 1.58 |
| OU-II | 0.30 | 0.2272 | 0.4405 | 2.4352 | 1.9778 | 1.863 | 0.625 | 104.3489 | 4.82 |
| OU-II | 1.20 | 0.9087 | 119.3794 | 141.2540 | 7.2844 | 4.154 | 109.691 | 104.0487 | 279.75 |
| OU-II | 2.40 | 1.8173 | 281.4343 | 288.0770 | 30.6351 | 25.256 | 274.030 | 103.2341 | 570.54 |
| OU-III | off | 0.0000 | 0.1792 | 0.4047 | 0.1293 | 0.117 | 0.035 | 0.5136 | 1.00 |
| OU-III | 0.15 | 0.1136 | 0.2609 | 0.9702 | 1.6485 | 1.602 | 0.331 | 49.1931 | 2.40 |
| OU-III | 0.30 | 0.2272 | 0.6418 | 3.5172 | 1.8882 | 1.801 | 0.631 | 90.0930 | 8.69 |
| OU-III | 1.20 | 0.9087 | 121.8994 | 127.9366 | 28.0020 | 25.594 | 118.940 | 103.9212 | 316.13 |
| OU-III | 2.40 | 1.8173 | 224.5525 | 231.8763 | 35.6877 | 30.428 | 217.392 | 103.5018 | 572.96 |
| TFG | off | 0.0000 | 0.1804 | 0.4169 | 0.0693 | 0.043 | 0.090 | 0.5825 | 1.00 |
| TFG | 0.15 | 0.1136 | 0.1938 | 0.5911 | 2.8397 | 2.623 | 0.332 | 9.3421 | 1.42 |
| TFG | 0.30 | 0.2272 | 0.2297 | 1.3902 | 8.4056 | 7.986 | 1.065 | 32.6559 | 3.33 |
| TFG | 1.20 | 0.9087 | 26.0578 | 112.5887 | 19.5788 | 13.843 | 20.904 | 103.5013 | 270.07 |
| TFG | 2.40 | 1.8173 | 7.8899 | 48.2449 | 16.3615 | 13.350 | 10.623 | 103.8668 | 115.73 |

## Sensor anti-alias bandwidth

Engine speed fixed at 2400 rpm and the hull level at
0.60 m/s².  A narrower filter removes power *and*
stops the high crank orders from folding below Nyquist.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2890 | 0.5049 | 0.1794 | 0.171 | 0.041 | 0.5879 | 1.00 |
| OU-II | 20 | 0.1865 | 0.2836 | 0.7478 | 1.4507 | 1.311 | 0.402 | 34.1251 | 1.48 |
| OU-II | 40 | 0.2939 | 0.3965 | 2.2380 | 0.9451 | 0.801 | 0.669 | 104.6273 | 4.43 |
| OU-II | 160 | 0.6034 | 2.4417 | 11.3922 | 1.6713 | 0.988 | 2.123 | 104.0589 | 22.56 |
| OU-III | off | 0.0000 | 0.1792 | 0.4047 | 0.1293 | 0.117 | 0.035 | 0.5136 | 1.00 |
| OU-III | 20 | 0.1865 | 0.3214 | 0.9884 | 1.5109 | 1.418 | 0.460 | 66.7209 | 2.44 |
| OU-III | 40 | 0.2939 | 0.6615 | 3.6447 | 0.8024 | 0.545 | 0.613 | 90.2023 | 9.01 |
| OU-III | 160 | 0.6034 | 5.0245 | 21.0117 | 1.8682 | 1.157 | 4.514 | 90.1981 | 51.92 |
| TFG | off | 0.0000 | 0.1804 | 0.4169 | 0.0693 | 0.043 | 0.090 | 0.5825 | 1.00 |
| TFG | 20 | 0.1865 | 0.2090 | 0.8721 | 7.4914 | 7.041 | 0.620 | 27.4444 | 2.09 |
| TFG | 40 | 0.2939 | 0.3191 | 2.6817 | 13.2117 | 12.923 | 2.231 | 76.7433 | 6.43 |
| TFG | 160 | 0.6034 | 2.7038 | 23.7931 | 7.3491 | 6.742 | 7.107 | 101.1553 | 57.07 |

## Matched-power control

The same bandwidth sweep with the hull level rescaled so every cell
records the same vibration RMS as the nominal cruise cell.  The folded
line frequencies still differ from cell to cell; the recorded power no
longer does.

| Family | Bandwidth [Hz] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2890 | 0.5049 | 0.1794 | 0.171 | 0.041 | 0.5879 | 1.00 |
| OU-II | 20 | 0.4543 | 0.2976 | 1.6881 | 1.5148 | 1.404 | 0.775 | 103.9128 | 3.34 |
| OU-II | 40 | 0.4543 | 1.2073 | 6.8226 | 1.0523 | 0.690 | 1.200 | 103.9424 | 13.51 |
| OU-II | 160 | 0.4543 | 0.9499 | 5.2654 | 1.4416 | 1.117 | 0.758 | 103.8940 | 10.43 |
| OU-III | off | 0.0000 | 0.1792 | 0.4047 | 0.1293 | 0.117 | 0.035 | 0.5136 | 1.00 |
| OU-III | 20 | 0.4543 | 0.3900 | 2.3237 | 0.9526 | 0.764 | 0.510 | 104.4088 | 5.74 |
| OU-III | 40 | 0.4543 | 2.3143 | 11.6680 | 1.4241 | 0.864 | 2.108 | 91.1990 | 28.83 |
| OU-III | 160 | 0.4543 | 1.9631 | 10.6534 | 1.7811 | 1.489 | 1.800 | 90.0535 | 26.32 |
| TFG | off | 0.0000 | 0.1804 | 0.4169 | 0.0693 | 0.043 | 0.090 | 0.5825 | 1.00 |
| TFG | 20 | 0.4543 | 0.5155 | 4.9095 | 9.2257 | 8.180 | 4.256 | 104.7454 | 11.78 |
| TFG | 40 | 0.4543 | 0.7339 | 10.2450 | 3.0585 | 2.614 | 7.568 | 103.9372 | 24.58 |
| TFG | 160 | 0.4543 | 0.3499 | 9.5259 | 9.9350 | 8.403 | 8.944 | 93.3385 | 22.85 |

## Sensor attribution

The nominal cruise cell rerun with the model's gyroscope terms switched
off, so the accelerometer is the only perturbed sensor.  Compare against
the 2400 rpm row of the engine-speed table.

| Family | Engine speed [rpm] | Recorded vib. [m/s²] | Z [m] | 3-D [m] | Pitch RMS [deg] | Pitch offset [deg] | 3-D offset [m] | Yaw [deg] | 3-D / baseline |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | off | 0.0000 | 0.2890 | 0.5049 | 0.1794 | 0.171 | 0.041 | 0.5879 | 1.00 |
| OU-II | 2400 | 0.4543 | 3.9010 | 13.2921 | 2.5894 | 2.309 | 3.412 | 103.8060 | 26.33 |
| OU-III | off | 0.0000 | 0.1792 | 0.4047 | 0.1293 | 0.117 | 0.035 | 0.5136 | 1.00 |
| OU-III | 2400 | 0.4543 | 4.8664 | 19.7558 | 1.9233 | 1.389 | 4.395 | 90.3046 | 48.82 |
| TFG | off | 0.0000 | 0.1804 | 0.4169 | 0.0693 | 0.043 | 0.090 | 0.5825 | 1.00 |
| TFG | 2400 | 0.4543 | 0.6147 | 9.5440 | 7.4058 | 6.085 | 7.917 | 104.9288 | 22.89 |

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
