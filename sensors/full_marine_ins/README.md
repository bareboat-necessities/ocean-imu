# Full marine INS examples

[ocean-imu](../../README.md) / [sensors](../README.md) / **full_marine_ins**

<p align="center">
  <img src="../../img/devices/AtomS3R_device.svg" width="300" alt="M5Stack AtomS3R">
</p>

These are the most complete AtomS3R applications in the repository. They combine calibrated IMU input, attitude/heading estimation, marine-motion processing, on-device UI, and serial/NMEA output.

They are intended for **marine motion estimation**, where accelerometer measurements contain both gravity and wave-induced translational acceleration and the system may need to initialize while already moving.

## Available estimators

| Example | Method | Main character |
| --- | --- | --- |
| [atomS3R_ins_kalman_ou2](atomS3R_ins_kalman_ou2/README.md) | Adaptive OU-II Kalman INS | More direct integral-drift correction; responsive sea-state adaptation |
| [atomS3R_ins_kalman_ou3](atomS3R_ins_kalman_ou3/README.md) | Adaptive OU-III Kalman INS | Higher-order integral regularization; main 3-D navigation implementation |
| [atomS3R_ins_pii_observer](atomS3R_ins_pii_observer/README.md) | Adaptive PII observer + Mahony | Much lighter computation; primarily vertical/heave-oriented motion estimation |

The OU-II and OU-III sketches learn tilt during startup, perform a one-shot magnetic north lock in the sea-state fusion filter, and then use filter yaw as the primary heading output. The PII application uses Mahony attitude with a lightweight adaptive vertical observer.

## Installation and bring-up

Follow the common [Arduino installation instructions](../README.md#arduino-installation). Before debugging one of these larger applications, first verify the same board with [Basic IMU](../imu_basic/README.md) and preferably one of the [Compass / AHRS](../compass_ahrs/README.md) examples. That separates sensor/calibration problems from full estimator behavior.

The full sketches enable the shared IMU calibration wizard by default through `SEA_STATE_ENABLE_WIZARD=1`. They also default to graphical UI and NMEA serial output; those behaviors can be changed with the compile-time switches near the top of each sketch.

## Which should I use?

Start with **OU-III** if you want the project’s primary full 3-D marine navigation estimator. Try **OU-II** when you prefer its more direct drift regularization and adaptation response. Use the **PII observer** when MCU cost and simplicity matter more than the full Kalman-state solution.

## Related material

The mathematical descriptions and validation links are collected in the [project README](../../README.md). Source implementations live under [`../../src/kalman_ou_ii/`](../../src/kalman_ou_ii/), [`../../src/kalman_ou_iii/`](../../src/kalman_ou_iii/), and [`../../src/pii_observer/`](../../src/pii_observer/).

[Back to all sensor examples](../README.md)
