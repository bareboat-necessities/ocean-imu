# AtomS3R qMEKF compass

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [compass_ahrs](../README.md) / **atomS3R_compass_qmekf**

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="240" alt="M5Stack AtomS3R">
</p>

This sketch implements a tilt-compensated AtomS3R compass using the project’s **Quaternion MEKF (qMEKF)** attitude estimator and the same shared AtomS3R application/calibration framework as the Mahony example.

Open the sketch: [`atomS3R_compass_qmekf.ino`](atomS3R_compass_qmekf.ino).

## What it provides

- quaternion attitude estimation with explicit accelerometer, gyro, and magnetometer noise parameters;
- initialization from accelerometer plus magnetometer when magnetic data are available, otherwise accelerometer-only initialization;
- accelerometer correction each filter cycle and fresh-magnetometer updates when available;
- graphical compass UI by default;
- NMEA 0183 serial output by default through the shared compass application;
- the common AtomS3R IMU calibration workflow.

Use this version when you want the attitude solution to be expressed as a probabilistic error-state filter rather than the lighter Mahony feedback filter.

## Install and upload

1. Complete the common [`sensors/` Arduino setup](../../README.md#arduino-installation).
2. Verify the device first with [Basic IMU](../../imu_basic/atomS3R_imu_m5_basic/README.md) if this is a new setup.
3. Open `atomS3R_compass_qmekf.ino`, select the AtomS3R-compatible board/port, compile, and upload.
4. Read heading/attitude from the display and NMEA/debug output from USB serial.

The compile-time switches at the top of the sketch control graphics mode, NMEA mode, and the NMEA talker ID.

## Compare

For the lighter feedback-filter implementation using the same application shell, see [AtomS3R Mahony compass](../atomS3R_compass_mahony/README.md).

[Back to Compass / AHRS](../README.md) · [All sensor examples](../../README.md)
