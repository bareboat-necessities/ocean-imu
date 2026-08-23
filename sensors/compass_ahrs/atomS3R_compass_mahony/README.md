# AtomS3R Mahony compass

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [compass_ahrs](../README.md) / **atomS3R_compass_mahony**

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="240" alt="M5Stack AtomS3R">
</p>

This sketch implements a tilt-compensated AtomS3R compass using the project’s **Mahony AHRS** backend and shared AtomS3R compass application framework.

Open the sketch: [`atomS3R_compass_mahony.ino`](atomS3R_compass_mahony.ino).

## What it provides

- fused roll, pitch, and yaw from gyro/accelerometer data;
- calibrated magnetometer updates when healthy magnetic data are available;
- graphical compass UI by default;
- NMEA 0183 serial output by default (`HDM`, `XDR`, and `ROT` through the shared app);
- the common AtomS3R IMU calibration workflow.

The Mahony backend is a good embedded baseline when you want a small and responsive attitude estimator without the covariance machinery of the qMEKF.

## Install and upload

1. Complete the common [`sensors/` Arduino setup](../../README.md#arduino-installation).
2. If the board has not been checked yet, run the [basic IMU example](../../imu_basic/atomS3R_imu_m5_basic/README.md) first.
3. Open `atomS3R_compass_mahony.ino`, select the AtomS3R-compatible board/port, compile, and upload.
4. Use the on-device UI for heading/attitude and the USB serial connection for NMEA or debugging.

The compile-time switches near the top of the sketch let you change the default graphics mode, serial NMEA mode, and NMEA talker ID.

## Compare

For the same hardware/application shell with a quaternion multiplicative EKF backend, see [AtomS3R qMEKF compass](../atomS3R_compass_qmekf/README.md).

[Back to Compass / AHRS](../README.md) · [All sensor examples](../../README.md)
