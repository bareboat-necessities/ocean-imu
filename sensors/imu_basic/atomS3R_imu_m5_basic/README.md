# AtomS3R basic IMU

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [imu_basic](../README.md) / **atomS3R_imu_m5_basic**

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="240" alt="M5Stack AtomS3R">
</p>

This is the recommended **first sketch to run on an AtomS3R**. It exercises the M5Unified IMU path and the shared `ocean-imu` AtomS3R helpers without adding an attitude or navigation estimator.

## What it does

The sketch:

- initializes M5Unified and the internal IMU;
- runs the IMU loop at the project default IMU rate;
- maps accelerometer, gyro, and magnetometer readings into the project frame;
- loads the saved calibration blob from NVS when one is present;
- applies saved magnetometer calibration for a calibrated diagnostic channel;
- prints timing, update status, and sensor values over serial.

Open the sketch here: [`atomS3R_imu_m5_basic.ino`](atomS3R_imu_m5_basic.ino).

## Install and upload

1. Complete the common [`sensors/` Arduino setup](../../README.md#arduino-installation).
2. Open `atomS3R_imu_m5_basic.ino` in Arduino IDE.
3. Select the AtomS3R-compatible ESP32-S3 board profile and the USB serial port.
4. Compile and upload.
5. Open Serial Monitor at **115200 baud**.

A successful run prints an initialization message and then periodic IMU diagnostic rows. If no saved calibration exists, the sketch says so explicitly and leaves the magnetometer correction at identity/zero for the diagnostic output.

## What to verify

Before moving to a filter, check that acceleration reacts on the expected axis when the board is tilted, gyro rates react to rotations, and magnetometer values change smoothly when the board is turned away from nearby ferrous objects.

## Continue

- [Mahony compass](../../compass_ahrs/atomS3R_compass_mahony/README.md)
- [qMEKF compass](../../compass_ahrs/atomS3R_compass_qmekf/README.md)
- [All sensor examples](../../README.md)
