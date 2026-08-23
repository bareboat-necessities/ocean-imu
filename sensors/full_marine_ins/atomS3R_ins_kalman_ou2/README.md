# AtomS3R marine INS — Kalman OU-II

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [full_marine_ins](../README.md) / **atomS3R_ins_kalman_ou2**

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="240" alt="M5Stack AtomS3R">
</p>

This sketch runs the adaptive **SeaStateFusion OU-II** marine estimator on the M5Stack AtomS3R.

Open the sketch: [`atomS3R_ins_kalman_ou2.ino`](atomS3R_ins_kalman_ou2.ino).

## What it does

- reads and calibrates the AtomS3R IMU through the shared hardware layer;
- runs the sea-state fusion loop at 200 Hz;
- learns tilt during startup rather than requiring a perfectly stationary boot;
- performs the filter’s one-shot magnetic north lock;
- uses filter yaw as the main heading output after north lock;
- estimates vessel motion using the adaptive OU-II process and integral-drift regularization;
- provides graphical UI plus serial/NMEA output;
- includes the shared IMU calibration wizard by default.

OU-II uses the repository’s more direct integral drift correction and is useful when you want a full Kalman marine estimator with relatively responsive adaptation to changes in sea state.

## Install and upload

1. Complete the common [`sensors/` Arduino setup](../../README.md#arduino-installation).
2. Verify [Basic IMU](../../imu_basic/atomS3R_imu_m5_basic/README.md) and, ideally, a [Compass / AHRS](../../compass_ahrs/README.md) example first.
3. Open `atomS3R_ins_kalman_ou2.ino` in Arduino IDE.
4. Select the AtomS3R-compatible board and USB port, then compile and upload.
5. Keep the device away from strong local magnetic disturbances while evaluating north lock and heading.

The top-of-file switches control the calibration wizard, graphics UI, NMEA output, and talker ID.

## Compare

- [OU-III marine INS](../atomS3R_ins_kalman_ou3/README.md) — higher-order integral regularization and the main 3-D navigation implementation
- [PII observer](../atomS3R_ins_pii_observer/README.md) — lighter observer-based alternative

[Back to Full marine INS](../README.md) · [All sensor examples](../../README.md)
