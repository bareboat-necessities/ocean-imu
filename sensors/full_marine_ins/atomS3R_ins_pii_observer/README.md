# AtomS3R marine motion estimator — PII observer

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [full_marine_ins](../README.md) / **atomS3R_ins_pii_observer**

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="240" alt="M5Stack AtomS3R">
</p>

This sketch combines the project’s **Mahony attitude solution** with the adaptive **vertical PII observer** for a lower-cost marine motion estimator on the M5Stack AtomS3R.

Open the sketch: [`atomS3R_ins_pii_observer.ino`](atomS3R_ins_pii_observer.ino).

## What it provides

- magnetic compass heading derived from Mahony yaw;
- AHRS roll, pitch, and yaw;
- heave / vertical displacement estimate;
- wave-envelope estimate;
- common AtomS3R calibration support and optional calibration wizard;
- on-device UI and serial/NMEA output.

The sketch explicitly maps the calibrated device/body NED convention into the Z-up convention used by its Mahony implementation. It also seeds Mahony once from accelerometer plus magnetometer data at startup so heading does not have to converge slowly from an identity quaternion.

## Install and upload

1. Complete the common [`sensors/` Arduino setup](../../README.md#arduino-installation).
2. Verify [Basic IMU](../../imu_basic/atomS3R_imu_m5_basic/README.md) first on new hardware.
3. Open `atomS3R_ins_pii_observer.ino`, select the AtomS3R-compatible board and USB port, then compile and upload.
4. Check attitude and magnetic heading before relying on the vertical motion output.

The calibration wizard is enabled by default. Graphics/NMEA behavior and true-versus-magnetic heading options are controlled by compile-time switches near the top of the sketch.

## When to use it

Choose this example when you want a computationally light embedded marine motion solution and primarily care about attitude plus vertical/heave behavior. For the repository’s full covariance-based 3-D estimators, use [OU-II](../atomS3R_ins_kalman_ou2/README.md) or [OU-III](../atomS3R_ins_kalman_ou3/README.md).

[Back to Full marine INS](../README.md) · [All sensor examples](../../README.md)
