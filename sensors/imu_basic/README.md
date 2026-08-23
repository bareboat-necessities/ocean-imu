# Basic IMU examples

[ocean-imu](../../README.md) / [sensors](../README.md) / **imu_basic**

<p align="center">
  <img src="../../img/devices/AtomS3R_device.svg" width="300" alt="M5Stack AtomS3R">
</p>

Use this directory for **first hardware bring-up**. These sketches stay close to the M5Unified sensor API and the AtomS3R calibration/mapping helpers so you can verify the device before adding an AHRS or marine motion filter.

## Examples

- [**atomS3R_imu_m5_basic**](atomS3R_imu_m5_basic/README.md) — reads the AtomS3R IMU through M5Unified, maps the measurements into the project frame, loads saved calibration from NVS when available, and prints accelerometer, gyro, magnetometer, and calibrated magnetometer diagnostics.

## Before running

Follow the common [Arduino installation instructions](../README.md#arduino-installation). For a new setup, confirm that the sketch uploads successfully and that serial output appears before moving on to the compass or full INS examples.

## Next steps

Once raw/mapped sensing looks reasonable:

- continue to [Compass / AHRS](../compass_ahrs/README.md) for roll, pitch, and heading;
- continue to [Full marine INS](../full_marine_ins/README.md) for vessel motion and wave estimation.

[Back to all sensor examples](../README.md)
