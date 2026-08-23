# Compass and AHRS examples

[ocean-imu](../../README.md) / [sensors](../README.md) / **compass_ahrs**

<p align="center">
  <img src="../../img/devices/AtomS3R_device.svg" width="300" alt="M5Stack AtomS3R">
</p>

These examples turn the AtomS3R into a **tilt-compensated electronic compass / AHRS**. Both applications share the same AtomS3R calibration, sensor conditioning, display, and compass application framework; the attitude estimator is the part that changes.

## Examples

| Example | Estimator | Good choice when |
| --- | --- | --- |
| [atomS3R_compass_mahony](atomS3R_compass_mahony/README.md) | Mahony AHRS | You want a compact, computationally light attitude solution |
| [atomS3R_compass_qmekf](atomS3R_compass_qmekf/README.md) | Quaternion MEKF | You want covariance-based attitude estimation and explicit sensor-noise modeling |

Both sketches default to the graphical compass UI and can emit NMEA 0183 compass/attitude data over serial. The source currently enables `HDM`, `XDR`, and `ROT`-style output through the shared compass application.

## Installation

Follow the common [sensor installation instructions](../README.md#arduino-installation), then open the README for the estimator you want. If this is the first sketch you are running on the board, verify the hardware first with [Basic IMU](../imu_basic/README.md).

## Calibration and magnetic environment

Heading quality depends on magnetometer calibration and installation environment. Keep the device away from temporary magnetic disturbances during evaluation, and use the shared calibration support before treating the heading as a meaningful vessel reference.

## Next step

For heave, 3-D displacement, sea-state adaptation, and wave-motion outputs, continue to [Full marine INS](../full_marine_ins/README.md).

[Back to all sensor examples](../README.md)
