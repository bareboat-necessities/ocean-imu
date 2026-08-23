# AtomS3R marine INS — Kalman OU-III

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [full_marine_ins](../README.md) / **atomS3R_ins_kalman_ou3**

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="240" alt="M5Stack AtomS3R">
</p>

This sketch runs the adaptive **SeaStateFusion OU-III** estimator on the M5Stack AtomS3R. It is the hardware application corresponding to the repository’s main higher-order 3-D marine navigation filter.

Open the sketch: [`atomS3R_ins_kalman_ou3.ino`](atomS3R_ins_kalman_ou3.ino).

## What it does

- reads and calibrates the AtomS3R IMU using the shared device layer;
- runs the estimator at 200 Hz;
- learns tilt during startup while the vessel may already be moving;
- performs a one-shot magnetic north lock inside the sea-state fusion filter;
- uses filter yaw as the primary heading after north lock while retaining tilt-compensated magnetic heading for diagnostics;
- estimates 3-D marine motion with the OU-III wave/process model and higher-order integral regularization;
- adapts the sea-state-dependent estimator parameters used by the full filter;
- provides on-device graphics and serial/NMEA output;
- enables the shared IMU calibration wizard by default.

## Install and upload

1. Complete the common [`sensors/` Arduino setup](../../README.md#arduino-installation).
2. Verify [Basic IMU](../../imu_basic/atomS3R_imu_m5_basic/README.md) and a [Compass / AHRS](../../compass_ahrs/README.md) example first when commissioning new hardware.
3. Open `atomS3R_ins_kalman_ou3.ino`.
4. Select the AtomS3R-compatible board profile and USB port, compile, and upload.
5. Evaluate startup and magnetic north lock in the intended mounting orientation and away from strong transient magnetic disturbances.

The compile-time controls at the top of the sketch include calibration-wizard enable, graphics mode, NMEA mode, and NMEA talker ID.

## More about OU-III

The root [project README](../../../README.md) links the OU-III paper, startup/stability/adaptation studies, and the simulation validation material. The implementation used by this sketch lives under [`../../../src/kalman_ou_iii/`](../../../src/kalman_ou_iii/).

## Compare

- [OU-II marine INS](../atomS3R_ins_kalman_ou2/README.md) — more direct integral drift correction
- [PII observer](../atomS3R_ins_pii_observer/README.md) — much lighter observer-based alternative

[Back to Full marine INS](../README.md) · [All sensor examples](../../README.md)
