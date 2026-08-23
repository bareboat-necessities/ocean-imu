# Sensor and M5Stack examples

[ocean-imu](../README.md) / **sensors**

<p align="center">
  <img src="../img/devices/AtomS3R_device.svg" width="360" alt="M5Stack AtomS3R">
</p>

This directory contains the hardware-facing Arduino examples for `ocean-imu`. The current examples target the **M5Stack AtomS3R** and build upward from direct IMU access to tilt-compensated compass/AHRS applications and complete marine motion estimators.

The examples intentionally form a progression. Start with the basic IMU sketch to verify the board, coordinate mapping, serial connection, and saved calibration. Move to the compass/AHRS examples once raw sensing is healthy, then use a full marine INS example when you want attitude, heading, and wave-motion estimation together.

## Choose an example

| Area | What it demonstrates | Start here |
| --- | --- | --- |
| [IMU basics](imu_basic/README.md) | M5Unified IMU access, body-to-NED mapping, saved calibration, serial diagnostics | First hardware bring-up |
| [Compass / AHRS](compass_ahrs/README.md) | Tilt-compensated heading with Mahony or quaternion MEKF, display UI, NMEA output | Heading and attitude |
| [Full marine INS](full_marine_ins/README.md) | Adaptive marine motion estimation with OU-II, OU-III, or PII observer | Vessel motion / wave estimation |

## Hardware

The sensor sketches are written for the M5Stack AtomS3R and use `M5Unified` for board and IMU access. Connect the AtomS3R over USB-C for power, programming, and serial output.

The reusable filtering code under [`../src/`](../src/) is separate from these hardware examples; the `sensors/` tree shows how that library code is integrated on a real embedded target.

## Arduino installation

1. Install Arduino IDE 2.x or another Arduino-compatible build environment.
2. Install the current Espressif **ESP32 Arduino core** and select the board profile you normally use for the AtomS3R.
3. Install `ocean-imu` as an Arduino library. You can clone/copy this repository into your Arduino libraries directory, or download the repository as a ZIP and use **Sketch → Include Library → Add .ZIP Library…**.
4. Install the library dependencies. [`library.properties`](../library.properties) declares `M5Unified` and `Eigen`; allow the Arduino IDE to install them when prompted, or install them manually with Library Manager.
5. Connect the AtomS3R, select its serial port, open one of the `.ino` sketches linked below, and upload it.

If you are bringing up a new board or a new installation, use [`imu_basic/atomS3R_imu_m5_basic`](imu_basic/atomS3R_imu_m5_basic/README.md) first. It is deliberately much simpler than the full estimators and makes wiring, board-package, sensor, mapping, and calibration problems easier to isolate.

## Suggested bring-up order

1. [Basic AtomS3R IMU](imu_basic/atomS3R_imu_m5_basic/README.md)
2. [Mahony compass](compass_ahrs/atomS3R_compass_mahony/README.md) or [qMEKF compass](compass_ahrs/atomS3R_compass_qmekf/README.md)
3. [OU-II marine INS](full_marine_ins/atomS3R_ins_kalman_ou2/README.md), [OU-III marine INS](full_marine_ins/atomS3R_ins_kalman_ou3/README.md), or [PII observer](full_marine_ins/atomS3R_ins_pii_observer/README.md)

## Coordinate convention

The AtomS3R helper layer converts device readings into the project body/NED-oriented convention used by the filters. Do not bypass that mapping in a full filter unless you also update the corresponding frame assumptions. The basic IMU example is the easiest place to inspect the mapped accelerometer, gyro, and magnetometer values before running an estimator.

## Related repository areas

- [`../src/AtomS3R/`](../src/AtomS3R/) — shared AtomS3R calibration, display, and UI helpers
- [`../src/ahrs/`](../src/ahrs/) — attitude filters
- [`../src/kalman_ou_ii/`](../src/kalman_ou_ii/) — OU-II marine Kalman filter
- [`../src/kalman_ou_iii/`](../src/kalman_ou_iii/) — OU-III marine Kalman filter
- [`../src/pii_observer/`](../src/pii_observer/) — lightweight PII observer
- [Project README](../README.md) — algorithms, publications, tests, and repository-wide build information
