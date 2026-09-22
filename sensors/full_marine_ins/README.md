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
| [atomS3R_ins_tfg](atomS3R_ins_tfg/README.md) | Right-invariant two-frame Lie-group INS | Same state and wave model as OU-III on one group, so attitude corrections move the kinematic states with them |
| [atomS3R_ins_pii_observer](atomS3R_ins_pii_observer/README.md) | Adaptive PII observer + Mahony | Much lighter computation; primarily vertical/heave-oriented motion estimation |
| [atomS3R_ins_nlo](atomS3R_ins_nlo/README.md) | Time-varying-gain nonlinear observer | Observer-based 3-D estimator with published gains; no covariance propagation |

The OU-II, OU-III and TFG sketches learn tilt during startup, perform a one-shot magnetic north lock in the sea-state fusion filter, and then use filter yaw as the primary heading output. The PII application uses Mahony attitude with a lightweight adaptive vertical observer. The NLO application bootstraps attitude with a Mahony startup stage, hands over to the nonlinear observer, and takes the magnetometer as a continuous yaw-only magnetic reference.

TFG differs from the OU sketches in how it is assembled rather than in what it publishes. `SeaStateFusionFilter_TFG` carries the core of the OU-III orchestrator but not its wave-direction estimator or its displacement detrender, which live in the OU wrapper class; the TFG sketch composes those from the same shared components, so its outputs mean what they mean elsewhere on this page. It also has no IMU lever-arm correction — see its [README](atomS3R_ins_tfg/README.md) for the full list of differences.

## Installation and bring-up

Follow the common [Arduino installation instructions](../README.md#arduino-installation). Before debugging one of these larger applications, first verify the same board with [Basic IMU](../imu_basic/README.md) and preferably one of the [Compass / AHRS](../compass_ahrs/README.md) examples. That separates sensor/calibration problems from full estimator behavior.

The full sketches enable the shared IMU calibration wizard by default through `SEA_STATE_ENABLE_WIZARD=1`. They also default to graphical UI and NMEA serial output; those behaviors can be changed with the compile-time switches near the top of each sketch.

## Which should I use?

Start with **OU-III** if you want the project’s primary full 3-D marine navigation estimator. Try **OU-II** when you prefer its more direct drift regularization and adaptation response. Try **TFG** when you want the same state and wave model with a geometrically consistent attitude correction, and can mount close to the centre of gravity. Use the **PII observer** when MCU cost and simplicity matter more than the full Kalman-state solution. Use the **NLO** when you want the observer-based alternative to the Kalman filters: explicit injection terms with published, frequency-scheduled gains rather than a propagated covariance.

## Related material

The mathematical descriptions and validation links are collected in the [project README](../../README.md). Source implementations live under [`../../src/kalman_ou_ii/`](../../src/kalman_ou_ii/), [`../../src/kalman_ou_iii/`](../../src/kalman_ou_iii/), [`../../src/kalman_tfg/`](../../src/kalman_tfg/), [`../../src/pii_observer/`](../../src/pii_observer/), and [`../../src/nlo/`](../../src/nlo/). The TFG convention contract is in [`../../docs/tfg-design.md`](../../docs/tfg-design.md).

[Back to all sensor examples](../README.md)
