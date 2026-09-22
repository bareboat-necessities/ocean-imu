# AtomS3R marine INS — Two-Frame Lie-Group (TFG)

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [full_marine_ins](../README.md) / **atomS3R_ins_tfg**

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="240" alt="M5Stack AtomS3R">
</p>

This sketch runs the **right-invariant two-frame Lie-group** estimator, `SeaStateFusionFilter_TFG`, on the M5Stack AtomS3R. It carries the same 21-dimensional state and the same OU wave model as OU-III, but attitude, the world-frame kinematic states and the body-frame biases live on one group, so an attitude correction rotates velocity, position, integral displacement and wave acceleration coherently instead of leaving them behind.

Open the sketch: [`atomS3R_ins_tfg.ino`](atomS3R_ins_tfg.ino).

## What it does

- reads and calibrates the AtomS3R IMU using the shared device layer;
- runs the estimator at 200 Hz;
- learns tilt during startup through the filter's private Mahony proxy, so the boat may already be moving at power-on;
- performs a one-shot magnetic north lock inside the fusion filter, with continuous hard-iron estimation;
- uses filter yaw as the primary heading after north lock while retaining tilt-compensated magnetic heading for diagnostics;
- estimates 3-D marine motion on the two-frame group, with the OU integral chain and its spectral-MSE regularizer;
- adapts the sea-state-dependent estimator parameters from the canonical period statistics;
- arms the accelerometer vibration guard by default, so machinery noise does not rectify into a standing tilt error;
- provides on-device graphics and serial/NMEA output;
- enables the shared IMU calibration wizard by default.

## How it is assembled

`SeaStateFusionFilter_TFG` carries the core of the OU-III orchestrator, but **not** OU-III's wave-direction estimator, its displacement detrender or its alternative frequency trackers: those live in the `SeaStateFusion_OU_III` wrapper class, and TFG has no wrapper. This sketch therefore composes them itself, from the same shared components the OU wrapper uses and driven the same way:

| Output | How this sketch produces it |
| --- | --- |
| heading, roll, pitch | `fusion_.mekf().quaternion()`, the body→world rotation, world frame NED |
| heave, heave rate | `mekf().get_position()` / `get_velocity()`, Z flipped to "up" |
| detrended heave | [`AdaptiveWaveDetrender3D`](../../../src/detrend/AdaptiveWaveDetrender3D.h) on that displacement, given the filter's own period as external frequency guidance |
| wave axis and sign | `heading_frame_acceleration` → [`VesselRaoEqualizer`](../../../src/wave_dir/VesselRaoEqualizer.h) → [`KalmanWaveDirection`](../../../src/wave_dir/KalmanWaveDirection.h) → [`WaveDirectionDetector`](../../../src/wave_dir/WaveDirectionDetector.h) |
| wave envelope | `C_HS · σ_aw · τ² / 2`, the same scalar reduction OU-III publishes as `getDisplacementScale()` |

Every value the application publishes therefore means what it means in the OU-III sketch; the estimator underneath it is the only thing that changed.

## Differences from the OU sketches

- **No IMU lever-arm knob.** `Kalman3D_Wave_TFG` has no off-centre-of-gravity correction, so there is no `IMU_LEVER_ARM_*_M` block. Mount the device as close to the vessel centre of gravity as the installation allows, or use [OU-III](../atomS3R_ins_kalman_ou3/README.md) where the offset is large enough to matter.
- **Gyroscope noise is a density, not a per-axis sigma.** `Config::gyro_noise_density` replaces OU-III's `sigma_g` vector, and the sample-rate scaling the accelerometer and magnetometer get must not be applied to it a second time.
- **One frequency tracker.** TFG reports a wave period from the canonical estimator; there is no `TrackerType` choice and no `ZERO_CROSSINGS_*` tuning block.
- **It needs a bigger loop-task stack.** `Kalman3D_Wave_TFG` runs a 21-state covariance through `Phi P Phi^T` and a Joseph update, and Eigen builds each triple product through full 21×21 temporaries. The named locals are member scratch buffers, which cuts one live fusion step from ~27 kB of stack to ~19 kB without changing a single output bit, but the expression temporaries stay: removing them would change Eigen's GEMM path and with it the filter's arithmetic. Painted-stack high-water mark over a 400 s host run of the whole sketch-level step:

  | | `-O2` | `-Os` |
  | --- | --- | --- |
  | before | 27.2 kB | 27.2 kB |
  | after | 18.9 kB | 15.3 kB |

  The Arduino-ESP32 loop task gets 8 kB by default, so the sketch still raises it, with `SET_LOOP_TASK_STACK_SIZE(32 * 1024)`. Without that the filter overflows the stack the moment it leaves `Cold` and starts running the MEKF — tens of seconds after boot, not at startup, so it would not look like a stack problem. Lower the value only against a fresh measurement on the device.
- **Vibration guard is exposed instead of the OU clamp knobs.** `ACC_VIBRATION_GUARD_HZ` near the top of the sketch sets the front-end corner; zero removes the guard and restores the unconditioned measurement path.

## Install and upload

1. Complete the common [`sensors/` Arduino setup](../../README.md#arduino-installation).
2. Verify [Basic IMU](../../imu_basic/atomS3R_imu_m5_basic/README.md) and a [Compass / AHRS](../../compass_ahrs/README.md) example first when commissioning new hardware.
3. Open `atomS3R_ins_tfg.ino`.
4. Select the AtomS3R-compatible board profile and USB port, compile, and upload.
5. Evaluate startup and magnetic north lock in the intended mounting orientation and away from strong transient magnetic disturbances.

The compile-time controls at the top of the sketch include calibration-wizard enable, graphics mode, NMEA mode, and NMEA talker ID.

## More about TFG

The convention contract — rotation direction, error definition and correction side, stated together — is in [`docs/tfg-design.md`](../../../docs/tfg-design.md), and it is pinned by `tests/kalman_tfg/convention-test.cpp`. The root [project README](../../../README.md) links the TFG paper and the comparison studies against the OU families. The implementation used by this sketch lives under [`../../../src/kalman_tfg/`](../../../src/kalman_tfg/), on the Lie-group operations in [`../../../src/lie/`](../../../src/lie/).

## Compare

- [OU-III marine INS](../atomS3R_ins_kalman_ou3/README.md) — same state and wave model, quaternion-MEKF geometry, plus the lever-arm correction
- [OU-II marine INS](../atomS3R_ins_kalman_ou2/README.md) — more direct integral drift correction
- [PII observer](../atomS3R_ins_pii_observer/README.md) — much lighter observer-based alternative

[Back to Full marine INS](../README.md) · [All sensor examples](../../README.md)
