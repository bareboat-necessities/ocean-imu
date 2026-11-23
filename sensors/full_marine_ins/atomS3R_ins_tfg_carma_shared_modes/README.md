# AtomS3R marine INS — shared-mode TFG–CARMA(4,1)

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [full_marine_ins](../README.md) / **atomS3R_ins_tfg_carma_shared_modes**

This sketch runs the project’s primary marine INS estimator, the 27-state shared-mode TFG–CARMA(4,1) filter [`KalmanTfgCarmaSharedModes`](../../../src/kalman_carma_shared_modes/KalmanTfgCarmaSharedModes.h), on the M5Stack AtomS3R. It is [atomS3R_ins_tfg_carma3](../atomS3R_ins_tfg_carma3/README.md) with the three independent CARMA(4,1) motion priors replaced by **shared directional stochastic modes**. Each axis carries one complexified CARMA(4,1) mode in its six-state minimal realization, the axes share its dynamics, and a Hermitian excitation is estimated from the filter's own world-frame acceleration. A **coupling confidence** withdraws the coupling and falls back toward the independent-axis prior when it cannot be trusted. The error state lives on the right-invariant two-frame Lie group of [`src/lie/TwoFrameGroup.h`](../../../src/lie/TwoFrameGroup.h) ([design](../../../docs/tfg-carma-shared-modes.md)).

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="140" alt="M5Stack AtomS3R">
</p>

Open the sketch: [`atomS3R_ins_tfg_carma_shared_modes.ino`](atomS3R_ins_tfg_carma_shared_modes.ino).

> The paired simulation measured 5.9% lower heave RMS and 2.8% lower horizontal RMS than independent-axis TFG–CARMA. The study found poorly calibrated covariance. The firmware has been compiled for the AtomS3R. Sea-trial accuracy, on-device update time and measured stack high-water remain unverified. [Study and limitations](../../../reports/results/tfg_shared_modes/README.md)

## Differences from atomS3R_ins_tfg_carma3

- The filter type is `SeaStateFusion_TfgCarmaSharedModes`, the CARMA3 wrapper with `KalmanTfgCarmaSharedModes<float>` as its inner filter (27 error states instead of 21). Startup, magnetic north lock, tuner, wave-excitation evidence, outputs, sentences, display and calibration paths are the TFG-CARMA3 sketch's.
- `begin(cfg)` selects the deployable configuration: the coupling is estimated from the filter's own acceleration, with the coupling confidence and the default excitation cap. No oracle, RAO or truth information is used.
- The sketch uses the saved magnetometer calibration by default. Build with `-DOCEAN_IMU_SCREEN_MAG_INSTALLATION=1` to enable gyro-aided residual calibration. When enabled, the learner checks for distortion every five seconds. Three consecutive detections allow the live filter to apply the correction gradually once magnetic refinement and excitation checks pass. [Screening and replay results](../../../doc/marine_mag/README.md#production-replay-results)
- The display shows `Fusion: TFG-SHARED`.

Its quaternion counterpart is [atomS3R_ins_carma_shared_modes](../atomS3R_ins_carma_shared_modes/README.md), which runs the same shared modes on the CARMA3 quaternion MEKF error.

## Mounting away from the centre of gravity

Set `IMU_LEVER_ARM_X_M`, `IMU_LEVER_ARM_Y_M` and `IMU_LEVER_ARM_Z_M` to the sensor position relative to CoG in metres. Positive axes point forward, starboard and down. An IMU one metre above CoG uses `(0, 0, -1)`. The wrapper uses the measured gyro to correct specific force before startup and fusion. [Results for six mounting positions](../../../reports/results/tfg_shared_modes/degradation/)

## Memory and timing

The current 27-state filter occupies 16,088 bytes in its single heap allocation on ESP32-S3 (Arduino-ESP32 3.3.11, Xtensa GCC 14.2). Its fixed-size table holds the transition and process-noise coefficients for the tuned wave shape. Each sample uses its measured interval. The default compass display transfers at most 16 rows per loop and defers graphics when less than 1.2 ms remains in the 5 ms polling interval.

Build with `-DOCEAN_IMU_PROFILE_LOOP=1` to emit a `[TIMING]` line every 1,000 accepted samples on the actual device.

| Field | Measurement |
|---|---|
| `fusion_cycles_mean`, `fusion_cycles_max` | CPU cycle counter around the fusion update and any magnetic update due on that sample |
| `loop_work_max_us`, `loops_over_5ms` | Maximum active loop time before pacing, and number of loops that exceeded the 200 Hz budget |
| `sample_gap_min_us`, `sample_gap_max_us` | Shortest and longest interval between accepted IMU poll timestamps |
| `jitter_mean_abs_us`, `gaps_over_7_5ms` | Mean absolute departure from 5,000 µs and count of intervals above 7.5 ms |

Use steady `live=1` windows to assess the operating workload. Keep the clock, graphics and serial settings the same when comparing firmware. At 240 MHz, the whole 5 ms interval is 1.2 million cycles.

### Earlier static size and stack measurements

| | TFG-CARMA3 | shared-modes TFG |
|---|---|---|
| filter object (heap, one allocation at `begin`) | 5.9 kB | 14.5 kB |
| wrapper object (heap) | 3.2 kB | 3.2 kB |
| largest single stack frame | 1.3 kB | 1.7 kB |
| worst static call chain of one update | 3.3 kB | 3.7 kB |

These figures come from `tests/kalman_carma_shared_modes/arduino/csm_size_probe/build.sh`, a static analysis recorded in [`esp32s3_size_stack.csv`](../../../reports/results/tfg_shared_modes/esp32s3_size_stack.csv). They were taken before the excitation-cap change, which adds a few scalars to the filter object. The call-chain figure excludes indirect calls and library frames. The complete sketch reserves 12 KiB for the Arduino loop task, which also runs the startup observer and, when enabled, screened magnetic calibration. The ESP32-S3 compiler reports a 1,872-byte frame for the gyro–mag calibration update with Arduino-ESP32 3.3.11. On-device peak use remains to be measured. The static stack figures predate the current kernel optimizations. Physical board timing and stack high-water must be measured with the running firmware.

## Install and upload

1. Complete the common [`sensors/` Arduino setup](../../README.md#arduino-installation).
2. Verify [Basic IMU](../../imu_basic/atomS3R_imu_m5_basic/README.md) and a [Compass / AHRS](../../compass_ahrs/README.md) example first when commissioning new hardware.
3. Open `atomS3R_ins_tfg_carma_shared_modes.ino`.
4. Select the AtomS3R-compatible board profile and USB port, compile, and upload.
5. Evaluate startup and magnetic north lock in the intended mounting orientation and away from strong transient magnetic disturbances.

Releases published after this sketch was added also carry a prebuilt ZIP; see the [AtomS3R command-line flashing guide](../../../docs/atoms3r-release-flashing.md) and select `atomS3R_ins_tfg_carma_shared_modes`.

The compile-time controls at the top of the sketch include calibration-wizard enable, graphics mode, NMEA mode, and NMEA talker ID.

## Compare

- [TFG-CARMA3 marine INS](../atomS3R_ins_tfg_carma3/README.md) — the same geometry and wrapper with independent-axis priors
- [Shared-modes QMEKF marine INS](../atomS3R_ins_carma_shared_modes/README.md) — the same shared modes with a quaternion MEKF error
- [OU-III marine INS](../atomS3R_ins_kalman_ou3/README.md) — OU-family alternative with integral drift regularization

[Back to Full marine INS](../README.md) · [All sensor examples](../../README.md)
