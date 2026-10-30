# AtomS3R marine INS — TFG-CARMA3 (experimental)

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [full_marine_ins](../README.md) / **atomS3R_ins_tfg_carma3**

This sketch runs the **experimental** TFG-CARMA3 Kalman filter, **SeaStateFusion TFG-CARMA3**, on the M5Stack AtomS3R. It is the [CARMA3 sketch](../atomS3R_ins_carma3/README.md) with the estimator replaced: the same CARMA(4,1) motion priors on X, Y and Z (no `S` states), carried on the right-invariant two-frame Lie group of [`src/lie/TwoFrameGroup.h`](../../../src/lie/TwoFrameGroup.h) instead of CARMA3's quaternion MEKF error. The fusion wrapper, startup, magnetometer handling and adaptation are CARMA3's code, so the two sketches differ only in the error geometry ([`docs/tfg-carma3-research.md`](../../../docs/tfg-carma3-research.md)).

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="140" alt="M5Stack AtomS3R">
</p>

Open the sketch: [`atomS3R_ins_tfg_carma3.ino`](atomS3R_ins_tfg_carma3.ino).

> **Research sketch, not the shipping estimator.** The TFG-CARMA3 results are simulation-only; there is no hardware, on-target timing or sea-trial evidence and no stability analysis. Translation matches CARMA3 (including CARMA3's still-IMU and mixed-sea weaknesses); attitude and bias are somewhat better where the CARMA prior misfits. Use [OU-III](../atomS3R_ins_kalman_ou3/README.md) for deployment.

## What it does

- reads and calibrates the AtomS3R IMU using the shared device layer;
- stops at startup, before the estimator, when the saved calibration cannot be applied in full (none saved, accelerometer fitted against another gravity, magnetometer fitted on another data source, or a missing set), reports the reason on serial every 5 s and offers the wizard on tap; it never runs on raw sensor values;
- runs the estimator at 200 Hz, with the BMI270 accelerometer switched from its 100 Hz reset rate to the gyroscope's 200 Hz so every gyro sample is used;
- rejects sporadic single-sample IMU spikes before calibration and fusion with the shared [`ImuSpikeFilter`](../../../src/AtomS3R/AtomS3R_ImuSpikeFilter.h) (disable with `-DATOMS3R_IMU_SPIKE_FILTER=0`);
- learns tilt during startup while the vessel may already be moving;
- performs a one-shot magnetic north lock inside the sea-state fusion filter;
- uses filter yaw as the primary heading after north lock while retaining tilt-compensated magnetic heading for diagnostics;
- keeps north on the startup/refined magnetic reference: continuous hard-iron learning is disabled;
- estimates 3-D marine motion with CARMA(4,1) priors whose frequency and amplitude follow the operating point measured by the shared OU-III tuner;
- optionally corrects an IMU mounted away from the centre of gravity (`IMU_LEVER_ARM_*_M` at the top of the sketch, zero by default);
- provides on-device graphics and serial/NMEA output (the same sentences as the OU-III sketch);
- enables the shared IMU calibration wizard by default.

Differences from the OU-III sketch are CARMA3's:

- the IMU lever-arm correction is removed at the fusion wrapper's input, so the startup tilt, the magnetic reference, the tuner and the wave direction all see centre-of-gravity motion; OU-III removes it in its MEKF's accelerometer update only. The [off-CoG study](../../../docs/carma-imu-lever-arm.md) compares the two;
- no `S = 0` pseudo-measurement, so there is no pseudo-update cadence setting.

Unlike the TFG sketch, it runs on the default Arduino loop-task stack: the filter keeps its covariance scratch in members and needs no enlarged stack.

## Install and upload

Install a prebuilt release without compiling using the
[AtomS3R command-line flashing guide](../../../docs/atoms3r-release-flashing.md).
Select `atomS3R_ins_tfg_carma3`. For an Arduino IDE source build, follow the steps below.

1. Complete the common [`sensors/` Arduino setup](../../README.md#arduino-installation).
2. Verify [Basic IMU](../../imu_basic/atomS3R_imu_m5_basic/README.md) and a [Compass / AHRS](../../compass_ahrs/README.md) example first when commissioning new hardware.
3. Open `atomS3R_ins_tfg_carma3.ino`.
4. Select the AtomS3R-compatible board profile and USB port, compile, and upload.
5. Evaluate startup and magnetic north lock in the intended mounting orientation and away from strong transient magnetic disturbances.

The compile-time controls at the top of the sketch include calibration-wizard enable, graphics mode, NMEA mode, and NMEA talker ID.

## More about TFG-CARMA3

[`docs/tfg-carma3-research.md`](../../../docs/tfg-carma3-research.md) describes the filter, its verification, the paired benchmark against CARMA3 and OU-III, the ESP32-S3 build-size and stack figures, and the recommendation. The experimental article is [`doc/kalman_tfg_carma/kalman_tfg_carma.tex`](../../../doc/kalman_tfg_carma/kalman_tfg_carma.tex). The implementation lives under [`../../../src/kalman_carma/`](../../../src/kalman_carma/) (`Kalman3D_Wave_TFG_CARMA3.h`, `SeaStateFusionFilter_TFG_CARMA3.h`).

## Compare

- [CARMA3 marine INS](../atomS3R_ins_carma3/README.md) — the same prior and wrapper with a quaternion MEKF error
- [TFG marine INS](../atomS3R_ins_tfg/README.md) — OU-III's state and wave model on the two-frame Lie group
- [OU-III marine INS](../atomS3R_ins_kalman_ou3/README.md) — the shipping 3-D navigation estimator

[Back to Full marine INS](../README.md) · [All sensor examples](../../README.md)
