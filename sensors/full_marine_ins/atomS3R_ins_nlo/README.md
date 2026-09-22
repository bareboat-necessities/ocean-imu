# AtomS3R marine INS — time-varying-gain NLO

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [full_marine_ins](../README.md) / **atomS3R_ins_nlo**

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="240" alt="M5Stack AtomS3R">
</p>

This sketch runs the **time-varying-gain nonlinear observer** of Bryne, Fossen and Johansen on the M5Stack AtomS3R. It is the hardware application of the repository's observer-based (non-Kalman) full 3-D estimator: attitude and translational states are driven by explicit injection terms with scheduled gains instead of by a propagated covariance.

Open the sketch: [`atomS3R_ins_nlo.ino`](atomS3R_ins_nlo.ino).

## What it does

- reads and calibrates the AtomS3R IMU using the shared device layer;
- runs the observer at 200 Hz through [`TimeVarGainNloAdapter`](../../../src/nlo/TimeVarGainNLO_Adapter.h);
- bootstraps attitude with the adapter's Mahony startup stage, then seeds and hands over to the observer;
- uses the calibrated 3-D magnetometer as a yaw-only magnetic reference, so heading is observable;
- aids all three translational axes with the paper's virtual zero-mean integrated-position measurement, since the device carries no position reference;
- schedules the observer's scalar translational tuning parameter `theta` on the tracked dominant wave frequency, which keeps the aiding loop's high-pass corner below the wave band;
- reports heave, heave speed, roll, pitch and heading;
- provides on-device graphics and serial/NMEA output;
- enables the shared IMU calibration wizard by default.

## Configuration

The observer is instantiated as

```cpp
TimeVarGainNloAdapter<false, NloMagType::Magnetometer, float, TrackerType::PLL>
```

- `WithGNSS = false` — no position reference on this device.
- `Mag = NloMagType::Magnetometer` — the body-frame magnetometer supplies the yaw reference. Without it yaw and z-gyro bias are unobservable.
- `TrackerType::PLL` — the frequency tracker the committed NLO simulation study runs.

All observer gains stay at the values published in the paper's tuning section; the adapter's own defaults are those values. The compile-time switches at the top of the sketch are calibration-wizard enable, graphics mode, NMEA mode, NMEA talker ID, true-versus-magnetic heading, magnetic declination, heading user offset, magnetometer gate strictness, and `SEA_STATE_NLO_FIXED_THETA`.

`SEA_STATE_NLO_FIXED_THETA=1` replaces the wave-frequency schedule with the paper's fixed `theta = 1`. It is off by default: at `theta = 1` the aiding loop's high-pass corner sits at 0.41 rad/s, a 15 s period, which is inside the ocean wave band and attenuates and phase-shifts the wave the device is meant to measure.

## Outputs and their limits

Use the vertical channel as the navigation product. Horizontal position and velocity states exist so that the estimated specific force used as the attitude reference is a genuine estimate; with only a virtual zero-mean constraint and no horizontal position aid they are not a calibrated surge/sway output. This is the observer's own documented caveat — see the header comment of [`../../../src/nlo/TimeVaryingGainNLO.h`](../../../src/nlo/TimeVaryingGainNLO.h).

Reported position and velocity have the unobservable DC component removed by the adapter's report high-pass, so heave is wave-band motion about the mean sea surface. The sketch additionally runs the shared [`AdaptiveWaveDetrender`](../../../src/detrend/AdaptiveWaveDetrender.h) on heave, and it is the detrended value that is sent as the NMEA heave transducer sentence, matching the other sketches in this directory.

Gyro bias is the observer's own state, driven by its integral gain. Unlike the Kalman and PII sketches, this sketch runs no separate stillness-gated bias average, and the rate-of-turn output is corrected with the observer's estimate. There is no accelerometer-bias state, and no lever-arm input: `TimeVaryingGainNLO` takes none, so mount the device as close to the vessel centre of gravity as practical.

## Install and upload

1. Complete the common [`sensors/` Arduino setup](../../README.md#arduino-installation).
2. Verify [Basic IMU](../../imu_basic/atomS3R_imu_m5_basic/README.md) and a [Compass / AHRS](../../compass_ahrs/README.md) example first when commissioning new hardware.
3. Open `atomS3R_ins_nlo.ino`.
4. Select the AtomS3R-compatible board profile and USB port, compile, and upload.
5. Let the startup stage finish before judging the output: the observer's attitude gains ramp down over the first 100 s, and the `theta` schedule needs the wave tracker to lock.

## More about the NLO

The method is described in [Nonlinear Observer with Time-Varying Gains for Inertial Navigation Aided by Satellite Reference Systems in Dynamic Positioning](https://torarnj.folk.ntnu.no/TimeVarGain.pdf), linked from the [project README](../../../README.md). The implementation lives under [`../../../src/nlo/`](../../../src/nlo/), its simulation harness under [`../../../tests/nlo/`](../../../tests/nlo/), and its charts under [`../../../doc/nlo/`](../../../doc/nlo/).

## Compare

- [OU-III marine INS](../atomS3R_ins_kalman_ou3/README.md) — the main 3-D Kalman navigation estimator
- [OU-II marine INS](../atomS3R_ins_kalman_ou2/README.md) — more direct integral drift correction
- [TFG marine INS](../atomS3R_ins_tfg/README.md) — OU-III's state and wave model on a two-frame Lie group
- [PII observer](../atomS3R_ins_pii_observer/README.md) — much lighter observer-based alternative

[Back to Full marine INS](../README.md) · [All sensor examples](../../README.md)
