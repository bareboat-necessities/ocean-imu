# AtomS3R magnetometer diagnostics

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [imu_basic](../README.md) / **atomS3R_mag_diagnostics**

Use this sketch when a compass reads wrong headings and you need to know whether the cause is the magnetometer wiring/axes, the saved calibration, the calibration location, or the alignment between the magnetometer and the accelerometer. It never changes the saved calibration (see [Magnetometer setting and drift](#magnetometer-setting-and-drift) for the sensor setting).

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="140" alt="M5Stack AtomS3R">
</p>

Open the sketch: [`atomS3R_mag_diagnostics.ino`](atomS3R_mag_diagnostics.ino).

## Steps

Open Serial Monitor at **115200 baud**. Tap the screen, or send `n`, to advance.

1. **Boot report** (`[DIAG]`): whether the BMM150 answers on the main I2C bus or only behind the BMI270 AUX bus, the BMI270 AUX polling rate, the magnetometer setting, and the saved calibration.
2. **STILL**, 10 s (`[STILL]`): lay the board flat, screen up, away from metal and do not touch it. Reports the magnetometer update interval, raw and calibrated noise, the heading and dip scatter of single readings, the IMU temperature and update mask counts, and the dip angle with the saved calibration.
3. **ROTATE** (`[ROT]`): turn and tilt slowly (under about 45°/s) through every direction, including upside down and on each edge. This uses the calibration wizard's capture and fit code, and also records slow-motion poses and gyro-measured rotations. It ends by itself when the capture is complete (or tap to stop, or after 3 minutes).
4. **Analysis** (`[FIT]`, `[FIELD]`, `[AXIS]`, `[ALIGN]`, `[VERDICT]`):
   - a fresh calibration fitted here, compared with the saved one (and with the previous run's fit, so a second run at the same place measures calibration repeatability);
   - field-strength spread and dip angle across those poses for both calibrations (a calibration that fits keeps both nearly constant);
   - an axis test that checks every signed axis permutation of the magnetometer against the rotation the gyro measured;
   - the rotation between the magnetometer and accelerometer frames that an ellipsoid fit cannot see, and the heading error it causes.
5. **LIVE** (`[LIVE]`): headings from the saved, fresh, and fresh-plus-alignment models. Hold the board level, then tilt it about ±20° at each cardinal direction; a correct model keeps the heading steady. Tap to run again.

The verdicts compare the dip with `DIAG_EXPECTED_DIP_DEG` (default 66.5°, Fair Lawn NJ). Look up your local inclination (NOAA magnetic field calculator) and build with `-DDIAG_EXPECTED_DIP_DEG=<deg>`; use a negative value in the southern hemisphere.

## Magnetometer setting and drift

Every AtomS3R sketch, the calibration wizard and this one start the BMM150 in a low-noise setting (47 XY / 41 Z repetitions at M5Unified's 30 Hz output rate, instead of the single repetition M5Unified leaves after reset) and apply Bosch's trim-register compensation with RHALL, the sensor's own temperature compensation. A saved magnetometer calibration made before compensation keeps applying to uncompensated data (compensation off) until the magnetometer is recalibrated. While the sketch waits for a tap (or in LIVE), send over serial:

- `q`: read the BMM150 configuration (`[MAGCFG] current`);
- `o`: switch to the M5Unified driver default (1/1 repetitions) until the next reboot, to compare;
- `p`: switch back to the low-noise setting;
- `c`: toggle Bosch trim/RHALL compensation; the saved calibration applies only to the data source it was fitted on, and the choice made here is kept until the next reboot;
- `d`: drift monitor. Leave the board untouched; every 5 s it prints the average raw field, its change since the start, noise, and IMU temperature (`[DRIFT]`). A change that follows temperature is sensor drift; a stable field means heading changes come from where the board is placed. Tap or send `d` again to stop.

Every `[STILL]`, `[FIT]` and `[VERDICT]` block names the setting in use. Two consecutive runs in the same setting show calibration repeatability in `[FIT] previous-run-vs-fresh`.

## Install and upload

Install a prebuilt release without compiling using the
[AtomS3R command-line flashing guide](../../../docs/atoms3r-release-flashing.md).
Select `atomS3R_mag_diagnostics`. For an Arduino IDE source build, complete the common [`sensors/` Arduino setup](../../README.md#arduino-installation), open `atomS3R_mag_diagnostics.ino`, select the AtomS3R board profile and port, then compile and upload.

Copy the whole serial log, from the `[DIAG]` header through the `[VERDICT]` block and a few `[LIVE]` lines, when reporting a compass problem.

## Continue

- [Basic IMU](../atomS3R_imu_m5_basic/README.md)
- [Mahony compass](../../compass_ahrs/atomS3R_compass_mahony/README.md)
- [qMEKF compass](../../compass_ahrs/atomS3R_compass_qmekf/README.md)

## Example Calibration Report

```
[MAGCFG] BMM150 low-noise: 47/41 repetitions (was 1/1), mode=0x38
[MAGCFG] BMM150 Bosch compensation on: x1=0 y1=0 x2=26 y2=26 xy1=29 xy2=-3 z1=23337 z2=722 z3=0 z4=0 xyz1=6574
[DIAG] magnetometer: low-noise (47/41 reps, 30 Hz), Bosch-compensated
[MAGTIME CAL] valid=1 delay_ms=36.0 allowance_ms=27.0 applied=0

[DIAG] ===== AtomS3R magnetometer diagnostics =====
[DIAG] imu_type=6 (bmi270=6) expected_dip=66.5 deg
[DIAG] BMI270 addr=0x68 AUX_CONF=0x46 (aux poll 25.00 Hz) AUX_DEV_ID=0x20 (i2c 0x10) AUX_IF_CONF=0x4F (data mode) AUX_RD_ADDR=0x42 PWR_CTRL=0x0F
[DIAG] BMM150 main-bus instance: NONE (magnetometer is only behind the BMI270 AUX bus)
[DIAG] main-bus probe 0x10: no ack
[DIAG] main-bus probe 0x11: no ack
[DIAG] main-bus probe 0x12: no ack
[DIAG] main-bus probe 0x13: no ack
[DIAG] saved calibration:
  build_mode: m5_imu_api
  ok: A=1 G=1 M=1
  accel: g=9.8025608 (matches g_cal_local) T0=36.60 rms_mag=0.0034
    b0=[-0.10482 -0.14907 0.06570]
    k =[0.000000 0.000000 0.000000] clamp T=[-1000.0 1000.0]
    S (a_cal = S*(a_raw - bias(T))):
      [0.999178350, 0.001822500, 0.000486332]
      [0.001822500, 0.998607159, -0.003130334]
      [0.000486332, -0.003130334, 0.997308433]
      diag=[0.999178 0.998607 0.997308], offdiag_rms=0.002110
    thermal=UNLEARNED/SPAN_SMALL meta=bound
    holds=13 blocks=284 capture=87s T_seen=[36.31 37.12] k_range=[0.00 0.00]
    bias_sd=[0.0019 0.0019 0.0018] cross_sd=[0.00025 0.00034 0.00028] k_sd=[-1.00000 -1.00000 -1.00000]
    sigma_obs=0.0031 heldout_cv rms/max=0.0068/0.0152 verify rms/max=0.0054/0.0084
  gyro:  T0=36.83
    b0=[0.002320 -0.010273 -0.007242]
    k =[0.000000 0.000000 0.000000]
    S (w_cal = S*(w_raw - bias(T)) ; S=I (stationary bias-only fit)):
      [1.000, 0.000, 0.000]
      [0.000, 1.000, 0.000]
      [0.000, 0.000, 1.000]
      diag=[1.000000 1.000000 1.000000], offdiag_rms=0.000000
    thermal=BIAS_ONLY reason=0 bins=7 evidence T=[36.71 36.93] clamp T=[36.83 36.83]
    k_sigma=[-1.0000000 -1.0000000 -1.0000000] information=0.020 meta=bound
  mag: field_uT=47.715 rms=0.3593
    b=[26.087 10.853 -17.710]
    A (m_cal = A*(m_raw - b)):
      [0.938192546, 0.002844804, 0.000681775]
      [0.002844804, 0.945208788, -0.002084012]
      [0.000681775, -0.002084012, 0.979314506]
      diag=[0.938193 0.945209 0.979315], offdiag_rms=0.002074
[DIAG] saved mag calibration applied at runtime: YES
[DIAG] next: STILL. Lay the board flat, screen up, away from metal; tap and do not touch.
[DIAG] (serial: 'q' read mag setting, 'o' driver default, 'p' low-noise, 'd' drift monitor)
[MAGCFG] current: chip=0x32 power=0x01 mode=0x38 rep_xy=23 (n=47) rep_z=40 (n=41)
[STILL] magnetometer mode: low-noise (47/41 reps, 30 Hz), Bosch-compensated
[STILL] polls=2000 accel+gyro=990 gyro-only=977 accel-only=0 none=33 mag-bit=246
[STILL] distinct mag readings=247 interval mean=40.4 ms min=30.0 max=84.9 ms (=> 24.7 Hz)
[STILL] raw mag mean [39.63 33.57 25.69] noise sd [1.08 0.80 0.51]
[STILL] cal mag mean [12.80 21.43 42.46] noise sd [1.02 0.76 0.50] |m|=49.26
[STILL] accel mean [-0.060 0.132 -9.789] |a|=9.790 tilt=0.9 deg; gyro mean [-0.067 0.434 0.553] deg/s
[STILL] single-reading scatter with saved cal: heading sd=2.4 deg, dip sd=1.0 deg; temp=39.5 C
[STILL] dip with saved cal=59.1 (expected 66.5), heading=299.7; accel D is negative (screen up)
[DIAG] next: ROTATE. Turn and tilt slowly (under ~45 deg/s) through every direction, including upside down.
```
