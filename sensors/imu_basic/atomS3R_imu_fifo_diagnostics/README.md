# AtomS3R IMU FIFO diagnostics

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [imu_basic](../README.md) / **atomS3R_imu_fifo_diagnostics**

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="240" alt="M5Stack AtomS3R">
</p>

This sketch measures, on the device, how good the IMU data is. It covers two things:

- **timing:** whether every sample arrives, on time and evenly spaced;
- **accuracy:** whether the calibrated accelerometer, gyro and magnetometer read what physics says they should.

It reads the sensors the way the INS sketches are meant to move to: the BMI270 FIFO, with sample times from the sensor clock and a sensor-triggered magnetometer. It applies the saved calibration as the INS sketches do. It estimates nothing, never changes the saved calibration, and changes no INS sketch.

Open the sketch: [`atomS3R_imu_fifo_diagnostics.ino`](atomS3R_imu_fifo_diagnostics.ino).

## What it changes compared with the INS sketches

| | INS sketches today | This sketch |
|---|---|---|
| Accelerometer / gyro rate | 100 Hz / 200 Hz (driver defaults); only frames with both new are used | both 200 Hz |
| Reading | data registers, polled by the loop; a busy loop loses samples | BMI270 FIFO, drained every loop; a busy loop delays samples but loses none |
| Sample time | host `micros()` when the loop polled | the sample's place on the sensor's own clock (a frame every 128 sensortime ticks), mapped onto `micros()` |
| Magnetometer | BMM150 free-running at 30 Hz, read by the BMI270 at 25 Hz, so the age of each reading drifts | BMM150 in forced mode; the BMI270 starts a conversion every 40 ms and reads it 40 ms later, so every readout is one new conversion at a fixed point on the sensor clock |

Ranges stay at ±8 g and ±2000 °/s, the ranges M5Unified's conversion assumes, so samples are in exactly the units `readImuMapped()` delivers (a host test checks this bit for bit).

## Running it

Calibrate first with the wizard in any INS sketch; without a calibration the sketch reports the raw sensor. Then:

1. Upload. Open a serial capture at 115200 baud, for example `python3 -m serial.tools.miniterm --raw <port> 115200 | tee fifo.log`, and reset the board. Frame lines (`R`) are off, so the summary lines stay readable; send `q` to turn them on for a full capture.
2. **Still poses:** put the board on a firm, level surface away from metal, and leave it untouched for 5 s or more. Then turn it to another resting orientation and leave it again. Use at least six orientations: level, each side down, upside down, nose up and nose down, plus one turned flat on the table. Every pose produces an `S` line, and every turn between poses a `W` line.
3. **Rocking:** hold the board and rock it back and forth, about one rock per second, tilting 20–40° and changing direction, until the TIME page shows `dly OK` (usually 1–3 minutes). Rocking also checks |m| in motion.
4. Leave it running past one or two display bursts (every 10 s), then stop the capture.

Tap the screen to page through **TIME**, **ACC/GYR**, **MAG** and **ROT**.

## Timing criteria (`T` lines, every 10 s)

| Criterion | Pass when |
|---|---|
| rate | 200 Hz on the sensor's own clock, no frame missing. The host-clock rate differs by the BMI270 oscillator's error, reported as `clock_ppm` (one AtomS3R showed about +8000 ppm in its first second: frames every 5.04 ms, 198.4 Hz on the ESP32 clock) |
| no loss | no frame lost across the display bursts (400 ms and 1.5 s, alternating every 10 s; `s` on serial runs another 1.5 s burst) |
| spacing | every delivered sample spacing (after warm-up) between 4.9 and 5.1 ms |
| mag | 25 Hz readouts, each with the BMM150 data-ready bit set (a new conversion), each 40 ms after the previous on the sensor clock |
| delay | three magnetometer-delay fits, each at most 4 ms wide, within 4 ms of each other |

The delay is fitted by the calibration wizard's timing code (`MagTimingCal`) from the rotation the gyro measured. It is expressed in the convention the estimators' magnetic transport uses: a gyro sample is the rate over the interval ending at it. It therefore includes:

- the conversion aperture;
- the time from the conversion to the readout;
- the gyro's own filter delay.

That is the number a filter needs.

## Accuracy criteria (`A` lines, every 10 s)

Each check uses still poses (the board resting, detected automatically) or the rotation between two of them. The limits are diagnostic, not sensor specifications, and each can be changed with a `FIFO_DIAG_*` define.

| Criterion | What is measured | Reference | Limit |
|---|---|---|---|
| accel | \|mean a\| in every pose | the calibration's gravity target `g_cal_local` (normal gravity at the calibration site) | ±10 mg |
| accel noise | per-axis standard deviation at rest | | 5 mg |
| gyro bias | mean rate at rest (residual bias after calibration) | 0 | 0.2 °/s |
| gyro noise | per-axis standard deviation at rest | | 0.2 °/s |
| gyro scale | gyro scale factor fitted so the integrated rotation between two poses carries the first pose's gravity direction onto the second's (needs two fits) | the accelerometer | ±1.5 % |
| field | F and dip in every pose | WMM-2025 at the configured place and the build date | F ±5 %, dip ±3° |
| \|m\| in motion | standard deviation of \|m\| over all readouts, with the field seen along at least 4 of the 6 body directions | constant | 2 % |

`na` means not yet measured, for example before the first still pose.

**Accelerometer.** In every pose, the S line also gives |a| without the calibration (`a_raw_error_mg`). The difference shows what the calibration corrects.

**Gyro.** The bias line shows the uncalibrated bias next to the calibrated one. Each W line carries:
- the rotation and how much the gravity and field directions changed;
- the prediction error with the gyro as measured;
- the fitted scale error, from gravity and from the field direction. The field covers turns about the vertical, which the accelerometer cannot see.

**What the saved calibration scales to.** The magnetometer fit scales the field to the strength the sensor itself measured where it was calibrated (`mag_cal_field_uT` in the H line). So F in a pose matches that number, not necessarily WMM: if both are some percent below WMM, the field where you calibrated was that much weaker, or the raw sensor reads low. `F_raw_uT` in the S lines is the uncalibrated strength. The H line also shows the gyro bias the calibration subtracts and the temperature it was measured at; the S lines show the bias left over (`gyro_bias_dps`) next to the uncalibrated one (`gyro_raw_bias_dps`).

**Magnetometer.** The expected field comes from the World Magnetic Model WMM-2025 ([`src/imu_calibrate/GeomagWmm.h`](../../../src/imu_calibrate/GeomagWmm.h)). It is evaluated at the calibration site (`CalibrationSite` in `AtomS3R_ImuUnits.h`, Fair Lawn NJ) and the build date.

In each pose the sketch reports:
- the total field F;
- its vertical component Z (along gravity, positive down);
- the horizontal component H;
- the dip angle.

It reports each against the model. The model omits local disturbances. Indoors, near steel, vehicles or electronics, F and dip can be off by several percent and degrees, even with a perfect calibration. What a calibration controls is consistency across poses (the dip should not change between poses) and |m| in motion.

To use another place, do either of these:
- build with `-DFIFO_DIAG_LATITUDE_DEG=… -DFIFO_DIAG_LONGITUDE_DEG=… -DFIFO_DIAG_HEIGHT_M=…`, where height is above the WGS84 ellipsoid;
- send `L <lat> <lon> [height_m]` on serial.

For a full magnetometer investigation (axis signs, a fresh calibration fit, magnetometer/accelerometer frame alignment), use [`atomS3R_mag_diagnostics`](../atomS3R_mag_diagnostics/README.md).

## Output

USB serial carries one `key=value` line per record:

| Line | When | Content |
|---|---|---|
| `H` | boot, every 5 s | sketch, build date, which calibrations are applied, and what they hold: accelerometer and gyro calibration temperatures, the gyro bias subtracted (and its temperature slope), the field strength the magnetometer fit scales to |
| `G` | boot, every 30 s, after `L` | place, date, expected field (F, H, Z, dip, declination) and gravity |
| `F` | boot, every 5 s | the FIFO and AUX configuration read back from the BMI270 |
| `D` | 1 Hz | FIFO fill, frame-type counts, lost frames and grid restarts, sensor clock, still state, live \|a\| and \|m\|, temperature, gyro configuration |
| `E` | when frames are counted lost | the read that showed it: frames counted, sensortime, FIFO bytes, whether in warm-up |
| `T` | 10 s | timing criteria |
| `X` | each delay fit | delay, fit width, residuals, reason when rejected |
| `S` | each still pose | duration, tilt, \|a\| error (calibrated and raw), noise, gyro bias, F/Z/H/dip and their errors |
| `W` | each turn between poses | rotation, direction changes, prediction errors, fitted gyro scale |
| `A` | 10 s | accuracy summary over all poses and the criteria |
| `R` | every frame, after `q` | sequence, time (µs), calibrated a (m/s²), ω (rad/s), m (µT, `nan` when stale), temperature |

R lines wait in a queue while USB is busy, so a display burst delays them but does not drop them.

Serial commands, one per line:
- `s`: a 1.5 s display burst;
- `r`: restart the timing and accuracy evaluation;
- `q`: start or stop the R lines;
- `n`: toggle the gyro's noise-performance mode (GYR_CONF bit 6; the INS sketches run with it off, as M5Unified sets it) and restart the evaluation. Compare the `gyro_noise_dps` of still poses in each mode;
- `L <lat> <lon> [height_m]`: set the place.

## Settings

- `FIFO_DIAG_AUX_OFFSET` (default 0) is the BMI270 AUX trigger-to-readout offset, in 2.5 ms steps. 0 keeps a full 40 ms between the forced-mode command and the readout. A shorter offset lowers the delay, but must stay above the 28.3 ms conversion time of the low-noise setting; the data-ready count shows when it does not.
- `FIFO_DIAG_LOAD_EVERY_MS`, `FIFO_DIAG_LOAD_MS` and `FIFO_DIAG_HEAVY_LOAD_MS` set the display bursts.
- `FIFO_DIAG_RAW_LINES` (default 0) starts with the R lines on; `FIFO_DIAG_GYR_CONF` (default `0xA9`) sets the gyro configuration at start (`0xE9` for noise-performance mode).
- `FIFO_DIAG_LATITUDE_DEG`, `FIFO_DIAG_LONGITUDE_DEG`, `FIFO_DIAG_HEIGHT_M` and `FIFO_DIAG_YEAR` set the place and date for the field model.
- `FIFO_DIAG_ACCEL_ERROR_MG`, `FIFO_DIAG_ACCEL_NOISE_MG`, `FIFO_DIAG_GYRO_BIAS_DPS`, `FIFO_DIAG_GYRO_NOISE_DPS`, `FIFO_DIAG_GYRO_SCALE_PCT`, `FIFO_DIAG_FIELD_ERROR_PCT`, `FIFO_DIAG_DIP_ERROR_DEG` and `FIFO_DIAG_MAG_NORM_STD_PCT` set the accuracy limits.

## Host test

`tests/imu_calibrate/fifo_diagnostics_sketch-test` compiles this sketch unchanged and runs it against a simulated BMI270 and BMM150 with known errors in two runs.

- **Good run:** small errors (+3 mg, 0.1 °/s bias, +1 % gyro scale, +2 % field). Every timing and accuracy criterion passes, and each error is measured at its size.
- **Bad run:** errors beyond the limits, at a place set over serial. The matching criteria fail.

[Back to basic IMU examples](../README.md)
