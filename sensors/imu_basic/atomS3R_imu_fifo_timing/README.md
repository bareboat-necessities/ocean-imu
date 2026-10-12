# AtomS3R IMU FIFO timing diagnostics

[ocean-imu](../../../README.md) / [sensors](../../README.md) / [imu_basic](../README.md) / **atomS3R_imu_fifo_timing**

<p align="center">
  <img src="../../../img/devices/AtomS3R_device.svg" width="240" alt="M5Stack AtomS3R">
</p>

This sketch tries out, on the device, the sensor acquisition the INS sketches are meant to move to, and checks that it delivers what the estimators assume. It estimates nothing and never changes the saved calibration. The INS sketches are not changed by it: they still poll the IMU through M5Unified (see [`docs/magnetic-timing.md`](../../../docs/magnetic-timing.md)).

Open the sketch: [`atomS3R_imu_fifo_timing.ino`](atomS3R_imu_fifo_timing.ino).

## What it changes compared with the INS sketches

| | INS sketches today | This sketch |
|---|---|---|
| Accelerometer / gyro rate | 100 Hz / 200 Hz (driver defaults); only frames with both new are used | both 200 Hz |
| Reading | data registers, polled by the loop; a busy loop loses samples | BMI270 FIFO, drained every loop; a busy loop delays samples but loses none |
| Sample time | host `micros()` when the loop polled | the sample's place on the sensor's own clock (a frame every 128 sensortime ticks), mapped onto `micros()` |
| Magnetometer | BMM150 free-running at 30 Hz, read by the BMI270 at 25 Hz, so the age of each reading drifts | BMM150 in forced mode; the BMI270 starts a conversion every 40 ms and reads it 40 ms later, so every readout is one new conversion at a fixed point on the sensor clock |

Ranges stay at ±8 g and ±2000 °/s, the ranges M5Unified's conversion assumes, so samples are in exactly the units `readImuMapped()` delivers (a host test checks this bit for bit).

## Running it

You need a saved calibration (run the wizard in any INS sketch first). Without a magnetometer calibration the delay fit is likely to reject its pairs.

1. Upload, then run `python3 tools/ins_diag_capture.py --port <port> --out fifo.log` and reset the board.
2. Leave the board still for about 30 s.
3. Rock it back and forth by hand, about one rock per second, tilting 20–40° and changing direction. Keep away from metal. Keep going until the screen shows `dly OK` (usually 1–3 minutes).
4. Tap the screen once or twice: each tap blocks the loop with 1.5 s of display drawing. A 400 ms burst also runs every 10 s by itself.
5. Stop the capture and run `python3 tools/ins_diag_report.py fifo.log`.

The screen and the P lines show the five criteria:

| Criterion | Pass when |
|---|---|
| rate | no frame missing on the sensor clock, and the host-clock rate within 1 % of 200 Hz |
| no loss | no frame lost across the display bursts (at least one burst in the window) |
| spacing | every delivered sample spacing (after warm-up) between 4.9 and 5.1 ms |
| mag | 25 Hz readouts, each with the BMM150 data-ready bit set (a new conversion), each 40 ms after the previous on the sensor clock |
| delay | at least three magnetometer-delay fits within 4 ms of each other |

The delay is fitted by the calibration wizard's timing code (`MagTimingCal`) from the rotation the gyro measured. It is expressed in the convention the estimators' magnetic transport uses: a gyro sample is the rate over the interval ending at it. It therefore includes the conversion aperture, the time from the conversion to the readout, and the gyro's own filter delay. That is the number a filter needs.

Serial commands: `s` runs a 1.5 s display burst, `r` restarts the evaluation, `q` stops or restarts the R lines.

## Output

USB serial carries the [diagnostic stream](../../../docs/ins-device-diagnostics.md) (`H`, `I`, `C`, `D`, `M`, `R`), plus:

- `F`: the FIFO and AUX configuration read back from the BMI270 after setup;
- `X`: one line per delay fit (delay, fit width, residual allowance, pair counts, reason when rejected);
- `P`: the pass criteria every 10 s.

The `R` lines are every delivered frame, in the units and format the INS sketches' diagnostic stream uses. A capture from this sketch can therefore be replayed through any INS sketch on a PC (`ins_sketch_replay-test fifo.log`). That shows what each estimator does with this acquisition before any INS sketch is changed. R lines wait in a queue while USB is busy, so a display burst delays them but does not drop them.

## Settings

- `FIFO_TIMING_AUX_OFFSET` (default 0): the BMI270 AUX trigger-to-readout offset, in 2.5 ms steps; 0 keeps a full 40 ms between the forced-mode command and the readout. A shorter offset lowers the delay but must stay above the 28.3 ms conversion time of the low-noise setting; the data-ready count shows when it does not.
- `FIFO_TIMING_LOAD_EVERY_MS`, `FIFO_TIMING_LOAD_MS`, `FIFO_TIMING_HEAVY_LOAD_MS`: the automatic and tap display bursts.

[Back to basic IMU examples](../README.md)
