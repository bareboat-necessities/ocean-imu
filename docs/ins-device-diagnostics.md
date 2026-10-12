# INS device diagnostics, capture and replay

The five AtomS3R INS sketches (OU-II, OU-III, TFG, NLO, PII) can stream what
the estimator actually received and did on the device, so that device behaviour
can be attributed and replayed instead of guessed from synthetic tests.

## Build

Set `SEA_STATE_SERIAL_DIAG` to `1` near the top of the sketch (or pass
`-DSEA_STATE_SERIAL_DIAG=1`). The default is `0`, and with `0` the sketch is
unchanged. With `1` the USB serial port carries the diagnostic stream instead
of NMEA or plotter output. Display, calibration and estimator behaviour are
unchanged. At boot the sketch waits up to five seconds for a serial reader to
attach, so a capture started before a reset sees the stream from its first
frame.

## Stream

Each line starts with a one-letter tag. Every line is formatted in full and
written only when the USB buffer can take all of it, so a slow reader drops
whole lines, never parts of them. The `lines_dropped` count in each D line
reports such drops. Format: `src/AtomS3R/AtomS3R_InsDiag.h`.

| Tag | When | Content |
|---|---|---|
| `R` | every accepted IMU frame | frame sequence number, host sample time, the exact `readImuMapped()` accelerometer/gyro/magnetometer/temperature values (IEEE-754 bits) and the magnetometer acquisition record |
| `D` | 10 Hz | loop health (polls, frames, polls with accelerometer-only/gyro-only/no data, loop overruns, longest loop body, shortest/longest frame interval, peak gyro rate, frames above 30 deg/s) and estimator state |
| `M` | 1 Hz | the shared magnetic-input timing line (`[MAGTIME]`) |
| `H` `I` `C` | after boot, then every 5 s | sketch name and stream format; BMI270/BMM150 output data rates, ranges, AUX and BMM150 repetition registers as read back from the chip; the calibration blob the sketch loaded and the active magnetometer source |
| `E` | on filter reset | `reset_fusion` (boot, recalibration, erase) |

D-line fields common to all sketches include `roll`, `pitch`, `hdg`,
`hdg_ok`, `heave`, `heave_clean` and `hdg_cmp`. `hdg_cmp` is the
tilt-compensated magnetic heading computed from the same published roll and
pitch and the calibrated field. A large `hdg`/`hdg_cmp` disagreement means the
published heading left the magnetometer, for example a north/south swap. Per
estimator:

* OU-II, OU-III, TFG: `live`, `north` (magnetic reference locked), `refined`,
  `tuner`, the startup gravity gate (`gate_good_s`, `gate_aligned`,
  `gate_trusted`) replicated read-only from the filter's own inputs and
  configuration, attitude one-sigma split into `yaw_sd`/`tilt_sd` (sampled
  after the IMU update, before any magnetic correction), `reseeds`
  (sudden attitude-covariance re-seeds such as a tilt-watchdog re-lock),
  magnetic corrections offered/applied, `env`, `wave_hz`.
* PII: `seeded`, Mahony accelerometer `trust` and gain, the vertical
  observer's input and states (`a_vert`, `a_f`, `v`, `p`, `S`, `d`, `b`) and
  gains, and the envelope estimator's internals.
* NLO: `init`, `bootstrap`, `theta`, gyro-bias magnitude, wave-frequency lock.

The stream is about 20 kB/s at the device's accepted frame rate, well within
USB CDC.

## Capture

```
pip install pyserial
python3 tools/ins_diag_capture.py --port /dev/ttyACM0 --out bench.log
```

Start the script, then reset the device. The script reopens the port after the
reset and prints frames per second and any missing frames every five seconds.
Any serial logger that writes raw lines to a file also works. A capture for
replay must not include taps on the screen, because recalibration and erase
reset the filter.

A useful bench session, with a pause of about 30 s between steps:

1. 10 minutes still on the bench.
2. A slow 360 degree turn in about 30 s, flat.
3. Heel to about 80 degrees, hold 10 s, back; then the other side.
4. Pitch to about 90 degrees (on end), hold 10 s, back.
5. A full 180 degree roll flip, hold 10 s, back.
6. Hand rocking about 20 degrees at about 0.5 Hz for 2 minutes.
7. Rough handling for 30 s, then 2 minutes still.

Note the clock time of each step, or the D-line `t_ms` when it starts.

## Report

```
python3 tools/ins_diag_report.py bench.log [--csv state.csv]
```

For each boot, the report lists:

* the sensor rates read back from the chip and whether a calibration was
  loaded;
* frame interval statistics, frame gaps longer than the magnetic rotation
  history accepts (20 ms), and missing lines;
* magnetometer acquisition cadence and how many magnetic corrections were
  applied;
* time to Live and to north lock, and how often the startup gravity gate
  trusted the tilt before Live;
* attitude re-seeds and large-tilt episodes, with the heading disagreement
  that followed them;
* north/south-scale heading disagreements, and heave, envelope or observer
  values that stopped changing while frames kept arriving.

`--csv` writes the D lines for plotting.

## Replay

Each INS test directory builds `ins_sketch_replay-test`. It compiles the
unchanged sketch with the diagnostic stream on, against the host stand-ins in
`tests/common/ins_sketch_stubs/`:

```
cd tests/kalman_ou_iii
make -f Makefile -f ../common/InsSketchHeaveDrift.mk ins_sketch_replay-test
./ins_sketch_replay-test bench.log > replay.log          # replay a capture
./ins_sketch_replay-test                                  # round-trip self-test
```

The replay loads the captured calibration blob and magnetometer source. It
then feeds every R frame into the sketch at the frame's own host time and
prints the replayed sketch's D/M/H/I/C lines. Those lines can be reported or
compared with the device's. The capture must come from the same sketch, or
from the FIFO timing sketch (below). Replay stops at a missing R line
(`--allow-gaps` continues, inexactly) and at a filter reset logged after boot.
`--segment N` selects a later boot in the same file. `--expect-ready` fails
unless the fused output is ready at the end. The summary on stderr names the
first ready frame.

The replay reproduces the device's sensor inputs, magnetometer acquisition
records, calibration and frame times bit for bit. Two things differ. Host time
inside a frame equals the frame time, whereas the device spends some hundred
microseconds there. Polls that delivered no frame are not replayed. The
self-test captures a synthetic 240 s run, parses it back and replays it.
Every frame's published outputs must match bit for bit. It runs in each INS
family's `run_tests.sh`.

## FIFO acquisition captures

`sensors/imu_basic/atomS3R_imu_fifo_timing` reads the IMU the way the INS
sketches are meant to: 200 Hz through the BMI270 FIFO, sample times from the
sensor clock, and a BMI270-triggered magnetometer (`docs/magnetic-timing.md`).
It writes the same stream, with `acquisition=fifo` on its H line. Its R lines
are every delivered frame in `readImuMapped()` units, so a capture from it
replays through any INS sketch:

```
./ins_sketch_replay-test fifo.log --expect-ready > replay.log
```

This shows what each estimator does with the FIFO acquisition before any INS
sketch is changed. The sketch also writes:

- an F line: the FIFO setup read-back;
- X lines: magnetometer-delay fits;
- P lines: its pass criteria.

`ins_diag_report.py` summarises all three. The OU-III test directory runs the
sketch against a simulated BMI270 (`fifo_timing_sketch-test`). It then
replays that capture through OU-III and requires its fused output to be ready
at the end.

