# Magnetometer timing calibration

## Acquisition

The AtomS3R build pins M5Unified 0.2.13 (`a6256725481f1bc366655fa48cf03b6095e30ad1`).
Its `BMI270_Class::begin` sets BMM150 normal mode to 30 Hz and leaves BMI270
`AUX_CONF` at reset `0x46` (25 Hz). The project keeps those settings and sets
47 XY / 41 Z repetitions. Bosch's conversion-time formula gives
`145*47 + 500*41 + 980 = 28295 us`. This is a conversion aperture, not an
instantaneous vector or a known group delay. Neither the acquisition API nor
the board code exports a conversion timestamp, and the BMM150 DRDY pin is not
serviced.

With Bosch compensation active, `readMagBody_` reads the eight BMI270 AUX data
bytes (BMM150 X, Y, Z and RHALL) in one burst. Status bits are masked, and a
reading is new when those bytes differ from the previous burst. XYZ and RHALL
from the same burst are compensated together. Equal payloads cannot be told
apart, so a change detects a new reading, not a unique conversion. The source
repeats its last reading for at most `MagAcquisition::MAX_AGE_MS` (200 ms) and
returns NaN after that.

`sample_us` is the host time before `M5.Imu.update()`, not the sensor sample
time. Gyro and accelerometer values share that timestamp. For each new
magnetometer reading, `MagAcquisition` records the host frame time and the AUX
read bracket (`frame_us`, `read_start_us`, `read_end_us`) together with a
sequence count. These are host times, not BMM150 conversion times.

The estimators use the latest calibrated reading at the current attitude.
OU-II, OU-III and TFG apply a discrete correction when their 35 ms spacing
gate admits a reading. NLO, PII and the standalone Mahony compass use
continuous direction feedback on every IMU tick. No filter compensates for
magnetometer delay.

Primary sources: [pinned M5 driver](https://github.com/m5stack/M5Unified/blob/a6256725481f1bc366655fa48cf03b6095e30ad1/src/utility/imu/BMI270_Class.cpp),
[BMM150 datasheet, sections 4.2–4.3 and 4.6](https://www.bosch-sensortec.com/media/boschsensortec/downloads/datasheets/bst-bmm150-ds001.pdf),
[BMI270 datasheet, AUX interface and AUX_CONF](https://www.bosch-sensortec.com/media/boschsensortec/downloads/datasheets/bst-bmi270-ds000.pdf).

## Dynamic timing calibration

After the gyro and magnetic spatial calibration, the wizard can record two
separate motion segments. Each lasts at least 14 seconds and extends to
45 seconds when more evidence is needed. Tilt the device left and right in a
comfortable back-and-forth rhythm (about one rocking cycle per second), away
from metal. The display shows the usable pair count and what currently limits
the fit. A failed fit can be retried or skipped, or the wizard aborted;
skipping saves the spatial calibration without a timing record. The stage is
skipped when the gyro calibration or the BMI270/BMM150 configuration
fingerprint is unavailable.

The wizard keeps the fitter (`imu_cal::MagTimingCal`) and its result in fixed
wizard storage, not on the loop stack. Gyro samples go into a
`ocean_imu::magnetic::RotationHistory` (`src/util/MagneticRotation.h`). This
ring buffer integrates calibrated rates between host timestamps, handles
`micros()` wrap, and refuses gaps, out-of-order samples and stale queries.

`MagTimingCal.h` compares normalized calibrated fields, rotating them by the
gyro rotation between the two observation times shifted by each candidate lag
(0..80 ms in 2 ms steps). Intervals are 120..180 ms long and do not overlap.
A fit is accepted when all of these hold:

- at least 50 pairs per segment;
- an informative interior minimum;
- training and validation lags within 8 ms of each other;
- residual RMS below 0.05 rad;
- at least 5 % held-out improvement over zero lag.

These are engineering acceptance thresholds, not hardware-qualified confidence
bounds. The fitted lag is relative to the host gyro timestamps, not an absolute
conversion time.

## Storage

The v5 calibration record (`ImuCalBlobV5`, key `blob_m5v5`, 372 bytes) holds a
`MagTimingRecord`:

- the fitted lag, profile width and residual-to-time allowance;
- training and validation RMS, improvement and pair counts;
- the BMI270/BMM150 configuration fingerprint;
- a binding to the gyro and magnetic calibration and the device identity;
- its own CRC.

An all-zero record means no timing was saved. At boot the record is loaded only
when the device, magnetometer source and configuration match. It is reported as
`[MAGTIME CAL] valid=... delay_ms=... allowance_ms=... applied=0`.

**The estimators do not apply the record.** A missing, skipped or incompatible
record never blocks estimator operation.

## Tests

- `tests/imu_calibrate/mag_timing_cal-test` covers injected 20/40/60 ms lags,
  asynchronous jitter, noisy single-axis rocking with pauses, field-axis
  unobservability, validation-only retry, held-out disagreement, clock wrap,
  and v5 storage of both absent and present records.
- `tests/imu_calibrate/mag_rotation-test` checks the rotation history
  contracts with Eigen allocation disabled.

These are synthetic tests, not AtomS3R measurements.
