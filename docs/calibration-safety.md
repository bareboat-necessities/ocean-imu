# Calibration safety and registered sensor axes

External sensor calibration is applied once, before fusion. Sensor acceleration
uses `g_std = 9.80665 m/s^2` for nominal-g conversion. The separately configured
Fair Lawn `g_cal_local`, saved `accel_g`, and gravity-mismatch rejection are
unchanged. Accelerometer calibration retains ten approximate guided poses,
time blocks, a robust full symmetric SPD fit, information and held-out gates,
three later verification holds, float/runtime revalidation, and qualified
compatible thermal-prior handling. No exact right-angle poses are assumed.

## Gyroscope temperature qualification

A normal short stationary capture estimates the mean accepted angular rate at
the measured mean temperature. Its temperature slope is exactly zero unless
all thermal identification gates pass. Full gyro recalibration does not borrow
an unqualified historical slope; accelerometer-only recalibration carries the
existing validated gyro calibration unchanged.

The thermal regression uses double precision and equally weighted populated
temperature-bin means, centered on their measured mean temperature. Each bin
requires 20 accepted stationary samples. Qualification requires at least four
populated bins, at least 5 degrees C between their temperature centroids, and
centered temperature information `sum((T_bin - mean(T_bin))^2) >= 25 degC^2`.
The three-axis one-sigma slope limit is `0.00005 rad/s/degC` per axis. Scatter
includes both bin lack-of-fit and raw residual scatter, with a `0.0002 rad/s`
floor; raw sample count does not make correlated hold noise artificially precise.
Every axis must satisfy `abs(k) + 3*sigma_k <= 0.001 rad/s/degC`.

The slope ceiling is a conservative rejection policy, not a sensor accuracy
claim. Bosch specifies BMI270 gyro TCO of +/-0.02 deg/s/K, approximately
0.000349 rad/s/K (Bosch Sensortec BMI270 product specifications).

Qualified slopes are applied only within the populated-bin temperature
interval plus a 2-degree margin at each end. Temperatures outside it use the
nearest bound. Nonfinite temperatures use the reference bias without any
thermal increment. Float conversion is checked over the complete permitted
interval; stored qualification, uncertainty, coefficient binding and clamp
limits are validated again on load and runtime reconstruction.

A stationary, temperature-correlated physical drift cannot in general be
separated from a true thermal effect using this capture alone. Qualification
limits that ambiguity; it does not claim a traceable environmental-chamber
calibration. Short captures normally remain bias-only.

## Magnetometer direction and capture

The ellipsoid metric determines `A^T A`, not an arbitrary sensor-to-body
rotation. The correction is its unique SPD square root, with the fitted
hard-iron offset subtracted first. A Cholesky whitening factor gives the right
norm but not the registered-frame direction. Independent physical mounting
misalignment is not observable from these magnitude-only samples and is not
absorbed into the soft-iron correction. Absolute field scale also remains tied
to the fitted radius, not an independently measured geomagnetic reference.

The wizard observes freshness independently of sample retention. Meaningful
changes reset the stale timer even after the buffer is full; missing, invalid
or frozen data do not. At most one eligible observation per 80 ms enters a
fixed-seed Algorithm-R reservoir using the existing 400-sample buffer. This
retains observations across the full interval, including late rotations,
without allocating another sample array. Retention is deterministic for a
given incoming stream, not a guarantee that every rare orientation is kept.

The 45-second minimum, 360 retained-sample goal, span ratios 0.35/0.55,
two-axis direction range 1.05, and unit-direction covariance determinant
0.0002 are enforced on the retained observations. The subsequent robust
ellipsoid fit retains its separate, stricter coverage/inlier gates. A real
freeze or missing data still fails after 12 seconds without meaningful change.

## Persistence and failure reporting

`ImuCalBlobV4` is a 324-byte explicitly laid-out record under `blob_m5v4`.
Its outer CRC covers all metadata. A separate gyro binding covers coefficients,
thermal evidence, slope uncertainty, qualification and runtime bounds together.
Older keys/layouts are not migrated or interpreted: installation requires a
new full wizard calibration. This also prevents loading a legacy rotated
magnetometer whitening matrix as a registered-frame correction.

`saveVerified()` succeeds only when a read-back validates and byte-matches the
candidate. Invalid candidates and failed reads of an existing record do not
initiate a write. Following a failed write, the previous valid bytes are checked
or restored and checked again. The status distinguishes verified retention,
verified empty cleanup, and recovery failure. A short/corrupt rollback or a
failed/falsely successful removal is reported as `RECOVERY_FAILED`; the wizard
shows **Recovery failed**, never an assurance that the previous calibration
was restored. Callers of the lower-level `save()` only get write acceptance,
not a read-back guarantee.

This single-key backend is not a power-loss journal. If storage fails during
both writing and recovery, the previous record can be lost. CRC validation
prevents applying corrupt data, but software cannot promise successful
rollback after a hardware/storage failure. The Preferences API also reports
zero length for both an absent key and a failed namespace open.

## Regressions

`tests/imu_calibrate/calibration_safety-test.cpp` covers constant/narrow and
wide temperature data, information and plausibility rejection, float/NVS
readback and extrapolation, diagonal/off-diagonal SPD and hard-iron distortion,
noise/outliers, direction/heading, capture saturation/late coverage/staleness,
clock rollover, missing temperature, and injected storage/recovery faults.
`test_sketch_temperature.py` checks all external-calibration INS consumers and
the compass path while protecting TFG's separate internal reference-temperature
arguments. Both run in the calibration suite alongside the unchanged
accelerometer campaign, gravity override and capture replay.
