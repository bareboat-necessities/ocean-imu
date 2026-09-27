# Calibration safety and registered sensor axes

External sensor calibration is applied once, before fusion. Sensor acceleration
uses `g_std = 9.80665 m/s^2` for nominal-g conversion. The separately configured
Fair Lawn `g_cal_local`, saved `accel_g`, and gravity-mismatch rejection are
unchanged. Accelerometer calibration retains ten approximate guided poses,
time blocks, a robust full symmetric SPD fit, information and held-out gates,
three later verification holds, float/runtime revalidation, and qualified
compatible thermal-prior handling. No exact right-angle poses are assumed.

## Accelerometer pose progress

The bar counts accepted 250 ms blocks from the start of a pose. Placement,
motion and rejected samples do not fill it. A full bar means that all five
useful seconds have been retained, plus the two independent verification
seconds for a later recheck. Previously retained progress stays visible while
the device settles again after a disturbance. The hints indicate placement,
motion or a pose that needs correcting; an unsuccessful hold still times out
after 30 seconds. A successful hold shows `Captured` and advances after the
brief acknowledgement (about one second).

Normal accelerometer diagnostics are best effort: whole lines are skipped
when they do not fit in the available serial transmit buffer. The burst of
block records at pose completion therefore cannot wait for USB buffer space.
For a complete raw/block replay capture, build with `ATOMS3R_ICAL_RAW_LOG=1`
and keep the serial host reading. That explicit diagnostic mode uses blocking
output and can affect capture timing; ordinary calibration does not require
a serial host. Calibration observations and quality gates do not depend on
whether diagnostics are delivered.

## Gyroscope stillness and temperature qualification

The device asks the user to put it on a table. It waits for three quiet
250 ms blocks, then collects six continuous quiet seconds at at most 40 Hz.
Progress resets automatically after motion or a sample gap over 60 ms; samples
from the interrupted interval are discarded. This uses gyro and accelerometer
scatter and changes in their block means, not just small angular-rate magnitude.
Once magnetic readings are available, changes in their one-second means also
reject slow yaw. Only distinct magnetic register readings enter those means;
polling a slower magnetometer at the IMU rate does not create new observations.
The first full magnetic mean remains the reference for the entire hold, and
capture waits for it before retaining gyro samples. This averages stationary
magnetic noise without letting a rolling reference follow slow rotation.
Losing or freezing that magnetic stream cannot disable the check
mid-hold. Repeated IMU register values cannot qualify a hold. The stage times
out after 70 seconds and offers a retry without repeating the other stages.
Each restart logs its cause over serial as `[GYR] restart=... reason=...` so a
sample gap, sensor scatter, field change, or actual motion can be distinguished.

The default block limits are 0.005 rad/s total gyro standard deviation,
0.004 rad/s change in mean gyro, 0.08 m/s^2 total acceleration standard
deviation, and 0.06 m/s^2 change in mean acceleration. Magnetic change is
limited to the larger of 0.6 uT and 1.5% of the reference raw field magnitude.
Direct narrow-temperature gyro fits also reject per-axis scatter above
0.006 rad/s. These are acceptance policies, not guarantees of zero motion:
constant rotation about gravity without a useful independent magnetic reference
is indistinguishable from gyro bias, and sufficiently small motion can remain
below the noise thresholds. Calibration needs a stationary support.

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
rotation. The correction uses the SPD square root as an initializer, then
refines a symmetric positive-definite matrix with the fitted hard-iron offset
subtracted first. A Cholesky whitening factor gives the right
norm but not the registered-frame direction. Independent physical mounting
misalignment is not observable from these magnitude-only samples and is not
absorbed into the soft-iron correction. Absolute field scale also remains tied
to the fitted radius, not an independently measured geomagnetic reference.

The device flow is **collect -> refine -> check -> save**. Guidance uses the
retained 3D coverage: **Turn slowly**, **Tilt forward/back**, **Roll left/right**,
or **Flip over**. The weakest axis determines the requested movement, with a
three-second prompt interval to keep the screen readable. These are broad
motions, not exact poses or measured 90-degree turns. Level-only rotations
cannot qualify a full 3D calibration for equipment that pitches and rolls.
Insufficient coverage keeps capture open with guidance instead of forcing an
immediate restart. Progress includes elapsed time, sample count and coverage.

The wizard observes freshness independently of sample retention. Meaningful
changes reset the stale timer even after the buffer is full; missing, invalid
or frozen data do not. At most one eligible observation per 80 ms enters the
existing 400-sample buffer. When full, observations in sparse direction cells
replace observations in the most crowded cell; otherwise fixed-seed reservoir
replacement retains later data. The 26 broad cells use centered raw directions
during collection and corrected directions during verification. Retention is
deterministic for a given stream, without requiring a second sample array.

Collection requires at least 45 seconds, 360 retained samples, span ratios
0.35/0.55, two-axis direction range 1.05, unit-direction covariance determinant
0.0002, and at least 12 direction cells with three observations each. Capture
is bounded to 220 seconds. Missing or frozen data fail after 12 seconds without
meaningful change. Coverage calculations update only when retained data change.

### Geometric refinement and acceptance

The normalized robust algebraic ellipsoid fit initializes a bounded
35-iteration double-precision Levenberg-Marquardt refinement of the geometric
residual `norm(A * (m - b)) - B`. A log-Cholesky parameterization keeps `A`
positive definite; the applied matrix is `L * L^T`, not a rotating whitening
factor. The field scale `B` stays fixed at the initial raw-radius estimate to
remove the scale ambiguity. Direction-balanced Huber weights prevent a long
easy rotation or isolated outliers from dominating. Refinement is automatic
and runs on the existing fit task with fixed-capacity workspace.

The final stored-float correction must pass all of these policies:

| Check | Limit |
| --- | --- |
| Matrix eigenvalues / condition | 0.2 to 5 / at most 10 |
| Hard-iron offset norm | At most 150 uT |
| Inlier residual | At most `max(0.6 uT, 0.025 * B)` |
| Inlier fraction | At least 85% overall and in each populated time quarter |
| Inlier RMS | At most `max(0.35 uT, 0.015 * B)` |
| 95th-percentile absolute residual | At most twice the inlier residual limit |
| Corrected 3D coverage | At least 12 populated cells, both signs on every axis, direction-moment determinant at least 0.001 |
| Information | Positive, conditioned nine-parameter information matrix; bias SD proxy at most `max(0.5 uT, 0.02 * B)`, matrix SD proxy at most 0.05 |
| Field consistency in time | Mean inlier residual range across populated time quarters at most `max(0.35 uT, 0.01 * B)` |

Information is weighted as one observation per direction cell with a 0.2 uT
noise floor. Its SD values are conservative identification proxies, not a
traceable accuracy or absolute heading certificate. Time checks use observation
timestamps, not reservoir slot order. Failure reasons identify inadequate
directions, interference, a changing field, or an implausible correction.
`MagCalibration::rms` and the blob's `mag_rms` retain their trimmed-RMS reporting
convention, recomputed from the refined correction. Acceptance uses the separate
full-inlier `quality.rms`, tail, fraction and time gates; the smaller trimmed
report cannot make a poor fit qualify.

### Independent verification

After refinement, a new turn-and-tilt sweep checks frozen coefficients. None of
these observations enter a refit. Verification requires at least 12 seconds,
140 samples, full corrected coverage, and the same residual, plausibility,
information and temporal gates; its timeout is 60 seconds. A changed magnetic
environment or failed check offers a MAG retry before any write. Final candidate
and read-back coefficients are rebuilt through the runtime path and checked
again against this independent sweep. Stored matrix plausibility and RMS are
also checked on ordinary load; the blob layout remains unchanged.

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

`calibration_workflow-test.cpp` exercises slow arbitrary sweeps, uneven and late
3D coverage, coverage-driven prompts, robust refinement, noisy/outlier heading
and tilted-vector accuracy, magnetic field changes, frozen-coefficient
verification, and rejection of poor stored fits. Gyro cases include continuous
quiet capture, steady yaw/roll, the half-moving bias regression, movement after
progress, missing samples, late/frozen magnetic references, and clock rollover.
Stationary magnetic-noise cases use a 200 Hz IMU with 10, 20 and 100 Hz magnetic
readings (0.35 uT per-axis noise) over 30 seeds each, checking uninterrupted
completion and bias accuracy. Separate noisy 10 Hz cases require a new quiet
hold after slow yaw, a field step, or loss/freezing of the magnetic stream.
These are deterministic host regressions; physical device timing, noise and
interaction still require an on-device check.
