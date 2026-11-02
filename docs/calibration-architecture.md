# Unified calibration

The implementation evolves the existing capture, fitting, runtime and NVS paths. Intrinsic matrices act before physical rotations. With column vectors, the runtime output is in vessel coordinates:

```
a_V = R_VH R_HI S_a (a_raw - b_a(T))
w_V = R_VH R_HI S_g (w_raw - b_g(T))
m_V = R_VH R_HI R_IM A_m (m_raw - b_m)
```

`I` is the accelerometer intrinsic frame, `M` the magnetometer intrinsic frame, `H` the housing frame, and `V` the vessel frame. The fixed driver axis permutation and nominal-g unit conversion precede these expressions. Intrinsic accelerometer/soft-iron matrices do not encode mounting orientation. The full gyro matrix maps its calibrated triad into `I`; a rate-only rig is not required when diverse gravity-referenced rotations are available.

## Reuse assessment

| Component | Disposition | Mathematical behavior and coverage |
| --- | --- | --- |
| `AccelCalCapture.h`, `AccelCalFit.h` | Retain | Multi-pose quiet blocks, robust PolarSPD fit, conditioning/plausibility gates, temperature qualification, cross-validation and float revalidation; `accel_cal-test`, replay and gravity-override tests. |
| `CalibrateIMU.h` | Extend | Preserve bias/thermal and magnetic fit utilities. Populate the existing gyro `S` matrix from gravity-referenced rotation segments. |
| `GyroCalCapture.h` | Extend | Preserve stale-word, time-gap, scatter, gravity and independent magnetic stillness checks. Expose qualified endpoint means for rotation fitting. |
| `MagCalFit.h`, `MagCalCapture.h` | Retain | Robust hard/soft iron, independent verification and field-quality diagnostics remain authoritative. |
| `AtomS3R_MagDiagnostics.h::solveAlignment` | Extend | Preserve iterated constant-dip rotation fitting; reject unobservable geometry and validate a frozen result on separate poses. A ridge alone is not evidence of observability. |
| `AtomS3R_ImuCalWizard.h` | Refactor/extend | Reuse menus, retries, quiet capture, worker tasks, progress and save verification. Add gyro rotation and alignment stages. |
| `AtomS3R_ImuCalBlob.h`, Preferences backend | Refactor/extend | Preserve CRC, verified writes and verified rollback of current-schema records. New schema/key rejects old records; no migration or legacy fallback. |
| Estimator measurement updates | Extend | Keep each estimator's residual/Jacobian, Joseph update and error injection. Share the magnetic NIS decision, inflation and health state machine. |
| Stationary noise | Extend | Separate 40-minute acquisition and Allan analysis; never infer a marine-motion prior from stationary sensor data. |

Existing tests remain in place unless a schema assumption explicitly requires updating. The modified magnetic-service policy requires its own stability qualification; existing certificates are not automatically evidence for this version.

## Gyro matrix and independent validation

The existing quiet-interval capture supplies bias and qualified endpoint means.
Twenty subsequent turn-and-hold segments provide 15 training rotations and 5
validation rotations (every fourth segment). Timestamp gaps over 20 ms,
non-finite values, excessive rates, exhausted storage, and unqualified endpoints
discard the active segment. A retry establishes a fresh quiet start.

After each rotation prompt, keep the device still until **Turn now** appears.
The existing stationary capture qualifies a fresh start after button/menu
handling. Integration starts at that fresh sensor timestamp. You may wait
before turning; an unchanged pose is not a failed rotation. Proven quiet
prefixes are released while waiting, retaining the six-second endpoint means
and every unqualified/moving interval. A gap before motion restarts start-pose
qualification; a gap during motion rejects the segment. The retry screen and
`[GYRROT]` serial diagnostics identify the actual failure (sample gap, stale
stream, insufficient angle, or unsuccessful stillness).

For calibrated unit specific force u, inertial constancy gives
du/dt = -[omega]x u. Each measured midpoint rate is corrected by the existing
temperature-bias model, then integrated with Exp(-S_g omega dt). A robust Huber
objective minimizes endpoint direction errors over training segments. Analytic
SO(3) sensitivities propagate through the same exponentials. The optimizer
uses Eigen's existing fixed-size LDLT; no optimization dependency is added.
The qualified raw stationary gyro mean seeds the first midpoint interval;
the integration does not discard the time between the start pose and its
first fresh rotation sample. Finite spin-up trajectories test this coverage.
The undamped 9x9 information matrix must be observable and conditioned.
Singular values, determinant, and coefficient bounds reject implausible
solutions. A bias-only thermal model requires temperature to stay within 2 C
of its measured reference during the rotation acquisition.

Final float coefficients must pass the frozen, independent rotation segments:
training and validation RMS below 0.025 rad and maximum validation error below
0.06 rad. A nearly factory-perfect identity is acceptable; nonidentity is not
forced. This is the gravity-referenced IMU-TK method adapted to the existing
accelerometer frame and capture machinery.

## Physical registration

The 21 quiet endpoint poses also provide raw magnetic means. Sixteen poses fit
the existing constant-dip alignment solver; five independent poses validate
the frozen rotation and dip. The undamped information must have full rank
before the original numerical ridge is added. A repeated pose, an unobserved
axis, a rotation over 20 degrees, or independent dip RMS over 1.5 degrees fails.
A large wrong signed-axis mapping must be corrected in the driver.

Mounting uses known physical directions. The retained Eigen SVD solves Wahba's
problem in solveFrameAlignment() for fixture tooling and rejects collinear
references: gravity alone cannot determine yaw. At the MOUNT screen, tap only
when the documented IMU, housing and vessel axes actually coincide (forward,
starboard, down). Otherwise send measured proper rotations:

~~~
MOUNT <R_HI: nine row-major numbers> <R_VH: nine row-major numbers>
~~~

The rotations are validated separately and composed only by RuntimeCals.
A full intrinsic update requires renewed gyro/magnetic registration; the
retired accelerometer-only UI could otherwise retain dependent qualifications.
The retained library candidate helper invalidates those dependent qualifications.

## Separate stationary noise workflow

After a complete current-schema calibration exists, the start menu offers a
separate 40-minute stationary capture. Warm up the sensor and use a quiet,
thermally stable support. Checks cover fresh timing, frozen words, finite
data, acceleration/rate limits, fixed one-second pose means, fixed magnetic
means including loss/staleness, and temperature span below 2 C. At 200 Hz,
480001 samples span 2400 seconds. Missing frames invalidate the run;
samples are never concatenated across gaps.

The fixed-memory implementation computes non-overlapping-cluster Allan
deviation at 17 dyadic averaging times. Adjacent differences share a cluster
mean; pair counts are not independent degrees of freedom. Identification
requires four consecutive dyadic levels (about 0.9 decade), at least 32
differences at the longest averaging time, a matching log slope and bounded
fit residual. White noise, a bias-instability-like plateau, and rate random
walk use slopes -1/2, 0, and +1/2. Long-term regions must start at tau >= 1 s.
A compatible finite-data slope is evidence over the observed range, not proof
of flicker noise or a bound on arbitrarily long bias histories. Forty minutes
often cannot identify long-term terms; absent terms are marked unidentifiable.

The white coefficient N uses Allan deviation N/sqrt(tau). Continuous intensity
is N^2; independent sample variance is N^2/dt and one-sided spectral amplitude
density is sqrt(2) N. Bias instability uses plateau/0.664. Rate random walk uses
K = adev sqrt(3/tau). Bias instability is never substituted for random-walk Q.

The record configures measured accelerometer measurement scatter and gyro
white process-noise density on every reset, plus gyro-bias walk only on
identified axes. Unidentified terms retain the explicitly configured
estimator model. No old calibration record is consulted. OU/CARMA marine
frequency, variance, poles and pseudo-measurements remain separate.
The standalone gravity-only qMEKF compass retains its configured allowance
for unmodeled translational acceleration: its accelerometer covariance is
measured sensor covariance plus that separate model covariance. The full INS
variants already represent marine acceleration in their state model.

Build with ATOMS3R_NOISE_RAW_LOG=1 for a complete serial raw trace. A dropped
line aborts that diagnostic acquisition; keep the host draining USB. Ordinary
capture needs no host and reports the Allan table and identified terms.
Reproduce the six-axis analysis with the same implementation:

~~~
make -C tests/imu_calibrate sensor_noise-replay
tests/imu_calibrate/sensor_noise-replay capture.log > noise.json
~~~

NOISERAW lines contain uint32 microseconds, calibrated accel xyz, calibrated
gyro xyz and temperature. Device-side magnetic stillness must also pass;
the six-axis replay does not reconstruct that independent check.

## Shared magnetic update policy

Each filter retains its measurement model and tangent coordinates. The shared
gate uses r' (H P H' + R)^-1 r with the full positive-definite 3x3 R.
OU-II/III rotate body-frame covariance into their de-heeled residual frame;
TFG variants rotate it into world coordinates. Standalone qMEKF's compass
consumer uses the full calibrated uT vector; unit normalization would remove
a measurement degree of freedom.

Three-dimensional chi-square quantiles are 11.34486673 (99%) and 21.10751347
(99.99%). A moderate innovation inflates R, preserving its correlations and
leaving H P H' untouched. A gross innovation or field-norm anomaly is rejected.
Rejection changes neither state nor covariance and never resets the magnetic
reference. Ten clean innovations at the ORIGINAL covariance are required
before recovery; recovery starts with inflated R and decays it gradually.
Inflated covariance cannot qualify its own recovery. Gains, Joseph updates
and quaternion/Lie-group resets use the effective covariance.
magneticHealth() exposes nominal NIS, inflation, health and cumulative counts.
The standalone qMEKF consumer adapts its dip reference slowly only AFTER a
healthy accepted update. Rebuilding the reference from the current reading
before the gate would let a disturbance erase its own innovation.

The existing OU/CARMA and TFG reference refinements also run after the current
packet's Kalman update. Replacing a reference that has supported applied
magnetic service, or assigning its absolute yaw, requires a full healthy,
accepted window. A rejected packet cannot bypass NIS through reference learning.
At healthy recovery a disturbed partial window is discarded, including its
norm mean, so that old disturbed samples cannot veto the returning clean field.
Independent strength/dip repair remains available for a provisional reference
that has never supported an accepted magnetic update. The independent proxy
tilt and original averaging/fitting algorithms are retained.

## Persistence and resource validation

Only schema 6, key blob_m5v6, is accepted. The 680-byte record retains explicit
padding, CRC, device identity, coefficient bindings, verified writes, exact
readback and verified recovery. Version 5 records, even under the new key,
cannot be read or automatically relabeled by save(). Runtime rebuilding also
fails closed. The boot gate requires fresh calibration for unsupported records.

The existing runImuCalWizard entry is defined out of line so large fitting
and GUI instantiations compile separately from estimator sketches. The compass
uses that same entry instead of retaining another wizard object. This reduces
compiler peak RAM and keeps the large wizard off the application stack.
The public runAndSave dispatcher runs the complete UI, capture and validation
flow on a dedicated 32 KiB task on the caller's core. Moving only the object
to the heap does not protect Arduino's smaller loopTask stack from nested
calls. The caller waits and services its watchdog when registered; only the
wizard task accesses the UI during that wait. Serial logs report the task's
stack budget and actual free-stack high-water mark.

The rotation buffer is reserved explicitly in PSRAM before the first full
calibration capture and released after fitting; it never occupies the loop-task
stack. AtomS3R builds require **Tools > PSRAM > OPI PSRAM**, followed by a
rebuild/upload (`PSRAM=opi` in the FQBN). The Arduino AtomS3 board defaults to
PSRAM disabled. The wizard reports the configured mode, detected/registered
external heap, total/free bytes, largest block and required allocation. Missing
PSRAM produces a setup instruction before any pose captures; low or fragmented
memory produces a separate allocation error. The same reservation is reused
through capture, rather than freed and allocated again. Separate noise capture
does not reserve the gyro trajectory. M5Stack's [AtomS3R configuration](https://docs.m5stack.com/en/core/AtomS3R)
uses `qio_opi` with `BOARD_HAS_PSRAM`.
The UI and fitting tasks each have a 32 KiB stack, with up to
64 KiB allocated simultaneously during a fit. Stack sizes and high-water marks use bytes, as
required by ESP-IDF, rather than upstream FreeRTOS words. Logs report allocated
bytes, free PSRAM, fit microseconds and stack headroom. These measurements
must be collected on hardware; host timing is not an ESP32-S3 benchmark.
The separate atomS3R_calibration_benchmark sketch reports synthetic gyro-fit,
noise-update and magnetic-gate timing plus actual task stack headroom; it
never reads or writes calibration records. Compile it with the same OPI PSRAM
FQBN as the production sketches and retain its BENCH serial records:

```sh
arduino-cli compile --library . --jobs 1 \
  --fqbn esp32:esp32:m5stack_atoms3:CDCOnBoot=cdc,USBMode=hwcdc,PSRAM=opi,PartitionScheme=default_8MB \
  sensors/imu_basic/atomS3R_calibration_benchmark
```

Use this partition scheme for production builds as well. The board has 8 MB
flash; default_8MB provides two 3,342,336-byte OTA slots and 1.5 MB SPIFFS.
The full OU-III image exceeds the default 1,310,720-byte application slot.

Existing accuracy, gravity, quiet capture, magnetic fitting, storage fault
and replay tests remain. New tests cover SO(3) derivatives, gyro matrix
recovery, degenerate geometry, independent rotations, a 40-minute noise
simulation, unsupported records, noncommuting frame composition, fixture
yaw observability, correlated NIS, covariance PSD, and a 90-degree magnetic
heading disturbance across six Kalman implementations.
Run the repository-required make all. Build serially on constrained hosts:
simultaneous large CARMA template compilations can exceed host RAM.

## References

- Tedaldi, Pretto and Menegatti, ICRA 2014, DOI [10.1109/ICRA.2014.6907297](https://doi.org/10.1109/ICRA.2014.6907297): gravity-referenced multi-position gyro calibration.
- [INSLIB](https://github.com/jnz/INSLIB): reference implementation and ideas, not a replacement for the existing accelerometer fitter.

