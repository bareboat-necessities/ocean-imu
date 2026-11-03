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
Requirement-by-requirement implementation and remaining qualification are tracked
in [calibration-requirements.md](calibration-requirements.md).

## Automatic fallback and basic operation

Fresh measured accelerometer calibration, stationary gyro bias, verified
magnetometer intrinsics and physical mounting declarations remain required.
A failed advanced fit no longer discards those measured sets.

The gyro worker first fits all nine coefficients of the existing correction
matrix. If that fails, it automatically refits the same captured rotations
with only three diagonal scales free. It reuses the same integration,
analytic sensitivities, robust objective, physical limits and independent
validation. The three-parameter information matrix must be observable and
conditioned. The optimizer starts from identity; it does not use coefficients
from the failed full fit. No turns are repeated for this attempt. A successful
fallback is stored as `gyro_matrix_ok=1`, `gyro_scale_only=1`, with a qualified
diagonal `S`; cross-axis corrections remain unmeasured. A successful full fit
has `gyro_scale_only=0`. DONE, the GYRO SCALE notice and startup identify the
chosen model. The full-fit baseline accuracy and noisy campaign are unchanged.

At **GYRO TILT**, one tap starts the twelve-turn procedure and two taps choose
**Use bias only** before rotation capture or PSRAM allocation. When both full
and scale-only fits fail, one tap appends four turns when possible (up to twenty)
or retries the final fit. Two taps choose bias only. The same escape is offered
after a failed turn or missing rotation memory. Three taps open a confirmed exit.
Bias-only restores the latest strict stationary fit from this session, including
its thermal qualification. All gyro models retain that strict bias: the relaxed,
hand-supported endpoint means constrain gravity and rotation fitting only. They
are never averaged into the stationary bias. Bias-only sets
`S=I`, `gyro_matrix_ok=0`, `gyro_matrix_skipped=1`, `gyro_scale_only=0` and clears
rotation-fit quality metadata. Earlier accel/mag stages and accepted turns are
kept. The choice has an untimed explanation page and never reuses coefficients
from an earlier record.

Magnetic registration is attempted once on the retained endpoints. On any
failure it automatically continues with nominal driver axes: the already
verified hard/soft-iron correction is retained, `R_IM=I`, `mag_alignment_ok=0`,
`mag_alignment_skipped=1`, and alignment quality metadata is zero. A brief
MAG SKIP notice and serial diagnostic state the result. There is no retry menu
or extra magnetic pose sequence, including when gyro tilt was skipped before
providing endpoint pairs. A successful gyro fit, whether full or scale-only,
is retained. Magnetic timing acquisition is omitted when gyro tilt was skipped
or magnetic registration did not qualify.

DONE and the startup BASIC CAL notice distinguish **Gyro: full matrix**,
**Gyro: scale only**, **Gyro: bias only**, and **MAG: nominal axes** as applicable.
These profiles are usable measured calibrations with declared unmeasured
corrections, not full-fit accuracy claims. Cross-axis errors remain after a
scale-only fit; scale and cross-axis errors remain in bias-only mode. Skipped
magnetic registration leaves PCB rotation error unmeasured. Correct nominal
driver axis mapping remains required. The boot gate accepts these profiles
after verified save/readback, while continuing to reject missing measured bias,
unverified intrinsics, missing mounting, corrupt records or a missing gyro
qualification without an explicit bias-only selection.

The profile flags use reserved bytes in the same 680-byte schema-6 record.
CRC and gyro coefficient binding protect them. A scale-only model must be
qualified and diagonal; a skipped gyro fit must be identity with zero quality
metadata. Skipped magnetic registration must also be identity with zero alignment
metadata. Contradictory profiles are rejected. Old unsupported records are not
migrated or used as fallback.

## Repairing a saved gyro calibration

With a usable current-schema record, **CALIBRATE → two taps: saved fits →
one tap: gyro bias only** takes one strict table-rest capture without repeating
accelerometer, magnetic or mounting stages. The instruction card is untimed;
after the placement window, allow about eight quiet seconds. Existing capture,
stationary fitting and verified NVS writes are reused.

This repair preserves the saved accelerometer, magnetic and mounting fields,
sensor identity and unrelated quality metadata. It replaces the gyro model
with fresh measured bias and nominal scale/axes, explicitly marked as bias-only.
A matrix fitted using an incorrect bias cannot be requalified without rotations.
Measured noise and magnetic timing are cleared because their gyro correction
basis has changed. Cancellation or fit failure leaves the saved record intact.
Only an already usable current record admits this tool; unsupported old records
still require full calibration.

When a console is connected, startup prints the current coefficients and frame
matrices. OU-III and TFG-CARMA3 also emit a best-effort `[HEADING LIVE]` record on
the transition to fused heading: `q` is fused heading, `m` is the direct magnetic
heading using the current tilt, and `d` is their wrapped difference in degrees.
This helps distinguish a Live heading discontinuity from an existing compass
error; it does not certify either heading against a physical reference.

## Gyro matrix and independent validation

The existing quiet-interval capture supplies bias and qualified endpoint means.
The initial set is twelve turn-and-hold segments: eight training rotations and
four independent validation rotations. Validation segments 3, 6, 9 and 12
are assigned before acquisition and never used to optimize the matrix.
`GyroCalibrationPlan` shares this policy with acquisition, magnetic registration
and record validation. The wizard selects an eight-training-segment profile;
the standalone gyro fitter retains its default twelve-training-segment setting.
Magnetic alignment uses the initial quiet pose plus eight training endpoints
(nine training poses), and the same four held-out endpoints. Conditioning,
coefficient plausibility, and train/validation residual limits are unchanged;
there is less acquisition redundancy, so physical qualification remains open. Timestamp gaps over 20 ms,
non-finite values, excessive rates, exhausted storage, and unqualified endpoints
discard the active segment. A retry establishes a fresh quiet start.

Before capture, three instruction cards wait indefinitely for a tap. They
explain support, twisting and tipping before any acquisition timer starts.
The per-turn start prompt is also untimed. Introductory cards are shown once,
not repeated after a failed turn or when extra turns are appended.

Each turn has four separate screens with one action at a time:

1. **HOLD — Start n/12**: steady it against the table. Supporting it with a
   hand is allowed; sensor motion, rather than hand contact, is measured.
2. **TWIST**: lift it, keep the same side up, and twist about a quarter-turn in
   the displayed direction. Left/right instructions alternate across turns.
3. **TIP**: tip to the named face: SCREEN UP, USB UP, LEFT DOWN, or RIGHT
   DOWN, then set it down. These use the device frame measured during the
   accelerometer stage. USB-down would rest on the cable and screen-down
   would hide feedback, so neither is requested.
4. **HOLD — Finish n/12**: keep the requested side steady against the table
   until **Turn kept**. The target is written once, followed by the support cue.

The short twist is necessary excitation: a replay of direct tips between
faces alone exposed a rank-deficient nine-parameter information matrix.
The twist cue advances after about 60 degrees of integrated rotation about
start gravity; it does not require a precise manually measured angle. A tip
lands within a generous 35-degree cone about the named face. The fitter still
uses measured endpoints and the original information, physical and residual limits. An
early tip before the twist is complete can establish a new quiet start for
that same turn rather than trap the operator in the original orientation.

The initial bias fit retains its six-second useful hold. Subsequent endpoints
use 2 seconds of useful samples plus magnetic settling, and at least 64
retained samples. A normal endpoint takes about 3.5 seconds and spans two
independent 1.5-second magnetic means. The initial flat-table bias hold retains its strict rate-scatter/change limits
(0.005/0.004 rad/s). Tilt endpoints instead use 0.08/0.04 rad/s: small braced
hand tremor can have appreciable rate scatter while preserving an accurate
mean gravity direction. Acceleration scatter/change, raw rate limits, magnetic
stability, freshness, final mean rate and independent fit validation still apply.
Physically coupled 0.1-degree rocking at 4, 6 and 8 Hz, including lever-arm
acceleration and measured-rate integration, passes the unchanged matrix and
independent-validation bounds. This is an endpoint profile, not a relaxation
of stationary bias or forty-minute noise acquisition. The visible bar includes settling and starts
with the first qualified block. Motion screens have no stillness bar.
Movement restarts only the current quiet hold, with a plain-language explanation.
The finish screen stays latched after entering the target cone; a transient
reading cannot silently send it back to TIP. The finish pose starts a fresh complete quiet window after entering its target.
A rejected final mean or misplaced endpoint restarts only that window on the
same target. Slow placement continues waiting; it neither qualifies a moving
endpoint nor opens the retry menu.
Every capture screen offers **Tap: retry turn**; this stops only the active
turn and opens its retry/bias-only/confirmed-exit menu, preserving accepted turns.

The gyro LCD uses the same pose rotations as the accelerometer guide: USB up
adds rotation 2, left edge down adds 1, right edge down adds 3, and screen up
uses 0. A dominant gravity direction must persist for 100 ms before changing
rotation. Confirmation, next-turn and retry menus retain this orientation.

A separate, 32768-byte capture task on core 0 reads the IMU at a 200 Hz target
cadence and owns the sensor-clock controller. The UI task reads atomic status
snapshots; LCD refreshes, button handling, and display wakeup cannot create
integration gaps. Actual delivery gaps over 20 ms still reject the active turn.
The target does not change during a valid turn; quiet waits, pauses, and every
unqualified/moving interval remain in the same controller. Each attempt is
bounded to three minutes. Periodic `[GYRROT]` diagnostics include phase, target,
quiet progress, restart count, and the sensor-side failure reason.

Retries discard only the active turn. Accepted turns, fitted bias, and completed
accelerometer/magnetometer stages remain. The retry screen names the current
turn and retained count. After twelve turns, fitting either qualifies the set
or offers four additional turns, keeping every prior capture. A further block
of four can be added if necessary, up to twenty total. Extra held-out flags are
also predetermined. A failed final fit never automatically erases the set:
one tap appends four turns when possible or retries the final fit; two taps
continue with measured gyro bias only. Exiting requires three taps and a
separate three-tap discard confirmation. There is no double-tap new-gyro-set shortcut.
The complete accelerometer and magnetometer stages are not rerun by this flow.
Captures are retained in this wizard session, not checkpointed across power loss. The 131072-sample PSRAM reservation (about 2.5 MiB)
includes ordinary pauses and slow motions across twenty turns and is checked
before gyro turns. Each turn can consume at most one twentieth of this buffer
(about 32 seconds of retained trajectory); a long failed attempt cannot use
the storage needed by later turns. Proven stationary waiting prefixes are
released. Fresh finite readings update sensor health before quiet-motion
gates, so a slow turn lasting over twelve seconds is not falsely classified
as a stalled IMU.

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

The initial twelve-turn set provides thirteen quiet magnetic endpoint poses.
The initial pose and eight training endpoints fit the existing constant-dip
alignment solver; four independent endpoints validate the frozen rotation
and dip. Appended turns keep their predetermined train/validation roles.
The undamped information must have full rank before the original numerical
ridge is added. The unified wizard reuses the existing iterative solve from
24 proper signed-axis seeds and allows full proper rotations. Each iteration
has a 0.35-radian step cap. The training poses select the best qualified dip
scatter; near-tied solutions separated by more than five degrees are rejected.
Independent poses never select a seed or adjust the answer. The fitted matrix
is always a proper rotation. Unobservable axes, training scatter over 1.5 degrees
and independent dip RMS over 1.5 degrees remain unacceptable. Existing diagnostic callers retain the
original six-iteration, small-trim profile unless they explicitly request full
registration. Driver axis conventions remain upstream and must establish consistent
right-handed frames: with unknown dip, constant-dip data alone cannot identify
a global field reversal caused by an incorrect reflected driver mapping.

Gyro fitting and magnetic registration are separate worker jobs. A magnetic
failure never reruns gyro fitting or requests more turns. Its diagnostic
identifies field-strength change, insufficient geometry, ambiguity, numerical
failure, dip inconsistency or failed independent validation. Serial ALIGN
records include counts, rotation, training scatter, independent RMS,
conditioning, time and stack headroom. On failure the wizard marks fine
registration skipped, retains magnetic intrinsics and all other fits, and
continues with nominal axes. Session data remains volatile across power loss.

### INSLIB comparison

INSLIB's [IMU calibration](https://github.com/jnz/INSLIB/blob/main/tools/inslib_imu_calib.py)
uses an initial stationary recording and varied stationary poses separated by
recorded motion. Its IMU-TK fit estimates accelerometer and gyro scale and
cross-axis terms. It does not require a named USB-up pose or a precision fixture.
Its generic CSV path supports other IMUs; better hardware does not change the
underlying endpoint-gravity principle.

Its [magnetic fitter](https://github.com/jnz/INSLIB/blob/main/tools/inslib_mag_calib.py)
first solves hard/soft iron, then calls
[align_by_dip](https://github.com/jnz/INSLIB/blob/main/tools/inslib_frame_align.py#L338)
with paired quiet gravity/field directions. This fits magnetometer-to-IMU PCB
registration by making dip constant across poses. Five poses are its minimum;
more varied poses improve observability. It warns at weak geometry, rotation
over 20 degrees and dip scatter over three degrees. These are warnings rather
than independent validation gates in that function. The common physical-field
assumption cannot distinguish a frame error from a changing magnetic field.

INSLIB makes magnetic calibration optional. In `solve_session()`, an exception
from magnetic fitting is logged and the successful IMU calibration is returned.
Housing alignment is separately optional: a level placement determines tilt,
and nonparallel known directions determine full rotation. It folds rotations
into output matrices. This project retains separately stored intrinsic,
mag-to-IMU, IMU-to-housing and housing-to-vessel transforms and supports qualified full/scale-only gyro fits, explicit bias-only operation,
and automatic nominal magnetic axes after registration failure. Fresh measured
bias and intrinsic sets are required in either case. Failed registration displays
a skip notice and continues with the successful IMU work.

Mounting declares how the device will be installed, not its current pose during
calibration. MOUNT explains the default: **screen up, USB points backward**.
One tap uses those device axes; two taps change the installation; three taps
open untimed help. No serial connection is required.

The button editor first names the case face that will point up, then a
perpendicular case face that will point toward the front of the boat. A tap
cycles faces, two taps choose, and three taps go back. The confirmation page
lists both choices and waits without a timeout. Left/right refer to the screen
viewed upright with USB at the bottom, as explained before editing. All 24
proper right-angle installations are supported; parallel directions and
reflections cannot be selected. The housing is the nominal AtomS3R case
(`R_HI=I`); the selected case-to-vessel rotation is stored in `R_VH`.

These declarations do not measure fine PCB-to-case tilt or arbitrary installation
angles. The retained Eigen SVD solves Wahba's problem in `solveFrameAlignment()`
for measured fixture directions and rejects collinear references: gravity alone
cannot determine unknown mounting yaw. The data model retains separate measured
`R_HI` and `R_VH` fields. Rotations are validated separately and composed only by
`RuntimeCals`.
A full intrinsic update requires renewed gyro/magnetic registration or an
explicit choice of the basic profiles; the
retired accelerometer-only UI could otherwise retain dependent qualifications.
The retained library candidate helper invalidates those dependent qualifications.

## Separate stationary noise workflow

After a usable current-schema calibration exists, **CALIBRATE → two taps:
saved fits → two taps: noise 40m** offers a separate 40-minute stationary capture.
Warm up the sensor and use a quiet,
thermally stable support. Checks cover fresh timing, frozen words, finite
data, acceleration/rate limits, fixed one-second pose means, fixed magnetic
means including loss/staleness, and temperature span below 2 C. At 200 Hz,
480001 samples span 2400 seconds. Missing frames invalidate the run;
samples are never concatenated across gaps. The worker also applies a 20 ms
freshness deadline when samples stop completely; the analyzer cannot detect
that case without the worker's independent clock.

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
leaving H P H' untouched. A gross innovation or field-norm anomaly is rejected, including when a large
measurement covariance would make its NIS small. All six cores supply the
field-strength ratio to the same gate; TFG reconstructs the world measurement
as residual plus reference.
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
Full recalibration supplies no previous temperature slope or stationary gyro
reference to the accelerometer capture/fitter. An unidentifiable thermal slope
remains explicitly unlearned. The retained numerical library can still validate
externally supplied thermal evidence; the wizard does not use it to substitute
previous coefficients. An existing current-schema record is used for the
separate noise workflow and verified write recovery, rather than to fill a
failed full fit.

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

The optional rotation buffer is reserved in PSRAM before gyro tilt capture
and released after fitting; it never occupies the loop-task stack. Bias-only
operation skips that allocation. Advanced tilt needs **Tools > PSRAM > OPI PSRAM**, followed by a
rebuild/upload (`PSRAM=opi` in the FQBN). The Arduino AtomS3 board defaults to
PSRAM disabled. The wizard reports the configured mode, detected/registered
external heap, total/free bytes, largest block and required allocation. Missing
PSRAM produces a setup instruction and a bias-only choice before gyro turns; low or fragmented
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

