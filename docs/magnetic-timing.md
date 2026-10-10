# Magnetic timing and relative rotation

## Acquisition audit

The AtomS3R build pins M5Unified 0.2.13 (`a6256725481f1bc366655fa48cf03b6095e30ad1`).
Its `BMI270_Class::begin` sets BMM150 normal mode to 30 Hz and leaves BMI270
`AUX_CONF` at reset `0x46` (25 Hz). The project preserves those settings and
sets 47 XY / 41 Z repetitions. Bosch's conversion-time formula gives
`145*47 + 500*41 + 980 = 28295 us`. This is a conversion aperture, **not an
instantaneous vector or a known group delay**. Axis effective times need not
coincide. Neither the acquisition API nor the board code exports a conversion
timestamp. There is no wired/serviced BMM150 DRDY interrupt in this path.

BMM150 has a DRDY pin and register `0x48` bit 0; reading data over AUX clears
the sensor flag. BMI270 has separate AUX-ready flags; the pinned driver reads
`INT_STATUS_1`, then a 20-byte AUX/accel/gyro burst. It updates its magnetic
cache on AUX-ready. A second eight-byte AUX read supplies RHALL for Bosch
compensation. These two reads can straddle an AUX refresh. A change in the
cached XYZ detects a new candidate, not a unique conversion identifier.
Identical quantized conversions are indistinguishable in the existing API.

`sample_us` is host **read-start**, before `M5.Imu.update()`, not sensor sample
time. `readImuMapped` then reads temperature and maps/calibrates the sensors.
The magnetic source repeats its last successful reading, for at most 200 ms.
The OU-III sketch's 35 ms spacing gate can accept a held sample repeatedly;
it does not track source identity. The filter performs gyro prediction,
scheduled S/accelerometer work, then magnetic correction at the latest
attitude. No magnetic rotation transport or latency covariance was present.

The private Mahony observer is accel/gyro-only. Initial magnetic acquisition
pairs cached magnetism with current yaw-removed attitude; refinement and the
optional continuous hard-iron tracker use current proxy tilt. Their magnetic
clocks are accumulated filter time at delivery, not conversion time. The
shipping OU-III sketch disables continuous hard iron after refinement; the
facade default and simulations still support it. NVS calibration already
correctly applies `A*(raw-b)`, Bosch trim/RHALL compensation is already present,
and the MEKF already uses the correct world-to-body sign and counts actually
applied corrections as service. None of these is a synchronization mechanism.

Physical acquisition-to-update latency and jitter require hardware capture.
At normal polling, conversion completion phase alone spans a 33.3 ms sensor
period; AUX has a separate 40 ms period. Host cache age and read-to-correction
time can be instrumented. A software timestamp must never be advertised as
the BMM150 measurement epoch. No measured AtomS3R timing is claimed here.

Primary sources: [pinned M5 driver](https://github.com/m5stack/M5Unified/blob/a6256725481f1bc366655fa48cf03b6095e30ad1/src/utility/imu/BMI270_Class.cpp),
[BMM150 datasheet, sections 4.2–4.3 and 4.6](https://www.bosch-sensortec.com/media/boschsensortec/downloads/datasheets/bst-bmm150-ds001.pdf),
[BMI270 datasheet, AUX interface and AUX_CONF](https://www.bosch-sensortec.com/media/boschsensortec/downloads/datasheets/bst-bmi270-ds000.pdf).

## Relative rotation and uncertainty

Shipping `qref` maps world to unheeled body B'; prediction left-multiplies
`Exp(-(gyro-bias)*dt)`. In physical body axes, integrating these increments
chronologically gives `D = R_bk,bm = R_wb(tk) R_wb(tm)^T`. Thus
`m_k = D (m_cal(tm)-b_HI)`, never `D m_cal - b_HI`. NVS soft iron and its
body-fixed offset must also be applied before transport. De-heeling occurs
inside the MEKF; transport is in physical sensor/body axes.

`RotationHistory` stores fixed-capacity time intervals, bias-corrected rates
and their negative-angle quaternion increments. It uses the same end-sample
piecewise-constant rate convention as shipping prediction, composes rotations
rather than adding rotation vectors, and subdivides boundary intervals. It
does not estimate attitude or take magnetic/accel corrections. There is no
extrapolation. A gap invalidates the chain; a clock wrap is safe for intervals
less than half the uint32 range. History and query latency are separately bounded.

For a transported vector `v`, the first-order rotation Jacobian has covariance
`[v]x Q_theta [v]x^T`. An orientation-independent upper bound is

`Q_theta <= I (q_g*T + lambda_b*T^2 + (|omega_m| sigma_tm + |omega_k| sigma_tk)^2)`.

Here `q_g` is the largest gyro noise-density variance, `lambda_b` bounds the
residual bias covariance eigenvalue, and timing errors may be fully correlated.
The bias term is quadratic in time because the same unknown bias persists;
it is not re-drawn at each IMU sample. The endpoint timing bound uses Cauchy,
not an independence assumption. This linearization requires small angular
uncertainty; finite aperture, angular acceleration, miscalibration and field
gradients are separate model errors, not automatically white noise.

For independent noises, `R_aligned = D R_m D^T + R_transport` and
`Sigma_r = R_k + D R_previous D^T + R_transport`. The live gyro-bias estimate
can depend on earlier magnetic samples. With unknown cross covariance the
library uses the conservative inequalities `Cov(x+y)<=2(Cov(x)+Cov(y))`
and `Cov(x+y+z)<=3(Cov(x)+Cov(y)+Cov(z))`. These do not make the transported
measurement independent of the MEKF state. A current-time update still
approximates the correlated update; timestamp-correct rewind/replay is the
validation reference, not a second deployed filter.

Adjacent residuals reuse magnetometer samples and sometimes gyro intervals;
they are not independent trials. `d^2=r^T Sigma_r^-1 r` uses the full 3x3
covariance. The chi-square(3) 99.9% value 16.266236 is an advisory reference,
not a claimed hardware false-alarm rate under bounded/model errors. No new
disturbance gate or covariance adaptation is driven by this statistic.

Normal rigid rotation should have a small residual. Timing errors tend to
scale with rate; integration/bias errors with interval and excitation;
hard/soft-iron errors with orientation; interference can change norm/dip or
produce abrupt residuals. These are clues, not uniquely identifiable causes.
A small residual during weak motion cannot validate heading, scale, a constant
world-field error or calibration. The diagnostic reports rotational
signal/noise rather than declaring unexcited data to be proof of correctness.

## Stability-study boundary

The 21-state OU-III, OU dynamics, integrated chain, SpectralMSE/R_S adaptation,
bias dynamics, magnetic Jacobian and quaternion reset are authoritative.
Timing preprocessing changes the delivered magnetic forcing and its effective
covariance, not the state architecture. MAGNETIC SERVICE still means actually
applied informative updates; candidates, duplicate packets and diagnostic
success do not establish it. Hardware conversion delay, proxy timing error,
calibration residual and shared gyro/state correlation remain proof-side
dependencies. This work does not close the stability theorem or its open
receipt/Schur, capture, physical qualification or all-time service obligations.

## Deployment stages and remaining synchronization limit

The first implementation commit is observation-only (Stage A). In the final
OU-III sketch, `SEA_STATE_MAG_HOST_ALIGNMENT=0` selects that behavior for an
A/B capture. The default Stage B transports only the **known host-cache age**
from the first observation frame to the present IMU frame. It does not subtract
an invented 14 ms or 33 ms sensor delay. Consequently physical BMM150 timing
mismatch is **not fully corrected**. The same transport API accepts externally
qualified measurement times, as exercised by the native replay tests.

The source keeps a read bracket and a monotonically increasing host observation
sequence. The sketch retains the 35 ms spacing limit, submits an observation
at most once, and rejects only invalid input or unavailable/too-old/gapped
rotation history (120 ms cap, 20 ms maximum gyro interval). Equal quantized
conversions and the two-read AUX race remain source-identity limitations;
these sequences must not be described as hardware conversion IDs. Typically
25 Hz source observations pass once each. The previous artificial repeated
corrections during a stall disappear, so the previous service count cannot be
carried over without auditing actually applied events.

At first observation, the frontend saves its actual private proxy quaternion,
attitude, accelerometer and gyro. Initial reference/refinement and optional
continuous hard-iron accumulation use those snapshots and observation time.
The proxy is never reconstructed using the MEKF bias estimate. After the
existing hard-iron/reference operations, the applied body offset is subtracted,
then the final MEKF input and its covariance are transported. No magnetic
reference gate, startup window, release count or continuous-estimator setting
is changed. The sketch still disables continuous hard iron; the identity
transport regression also covers the facade's enabled default.

For covariance transport, the history retains the largest bias uncertainty
encountered, rather than assuming that the smaller current posterior bounded
all past increments. With unknown correlations inside the rotation error its
standard deviations are added before squaring. The factor-two/factor-three
bounds above then cover correlation with the magnetic noises. This remains a
first-order statistical model; it is not a proof of the physical bias limits.
The live transport uses the measured host read bracket; unobserved conversion
phase/aperture is separate. The observation-only consistency score uses a
conservative 28.295 + 33.333 + 5 ms endpoint allowance for normal host polling,
not a claimed measured distribution or a bound through arbitrary host stalls.

Stage C is deliberately not enabled. Under constant angular velocity a fixed
magnetic delay can give a small consecutive rotation residual despite a biased
heading. Constant scale error, a constant world-field change and quiet motion
also need independent checks. High residuals identify unexplained change,
not its unique cause. Neither synthetic statistics nor the advisory chi-square
reference justifies removing real-vessel magnetic service.

The diagnostic serial row is emitted only by the existing low-rate serial
stream, never at gyro rate. `host_age_ms` and `read_to_update_ms` are host
measurements; `jitter_range_ms` is the running min/max range of the latter,
not physical sensor jitter. Rotation angles are in radians, residual in uT,
`d2` is dimensionless, and `applied` is the last attempted correction's actual
result. `duplicates` counts held polls, not rejected physical conversions.

Hardware qualification still needs logic-analyzer/DRDY evidence of conversion,
AUX and host read phases, per-axis aperture and timestamp uncertainty, plus
AtomS3R loop/stack measurements under display and USB load. No measured hardware
latency or 200 Hz timing guarantee is inferred from host benchmarks.

## Deterministic validation and results

`tests/kalman_ou_iii/mag_timing-test.cpp` drives the actual 21-state core with
one immutable motion/noise record and three paths: legacy delayed input,
current-time transport, and test-only rewind/replay. The reference restores
the pre-sample state/covariance and repeats every gyro/S/acc/magnetic operation
at its original epoch; it never consumes future samples. It is not deployed.
The test keeps the shipping gyro density (0.00135 rad/sqrt(s)) and magnetic
standard deviation (0.8 uT). Synthetic magnetic noise is 0.15 uT per axis; the
low resulting NIS is expected, not grounds for retuning shipping constants.

The 26 s fixtures score heading after 3 s. All paths receive identical
calibrated samples, the same known body offset, and identical noise. Magnetic
vector RMS compares the delivered vector to the correct field at that update
epoch (the replay row uses its historical epoch). Peak error and rejection
rate use the definitions in the test. The quiet and biased cases deliberately
have no improvement assertion. Results below are native GCC 13.3, `-O1`,
IEEE finite checks; they are **not measured vessel accuracy**. Full precision,
all delays, peaks, update rejection and post-disturbance tail results are in
[paired-replay.csv](../reports/results/magnetic_timing/paired-replay.csv).

| Scenario / delay | Legacy heading RMS, deg | Aligned RMS, deg | Replay RMS, deg | Aligned vs replay RMS, deg |
|---|---:|---:|---:|---:|
| yaw / 5 ms | 0.334 | 0.211 | 0.229 | 0.023 |
| yaw / 10 ms | 0.520 | 0.211 | 0.229 | 0.023 |
| yaw / 20 ms | 0.943 | 0.211 | 0.229 | 0.023 |
| yaw / 40 ms | 1.825 | 0.211 | 0.230 | 0.024 |
| combined variable / 40 ms | 1.070 | 0.219 | 0.201 | 0.043 |
| async jitter / 20 ms | 0.702 | 0.213 | 0.190 | 0.041 |
| gyro bias error / 40 ms | 0.828 | 0.799 | 0.690 | 0.131 |
| soft iron residual / 20 ms | 1.643 | 1.023 | 1.301 | 0.363 |
| interference recovery / 20 ms | 15.369 | 11.201 | 14.870 | 4.265 |
| quiet / 40 ms | 0.203 | 0.197 | 0.204 | 0.018 |

At 40 ms yaw delay, vector RMS is 0.693 / 0.254 / 0.253 uT and peak heading
error is 2.278 / 0.416 / 0.446 degrees (legacy / aligned / replay). No valid
magnetic update was suppressed in any paired fixture. Interference is a
20,-15,5 uT body-field step from 10 to 15 s. After 22 s its heading RMS is
2.096 / 0.751 / 2.953 degrees. Transport differs from replay by 4.265 degrees
RMS across the disturbed run: it is **not equivalent to a correlated delayed
Kalman update**, and the different covariance weighting affects recovery.
There is no claim that this covariance bound or recovery result generalizes
to hardware interference. The clean-motion comparison tolerance is 0.3 deg;
measured differences are below 0.044 deg.

Coverage includes pure yaw; combined roll/pitch/yaw with variable rates; all
5/10/20/40 ms delays; asynchronous 4–6 ms gyro and 33–37 ms magnetic schedules
with 0–5 ms delivery jitter; gyro-bias error; subtraction before rotation;
fixed soft-iron error and corrected calibration; sudden interference/removal;
quiet/no-excitation ambiguity; duplicate, missing, stale, future, out-of-order,
invalid/nonfinite samples; partial intervals; and uint32 clock wrap.
`mag_rotation-test.cpp` disables Eigen allocation throughout and advances two
million gyro samples (10,000 s), spanning two micros wraps. The independent
white-noise marginal check gives mean d2=3.034858 for 4,999 residuals, with four
above the chi-square advisory reference; adjacent values are not iid. A
constant-delay uniform-yaw case explicitly demonstrates the diagnostic's
timing blind spot. A 130 s facade comparison covers startup, reference
refinement and the enabled continuous-hard-iron default at identity transport.
Existing calibration/continuous-hard-iron and shipping-contract regressions
are retained. No estimator tuning constant is changed to meet these checks.

## Embedded resource assessment

The history is fixed at 32 intervals: 1,552 bytes with the native ABI, versus
zero history bytes on main. The observation snapshot is 176 bytes, previous
magnetic sample/covariance 60 bytes and diagnostic row 28 bytes. With metadata,
counters and configuration the added persistent storage is approximately
2 KiB; actual Xtensa alignment/linker totals must come from the board build.
No heap allocation occurs in history, transport or consistency calculations.
The source's extra sequence/read-bracket fields add 16 bytes per metadata
copy. No second 21-state estimator or replay buffer is shipped.

Every gyro sample adds one 3-vector bias subtraction, fixed-size bias-frame
access/covariance bound, one norm, sine/cosine, normalized quaternion and ring
write. Lookup inspects at most 32 intervals, multiplies at most 32 quaternions
and evaluates at most two partial-interval exponentials. It never scans an
unbounded queue. Magnetic work adds fixed 3x3 covariance transforms/LLT solves
and vector/quaternion operations; trigonometric diagnostic angles run only
on a new host observation. Serial output uses the existing low-rate stream.

The isolated `-O3` x86-64 benchmark recorded 20.5 ns/push, 206.7 ns/maximum-age
lookup, 115.5 ns/alignment and 295.8 ns/consistency-with-push. These are native
constant-rate microbenchmarks, **not MCU bounds or whole-loop timings**; the
raw output is in [native-resources.txt](../reports/results/magnetic_timing/native-resources.txt).
GCC `-fstack-usage` reports individual static frames of 48 bytes (increment),
160–208 (push), 160 (lookup), 240 (alignment) and 560 (consistency), plus called
helpers and the existing estimator frames. They are not a bound on the whole
AtomS3R call stack. The design keeps constant per-gyro work and bounded
per-magnetic work within the 5 ms schedule; actual loop worst-case time,
flash/RAM change and task-stack high-water marks remain hardware/board-build
qualification items. No physical 200 Hz guarantee is claimed.

## Shipping-faithful handoff

1. Preserved: all 21 states, literal prediction/S/acc/mag/reset ordering,
   OU chain, coupled tuner/SpectralMSE/R_S, bias dynamics, Jacobian and
   reference acquisition/release conditions. Only magnetic input preparation
   and host-observation identity change.
2. Relaxations: first-order rotation uncertainty, conservative unknown-cross-
   covariance bounds and current-time correction. Synthetic point timestamps
   omit the BMM150 finite aperture. These tests do not admit or certify the
   physical shipping motion/service family.
3. Failures: implementation/provenance issues are class E in the research
   ledger. Large-interference replay disagreement identifies the approximation's
   limit; it is not a theorem counterexample or a source-uniform bound.
4. No genuine shipping stability counterexample was found.
5. Retained results: original legacy operations, existing algebraic lemmas and
   finite controls remain valid in their stated scope. The new timed-input
   gyro-bias/covariance reverse ports are not covered by LR1's legacy zero-port
   result. No theorem status is promoted.
6. Next: qualify physical conversion/AUX/read timing and actually applied
   service on the same hardware history, then carry transport's bias/covariance
   derivatives and shared-noise terms into the existing finite-error supply
   inequality. Do not substitute a timestamp guess or an independent estimator.
