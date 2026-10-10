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
