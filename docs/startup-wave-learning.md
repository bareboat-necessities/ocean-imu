# Gravity and magnetic-field startup in waves

The experiment in [`tests/startup_learning`](../tests/startup_learning/README.md)
compares the previous Mahony-proxy handover with five full-history reference
fits and the fixed-memory implementation used by the live wave filters.
The base is `8c9f535709fe7ff4d87d75dc5c98f4d0a48daa7b`; both cases use the pinned
`oceanography-waves-lib` v1.4.0 archives, at 5 kn with ±3° / 20 s heading motion.

## Production handover

The bounded-memory solver seeds attitude, the local magnetic field and gyro
bias in the live wave filters. The table compares actual OU-III handovers with
the learner disabled and enabled. Each cohort contains 120 starts per
direction and sensor profile. Both policies use the same eligibility checks.
The figures below use calibrated magnetometers.

| Cohort | Seas | Handovers by 35 s | Previous → deployed tilt RMS | Previous → deployed yaw RMS |
|---|---|---:|---:|---:|
| Validation | Following | 109/120 | 0.905° → 0.333° | 2.202° → 0.779° |
| Validation | Head | 120/120 | 1.157° → 0.419° | 3.715° → 0.967° |
| Confirmation | Following | 95/120 | 0.598° → 0.362° | 1.916° → 0.665° |
| Confirmation | Head | 120/120 | 0.971° → 0.466° | 2.521° → 0.924° |

With residual magnetic calibration, following-sea yaw changed from 2.494° to
2.449° in validation and from 3.633° to 2.550° in confirmation. Head-sea yaw
changed from 4.680° to 1.661° and from 4.106° to 3.305°. All 888 matched
handovers across 960 starts used the joint estimate. Each policy missed the
same 72 deadlines. The full-history experiment below also scores the
35-second endpoint on every trial. Its optional hard-iron fit is experimental.

Run `tests/startup_learning/startup-learning --production-pair` once per input
with `--input`, `--output`, `--profiles calibrated,residual`, and the
`--seeds` and `--offsets` in the
[production scores](../reports/results/startup_learning/production-handover.json).

## Accuracy after the startup transient

Both v1.4.0 directions were rerun on the same 13 binaries with the pinned GCC
toolchain. The 544 filter/profile/record cells include 2,768 gated
measurements. Relative to the pre-integration locked measurements, 1,127 fell,
606 were equal and 1,035 rose. With the existing 0.5% gate margin as the
comparison threshold, 55 fell by more than 0.5% and 44 rose by more than 0.5%.
The largest increase was TFG pitch on the following 8.5 m JONSWAP record
(0.08603° to 0.09181°). For the eight production wave-filter families, gyro
bias RMS decreased in every family when pooled over the 16 records. Yaw RMS
decreased 0.044% for OU-III and rose 0.009–0.054% in the other seven families.
The ordinary wave-filter gate scores the final 900 s of each 1,200 s run.

[Paired RMS by family, profile and sea direction](../reports/results/startup_learning/live-replay.csv)
· [Source measurements and CI artifacts](https://github.com/bareboat-necessities/ocean-imu/actions/runs/44184184023)

## Results within 35 seconds

The validation set has five unused sensor-noise/calibration seeds, three later
power-on phases, eight wave records and both travel directions: 480 startups
across two sensor profiles. Every method receives identical noisy samples.
Only sensor synthesis and scoring can read truth. This is held-out replay on
the same wave files, not new independent ocean realizations.

| Equal physical handover epoch | Shipping yaw RMS | Candidate yaw RMS | Shipping tilt RMS | Candidate tilt RMS | Candidate |
|---|---:|---:|---:|---:|---|
| Following, calibrated magnetometer | 2.202° | 0.802° | 0.905° | 0.333° | joint |
| Head, calibrated magnetometer | 3.715° | 0.966° | 1.157° | 0.408° | joint |
| Following, residual magnetic calibration | 2.494° | 1.876° | 0.905° | 0.370° | joint_cal |
| Head, residual magnetic calibration | 4.680° | 1.640° | 1.157° | 0.437° | joint_cal |

Tilt is the angle between the estimated and true unit gravity directions in
the current body frame. It measures gravity orientation independently of yaw;
it is not a separate roll or pitch Euler-angle RMS.

All 240 calibrated-magnetometer trials meet yaw ≤2° and tilt ≤1° at 35 s.
The worst calibrated endpoint errors are 1.387° yaw / 0.574° tilt following
and 1.545° yaw / 0.610° tilt head. Independent fixed-deadline accuracy gates
include every trial, even when shipping never hands over; paired gates alone
would omit those deadline misses. The reports retain both RMS and worst errors.
The shipping wrapper reaches Live by 35 s in 458/480 trials. The 22 deadline
misses remain in the coverage file and the fixed-horizon results; they are
excluded only from the comparison requiring an actual shipping handover.
The baseline is stopped before any live magnetic correction, so these are
startup seed errors, **not subsequent live-navigation RMS**.

For calibrated magnetometers, matched-time tilt RMS falls from 0.905° to
0.333° following and from 1.157° to 0.408° head; gyro-bias vector RMS falls
from 0.0537 to 0.0244°/s and 0.0535 to 0.0251°/s respectively.

Full per-method tables, readiness failures, latency, field magnitude/direction
and replay provenance are in the
[validation report](../reports/results/startup_learning/validation/REPORT.md).
The [ideal/field-step controls](../reports/results/startup_learning/controls/REPORT.md)
keep adversarial failures visible, including numerically plausible fits that
must not be handed over as ready.

![Startup convergence](../reports/results/startup_learning/validation/convergence.svg)

## Why the experiment helps

The existing provisional magnetic acquisition averages in a yaw-removed
Mahony tilt frame. Its stored absolute yaw is the acquisition-window gauge,
which is then written at the later handover epoch. With heading oscillations,
the window's mean heading need not equal the current heading. A provisional
tilt error also changes the learned magnetic reference; current shipping
defaults start another reference refinement at 90 s.

The isolated candidates instead keep one gyro-transported initial frame from
the first sample. They learn a constant, unrestricted three-component local
magnetic vector there, without a gravity/attitude guess. Gravity is represented
by a **unit direction**; its known magnitude is never fitted. The tapered force
statistic from `doc/kalman_ou_iii/ins-startup.tex` suppresses finite-window wave
contamination. Finally, the magnetic-NED basis is composed with the terminal
relative rotation, producing a current-epoch attitude.

The `joint` fit combines magnetic constancy, a gyro-calibration prior and a
wave-tolerant constraint on the long-window transported force slope. The
additional force constraint helps the gyro-bias direction that one magnetic
vector observes poorly. The `joint_cal` ablation also fits a regularized
constant hard-iron offset. Every candidate bias change repropagates all prior
gyro samples; terminal-only bias correction would be inconsistent.

Fitting hard iron helps the residual-calibration profile, but worsens the
already-calibrated profile. The code does not choose a method using truth or
the dataset label. Both are explicit experiment arms; a future device policy
would need externally justified calibration uncertainty/quality.

## Force-slope weighting in low waves

The fixed slope uncertainty can give a noisy gyro-bias estimate too much
influence on terminal tilt in low waves. The experimental
`--scale-force-slope` arm measures the time-weighted RMS fluctuation of
specific force after gyro transport into the initial frame. It scales the
existing slope standard deviation by `clamp(RMS / (1 m/s²), 0.15, 1)`.
Lower force variation increases the weight of the slope observation.
The scale uses the calibration-prior gyro trajectory and stays fixed through
each solve. The reference amplitude and floor are experimental settings.

On the original 480 starts, the number of direction/profile/record groups
with lower matched-time tilt RMS than shipping increases from 27/32 to
30/32. The two remaining regressions occur in 0.270 m head seas with
residual magnetic calibration. The largest increases from shipping are
0.0851° with fixed weighting and 0.0065° with scaling.

On a further 480 starts using seeds 31–35 and offsets 41, 191 and 421 s,
scaling improves tilt over shipping in all 32 groups. These seeds and phases
reuse the same 16 wave records. Pooled calibrated errors at the actual
shipping handover are:

| Cohort | Seas | Matched starts | Fixed → scaled tilt RMS | Fixed → scaled yaw RMS |
|---|---|---:|---:|---:|
| Original | Following | 109 | 0.333° → 0.327° | 0.802° → 0.793° |
| Original | Head | 120 | 0.408° → 0.393° | 0.966° → 0.932° |
| Further starts | Following | 95 | 0.374° → 0.380° | 0.662° → 0.675° |
| Further starts | Head | 120 | 0.489° → 0.471° | 0.944° → 0.897° |

All 480 calibrated starts across both cohorts meet the 2° yaw and 1° tilt
targets with scaling at 35 s, including starts that shipping misses.
Both weightings reject the 32 field-step controls. With residual magnetic
calibration, scaled yaw RMS on the original cohort increases from 1.876°
to 1.892° following and from 1.640° to 1.678° head-on. The mixed results
support keeping the fixed weighting as the default and comparing the two
fits explicitly.

[Original per-record results](../reports/results/startup_learning/slope_study/validation.csv)
and [confirmation results](../reports/results/startup_learning/slope_study/confirmation.csv)
give matched and 35 s RMS for tilt, yaw and gyro bias. Each handoff row is
the shipping baseline shared by both arms.

## Handover contract

| Field | Meaning |
|---|---|
| `time`, `attitude` | Acquisition epoch and body→magnetic-NED quaternion at that epoch |
| `down_body` | Unit gravity direction in current body axes; multiply by known gravity magnitude |
| `field` | Learned magnetic-NED vector in µT, retaining magnitude and inclination |
| `gyro_bias`, `accel_bias`, `mag_bias` | Body-frame residual calibration seeds; acceleration bias remains its prior |
| `angular_rate` | Latest instantaneous body rate minus fitted gyro bias, rad/s |
| `last_raw_gyro` | Endpoint sample to prime the next interval without replaying the last one |
| `gyro_bias_cov` | Conditional gyro-bias covariance; includes fitted hard-iron nuisance elimination when enabled |
| `valid`, `ready` | Numerically usable estimate versus a model-consistency readiness heuristic |

Gravity direction plus the learned magnetic vector defines magnetic north.
True north requires external declination. Position, velocity and wave states
need their own initialization and uncertainty; this attitude experiment does
not claim to recover them. A full live handover also needs attitude/reference/
bias/wave cross-covariances, calibration identity and scheduler/history state.
The current conditional covariance is explicitly insufficient for that step.

## Boundaries and reproducibility

Residual magnetic calibration is still a material limitation: good aggregate
RMS does not mean every handover meets the 2° / 1° target. A constant-field
optimizer can also explain a magnetic step using a large invented bias, so a
fit beyond four prior standard deviations cannot be marked ready. Low residual
alone is not a certificate. Neither the bounded-displacement statistic nor the
force-slope prior proves that physical acceleration has zero mean in 35 s.

The hardware-neutral learner validates sample epochs, gaps, finite values,
magnetic freshness and the 35 s maximum window. Actual magnetic acquisition
timestamps must be supplied; unknown delivery latency is not solved here.
Recorded host solve time does not qualify ESP32 memory or execution time.

Run the commands in the [test README](../tests/startup_learning/README.md).
The focused `isolated startup learning` workflow reproduces the complete
held-out paired matrix. Unit tests cover frames, sample timing, gravity/field
recovery, gyro-bias recovery, quiet water, wave contamination, invalid packets,
deadline enforcement, stale/vertical magnetic fields and field-step rejection.
The experiment is in both Make and CMake host builds. The production-pair run
exercises the actual shared live wrapper with both startup policies.
