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

## Which handover values help

The previous handover already supplied proxy tilt, an acquired yaw, a local
field and the last raw gyro sample. The shared learner supplies tilt and yaw
at the handover epoch, a fitted magnetic reference and a body gyro-bias seed.
Across the 888 matched handovers, corrected angular-rate RMS fell from 0.1439
to 0.1311°/s in following seas and from 0.1227 to 0.1157°/s in head seas.
The next prediction uses the existing raw gyro endpoint and the fitted bias.

We restored one previous value at a time at the actual OU-III seed assignment.
Each row below pools the first second of live output over the same 32 cold
starts (16 records, two sensor-noise seeds). All seven policies received the
same samples. The joint fit qualified on 29 starts. The three fallback starts
contribute to every score.

| Handover policy | Tilt RMS (°) | Yaw RMS (°) | Gyro-bias RMS (rad/s) | Accel-bias RMS (m/s²) |
|---|---:|---:|---:|---:|
| Previous handover | 0.8439 | 4.8111 | 0.000782 | 0.07859 |
| Deployed handover | 0.5436 | 3.7367 | 0.000423 | 0.07859 |
| Restore previous tilt | 0.7148 | 3.8927 | 0.000423 | 0.07859 |
| Restore previous yaw | 0.6044 | 4.5959 | 0.000423 | 0.07859 |
| Restore previous magnetic field | 0.5496 | 3.7350 | 0.000423 | 0.07859 |
| Restore previous gyro bias | 0.5444 | 3.7374 | 0.000781 | 0.07859 |
| Add experimental gravity-axis accel-bias seed | 0.5459 | 3.7363 | 0.000423 | 0.06527 |

The learned yaw accounts for most of the pooled initial heading gain. Head-sea
yaw fell from 5.5992° to 3.3475°. Following-sea yaw rose from 3.8656° to
4.0891° in this smaller fixed-start cohort. Across both directions the
10–60 s window gave 0.4576° to 0.4522° tilt and 3.7332° to 3.6685° yaw.

The experimental accelerometer-bias seed uses the mean force along gravity to
estimate one bias component and keeps the calibration prior in the two
orthogonal directions. It reduced first-second bias RMS by 17%. In a second
comparison, 48 starts (16 records, three seeds) ran for the full 1,200 s under
each of three policies. The fit qualified on 44 starts.

| Final 900 s, pooled RMS | Previous | Deployed | With gravity-axis accel-bias seed |
|---|---:|---:|---:|
| Yaw (°) | 0.81778 | 0.81605 | 0.81776 |
| Accelerometer bias (m/s²) | 0.039893 | 0.039884 | 0.039833 |
| Y displacement (m) | 0.246066 | 0.246011 | 0.246079 |
| Heave (m) | 0.197872 | 0.197870 | 0.197949 |

The four restored-component policies were also run for all 48 full-length
starts, completing 336 long replays across the seven policies. Late Y RMS
was 0.246011 m with the complete handover, 0.246016 m with previous tilt,
0.246048 m with previous yaw, 0.246010 m with the previous field and
0.246022 m with the previous gyro bias. Late bias RMS across those five
policies was 0.039882–0.039896 m/s². The component gains are predominantly
transient; no component supplies a large sustained Y or bias improvement.

The extra bias seed reduced late accelerometer-bias RMS by 0.13% relative to
the deployed handover, while following-sea yaw rose from 0.80318° to 0.80812°.
The production filters therefore retain the calibrated accelerometer-bias prior.
The late scores use the standard validation window. Each pooled RMS is the
square root of the sample-weighted mean squared error.

[Per-start results and production-control checks](../reports/results/startup_learning/handover_components)
· [Reproduction commands](../tests/startup_learning/README.md#reproducing-the-handover-component-comparison)

## Accuracy after the startup transient

The same 13 production executables scored all 544 filter, profile and record combinations in head and following seas. Against the pre-integration baseline, 1,136 of 2,768 measured RMS channels fell, 620 were equal and 1,012 rose. Changes exceeded 0.5% on 74 improving and 38 worsening channels. These counts use measured RMS before the regression margin is applied. Wave-filter scores use the final 900 s of each 1,200 s record.

On the default 8.5 m following JONSWAP record, TFG pitch RMS was 0.08603° before startup integration, 0.09181° with the original 2° handover uncertainty and 0.07858° with the accepted-fit uncertainty below.

A second run using the checked-in limits passed all 544 cases. All 2,768 numerical values reproduced exactly. Both sea directions used the same 13 binaries. [Passing production matrix](https://github.com/bareboat-necessities/ocean-imu/actions/runs/44214928709) · [Archive and measurement hashes](../reports/results/startup_learning/tfg_tilt_uncertainty/verified-matrix.json)

[Paired results by family, profile and sea direction](../reports/results/startup_learning/tfg_tilt_uncertainty/compared_to_preintegration.csv) · [Production matrix and raw logs](https://github.com/bareboat-necessities/ocean-imu/actions/runs/44209357875)

## Tilt uncertainty after an accepted TFG fit

TFG uses a 0.5° standard deviation on each horizontal attitude-error axis when the wave fit supplies the attitude. The fallback Mahony proxy retains its own setting of roughly 2°. Both paths keep the configured yaw uncertainty. Over 458 validation handovers and 430 confirmation handovers with the fit accepted, radial tilt RMS was 0.384° and 0.446°, respectively. The corresponding maxima were 0.672° and 0.875°. These are measurements on the pinned wave records.

We compared the unmodified TFG executable at each of the two revisions on 16 wave records with the default sensor draw and seeds 7, 41 and 101. Sensor samples, coordinates and scoring windows match within each pair. The table pools all 64 cases using the Ubuntu system math library used by the pinned GitHub replay job. The windows start at fixed times after power-on. On the seed-7 low JONSWAP following-sea record, handover occurred at about 41.85 s, after the 40 s window began.

| Quantity | 40–60 s, 2° → 0.5° | 300–1,200 s, 2° → 0.5° | Late change |
|---|---:|---:|---:|
| Roll RMS (°) | 1.11949 → 1.10082 | 0.28084 → 0.27945 | -0.495% |
| Pitch RMS (°) | 2.62114 → 2.56580 | 0.28411 → 0.27995 | -1.463% |
| Yaw RMS (°) | 6.92899 → 6.81220 | 0.83171 → 0.82438 | -0.882% |
| X displacement RMS (m) | 4.962575 → 4.855786 | 0.301550 → 0.301660 | +0.037% |
| Y displacement RMS (m) | 2.117846 → 2.096079 | 0.245368 → 0.244848 | -0.212% |
| Z displacement RMS (m) | 2.224566 → 2.213637 | 0.145095 → 0.145081 | -0.010% |
| 3-D displacement RMS (m) | 5.836191 → 5.733445 | 0.414958 → 0.414726 | -0.056% |
| 3-D accelerometer-bias RMS (m/s²) | 0.066638 → 0.066638 | 0.068162 → 0.067475 | -1.008% |

Late pitch, Y displacement and accelerometer bias improved on each of the four pooled draws. The default 8.5 m following JONSWAP case moved from 0.09181° to 0.07858° late pitch RMS. Its late yaw increased from 0.43950° to 0.45946°. The pooled late X increase was 0.000110 m. The production change applies to TFG-OU. The OU-III experiment with a tighter prior had mixed results.

We tested whether reducing accepted-fit yaw uncertainty from about 5° to 3.5° would recover the JONSWAP yaw loss. Its late yaw changed from 0.45946° to 0.45963°. Across all 64 matched cases late yaw rose by 0.018% and late pitch rose by 0.085%, so production keeps the existing yaw setting. [Yaw-prior sensitivity scores](../reports/results/startup_learning/tfg_tilt_uncertainty/ubuntu/fitted_yaw_comparison.csv.gz)

The production replay uses its existing 0.5% margin for every record. The limits below compare the previous PR gates with the limits measured at `68f0978`. The global TFG pitch maximum occurs on the 8.5 m head-sea PM-Stokes record. Its change is smaller than the improvement on either high-wave following-sea record.

| TFG default regression limit | Previous | Current |
|---|---:|---:|
| Following 8.5 m JONSWAP pitch (°) | 0.09228 | 0.07897 |
| Following 8.5 m PM-Stokes pitch (°) | 0.12054 | 0.09599 |
| Following 8.5 m JONSWAP yaw (°) | 0.44170 | 0.46176 |
| Global pitch (°) | 0.23928 | 0.23880 |
| Global roll (°) | 0.53260 | 0.52151 |
| Global Y displacement (m) | 0.50981 | 0.50804 |
| Global accelerometer-bias vector (% of maximum true bias) | 107.73600 | 104.55015 |

Sixteen TFG per-record limits tightened by more than 0.5%. Two loosened by more than 0.5%. Both concern yaw on the same following-sea JONSWAP record, once for each gyro profile. Three limits in the separate legacy cubic TFG simulator also tightened by more than 0.5%.

For local reproduction, the system math library can affect the final digits. With the same compiled TFG executable and pinned inputs, Debian's system math library produced 0.09603627° pitch on the 8.5 m following PM-Stokes record. Using the Ubuntu math library from the repository's compiler qualification gave 0.09550522°, matching GitHub. All 288 gated numerical values across 32 TFG record/profile combinations matched GitHub with that library. [Compiler and math-library control](../reports/results/startup_learning/tfg_tilt_uncertainty/host-math-control.json)

[Global limits](../tests/common/CurrentReplayLimits.h) · [Per-record limits](../tests/common/RecordReplayLimits.inc) · [Exact cap changes and prior gate outcomes](../reports/results/startup_learning/tfg_tilt_uncertainty/quality-gates.json)

[Canonical four-draw scores and provenance](../reports/results/startup_learning/tfg_tilt_uncertainty/ubuntu) · [Independent host scores](../reports/results/startup_learning/tfg_tilt_uncertainty/windows.csv.gz) · [First live correction and wave-acceleration controls](startup-gravity-retention.md)

## Randomized waves and engine vibration

In a second following-sea experiment, we compared the preintegration TFG executable with the deployed startup on 90 wave histories. Ten recorded wave-phase, IMU and initialization seed triplets were reused across the eight stationary sea states and one transition. The generated input bytes matched across the two original archives for 80 histories. We regenerated the other ten once each and gave both executables the identical file. The locally generated files differed in digest from the historical files, including a separate control. On the eight original release records, the two executables independently reproduced 64 archived numerical fields per revision. Each score used 180,000 samples from 300–1,200 s.

| TFG late RMS | 90 matched following-sea histories, before → after integration |
|---|---:|
| Pitch (°) | 0.320540 → 0.320525 |
| Yaw (°) | 1.082355 → 1.082849 |
| 3-D displacement (m) | 0.487210 → 0.487127 |
| 3-D accelerometer bias (m/s²) | 0.077636 → 0.077653 |

[Matched randomized scores](../reports/results/startup_learning/tfg_tilt_uncertainty/ubuntu/randomized-paired.csv) · [Inputs, controls and binary hashes](../reports/results/startup_learning/tfg_tilt_uncertainty/ubuntu/randomized-receipt.json)

An additional stress check isolates the change from 2° to 0.5° fitted tilt uncertainty under engine vibration. We replayed the 8.5 m JONSWAP and PM-Stokes records in both directions with a 2,400 rpm input and 80 Hz bandwidth. Each row below pools four wave records at the default sensor draw with the production vibration guard. We also repeated every configuration with the guard disabled, as in the existing engine-degradation ablation. The four overlapping unguarded following-sea replays reproduced all 36 archived pitch, bias and motion measurements exactly.

| Simulated engine level | Late pitch (°), 2° → 0.5° | Late yaw (°) | Late 3-D displacement (m) | Late 3-D accel. bias (m/s²) |
|---|---:|---:|---:|---:|
| 0.60 m/s² | 1.01304 → 1.00138 | 2.77593 → 2.74699 | 0.83477 → 0.83392 | 0.20313 → 0.20098 |
| 1.20 m/s² | 1.08526 → 1.08538 | 3.17847 → 3.17697 | 0.84912 → 0.84907 | 0.22139 → 0.22126 |

[All 32 engine replays](../reports/results/startup_learning/tfg_tilt_uncertainty/ubuntu/engine-matched.csv) · [Reproduction inputs](../reports/results/startup_learning/tfg_tilt_uncertainty/ubuntu/engine-receipt.json) · [Archived engine-ablation scores](https://github.com/bareboat-necessities/ocean-imu/actions/runs/44214928822/artifacts/14391733010)

## Long-run displacement and accelerometer bias

Each row compares the pre-integration baseline with the current production configuration across all 16 records on the default sensor draw.

| Filter | 3-D displacement RMS (m), before → after | Accelerometer-bias vector RMS (m/s²), before → after |
|---|---:|---:|
| OU-II | 0.526912 → 0.526903 | 0.038953 → 0.038949 |
| OU-III | 0.433715 → 0.433715 | 0.035788 → 0.035782 |
| CARMA | 0.228164 → 0.228162 | 0.032648 → 0.032649 |
| Shared-mode CARMA | 0.205276 → 0.205282 | 0.031877 → 0.031880 |
| TFG-OU | 0.418260 → 0.418016 | 0.053934 → 0.053535 |
| TFG-CARMA | 0.228222 → 0.228237 | 0.032678 → 0.032676 |
| TFG shared-mode CARMA | 0.205369 → 0.205371 | 0.031887 → 0.031885 |
| Tangent shared-mode CARMA | 0.205354 → 0.205369 | 0.031882 → 0.031880 |

TFG-OU accelerometer-bias RMS decreased by 0.74% under the default draw. OU-III decreased by 0.016%. The shared-mode CARMA value rose by 0.009%. The large handover tilt improvement has translated into small long-run changes for most families.

## Where the X/Y advantage disappears

A separate controlled replay toggles the production startup flag on the same
16 default-seed records. All handovers occur between 31.84 and 38.29 s.
Each time window below therefore scores live output.

| Time since power-on | X displacement RMS (m), previous → deployed | Y displacement RMS (m), previous → deployed | Accel-bias RMS (m/s²), previous → deployed |
|---|---:|---:|---:|
| 40–60 s | 0.87622 → 0.84422 | 0.49928 → 0.45737 | 0.05465 → 0.05465 |
| 60–90 s | 0.53119 → 0.53552 | 0.81046 → 0.79562 | 0.05596 → 0.05596 |
| 90–120 s | 0.26636 → 0.26740 | 0.59260 → 0.59155 | 0.05655 → 0.05655 |
| 120–180 s | 0.48663 → 0.48681 | 0.39738 → 0.39736 | 0.03733 → 0.03733 |
| 300–1,200 s | 0.29893 → 0.29896 | 0.25029 → 0.25026 | 0.03579 → 0.03578 |

The previous startup also continues converging. Its roll/pitch RMS falls
from 0.28035° / 0.40726° at 40–60 s to 0.09523° / 0.13759° at 90–120 s.
The corresponding deployed values are 0.26908° / 0.39165° and
0.09522° / 0.13747°. The direction of gravity in the moving body is encoded
in the evolving attitude; known constant world gravity is used by both arms.
These measurements do not show the learned direction simply returning to its
larger initial error.

The initial X/Y advantage reaches the live estimator. The production
accelerometer-bias hold keeps the bias means identical before refinement.
The live attitude continues changing under prediction and corrections, with
the configured 2° tilt uncertainty. At about 120 s the proxy-based magnetic
refinement replaces the initial field and releases the bias hold.

We tested five mechanisms using the actual OU-III replay and shared startup
solver. These experimental policies are available in
[the reproducible driver](../tools/startup_live_retention.py).

| Final 900 s, all 16 records | 3-D displacement RMS (m) | Accel-bias RMS (m/s²) | Yaw RMS (°) |
|---|---:|---:|---:|
| Deployed startup and refinement | 0.433717 | 0.035782 | 0.780739 |
| Retain startup field, skip refinement | 0.433787 | 0.038382 | 0.778797 |
| Disable continuous magnetic calibration | 0.436112 | 0.035220 | 2.922357 |
| Joint fit at refinement, field and yaw | 0.433695 | 0.038607 | 0.769398 |
| 0.5° handover tilt uncertainty | 0.433719 | 0.035784 | 0.780580 |
| Release bias hold at accepted handover | 0.433910 | 0.037206 | 0.800998 |
| Bias uncertainty matched to configured sensor distribution | 0.433653 | 0.030832 | 0.730925 |
| Seed the private observer attitude too | 0.433711 | 0.035610 | 0.780898 |

The matched bias prior uses 0.045295 m/s² per axis, from the configured
uniform residual range ±0.0784532 m/s². The deployed OU-III prior is
0.004 m/s². This test uses distribution parameters available to the simulator
and supplies no realized bias to the estimator. It reduces pooled bias RMS
13.8% on this default-seed cohort, while X/Y changes remain mixed.
A further 96 full replays used independent seeds 7 and 41, all 16 records,
and an additional previous-startup control with the same larger bias prior.
On seed 7 the proposed prior increased bias RMS from 0.023843 to
0.029137 m/s². On seed 41 it reduced it from 0.054626 to 0.049804 m/s².
The previous-startup control with the same prior gave nearly identical late
scores. This effect is primarily a calibration-prior change, and it does not
establish a general way to preserve the gravity startup advantage.

A further 64 full replays seeded the private Mahony observer with the
accepted sensor-only startup attitude, gyro bias, or both. Its clocks, gains,
raw gyro endpoints and wave statistics were preserved. Seeding attitude
reduced late bias RMS by 0.51% following and 0.44% head-on, while 3-D RMS
changed by less than 0.002% and pooled yaw increased 0.020%. Seeding gyro
bias alone slightly worsened bias RMS. This is another mixed diagnostic,
not an adopted production change.

Earlier bias release improves the initial bias readout but increases late bias
RMS. A smaller tilt covariance alone gives negligible late changes. The
instrumented refinement's disabled control agrees with the original build
within 0.0000012 m in 3-D displacement and 0.000005° in yaw across every record.
Those differences are much smaller than the candidate regressions.

[All per-record metrics, time windows and provenance](../reports/results/startup_learning/live_retention)

## Magnetic-step consistency guard

A stationary sensor-only probe exposed a failure of the previous residual gate.
With true attitude and gyro bias both zero, a 5 µT X field step at 15 s was
accepted at 35 s and produced 2.391° tilt and 0.002384 rad/s false gyro bias
with the deployed slope weighting. A 2 µT Y step was also accepted, producing
−7.205° yaw. Rejecting the existing 10 µT control alone is insufficient.

The shared solver now compares adjacent half-second magnetic means in the
gyro-transported frame. The threshold combines a conservative estimate of
white measurement noise with the rotation permitted by the existing gyro-bias
acceptance bound. A detected discontinuity invalidates the startup fit until
the window resets; the ordinary proxy fallback remains available.

[The deterministic probe](../tests/startup_learning/field-step-probe.cpp) rejects
all 30 nonzero 1–10 µT steps and accepts the six clean controls under both
slope policies. A further 1,440 trials cover three axes, 12 step phases,
10 noise seeds, 0 or 0.8 µT white noise, and an injected gyro bias. All 720
clean controls pass and all 720 trials with a 5 µT step reject.

Across the same 960 production wave starts, the guarded solver preserves
all 888 joint handovers, all 72 deadline misses, and every measured seed
value exactly. [Measurements and source receipts](../reports/results/startup_learning/field_steps)
retain the previous failures as well as the guarded results. This is a
sensor-consistency heuristic, not a guarantee against every smooth or low-SNR
field disturbance. It repairs false acceptance; it supplies no evidence of
broad long-run displacement or accelerometer-bias gains. Final-source live
matrix, directional studies, builds and proof regeneration remain required.

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
