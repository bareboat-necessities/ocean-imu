# OU-III IMU lever-arm installation study

The wave records define motion at the vessel centre of gravity (CG).  A rigidly
mounted IMU displaced by body-frame vector $r$ additionally measures
`alpha x r + omega x (omega x r)`.  The simulator applies that term to the
accelerometer truth before sensor corruption, and applies the filter's own
lever-arm model after corruption, immediately before fusion.  Nothing else
changes: the noise realization, OU-III configuration, adaptation, startup
logic, pseudo-measurements, vibration guard, and scoring are identical in
every arm.

| Arm | What the filter receives |
| --- | --- |
| baseline | IMU at the CG |
| unmodeled | off-CG specific force, no filter-side model |
| gyro | off-CG specific force, compensated from the measured rate |
| exact | off-CG specific force, compensated from truth kinematics |

The canonical body directions are x = forward, y = port, and z = vertical.
Scoring uses the trailing **900 s** of each 1200 s record.

The `exact` arm is an oracle bound on what any lever-arm model can recover.
The `gyro` arm is the deployable one: it sees only the noisy, biased rate and
reconstructs `alpha` through a two-pole low-pass at 15 Hz
followed by a causal second-order difference.

Source commit: `5c19226291a2beb735906009d5eaf776cc5e0b0f`.

## Pooled results

| Arm | Axis | Offset [cm] | 3D disp [m] | 3D / CG | Max roll/pitch RMS [deg] | Tilt / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed | Tilt excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | cg | 0 | 0.4052 | 1.000x | 0.2515 | 1.000x | 0.0000 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 10 | 0.4052 | 1.000x | 0.2515 | 1.000x | 0.0298 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 20 | 0.4052 | 1.000x | 0.2515 | 1.000x | 0.0597 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 30 | 0.4052 | 1.000x | 0.2515 | 1.000x | 0.0895 | 0.0000 | n/a | n/a |
| exact | y-port | 10 | 0.4052 | 1.000x | 0.2515 | 1.000x | 0.0258 | 0.0000 | n/a | n/a |
| exact | y-port | 20 | 0.4052 | 1.000x | 0.2515 | 1.000x | 0.0516 | 0.0000 | n/a | n/a |
| exact | y-port | 30 | 0.4052 | 1.000x | 0.2515 | 1.000x | 0.0773 | 0.0000 | n/a | n/a |
| exact | z-vertical | 10 | 0.4052 | 1.000x | 0.2515 | 1.000x | 0.0386 | 0.0000 | n/a | n/a |
| exact | z-vertical | 20 | 0.4052 | 1.000x | 0.2515 | 1.000x | 0.0772 | 0.0000 | n/a | n/a |
| exact | z-vertical | 30 | 0.4052 | 1.000x | 0.2515 | 1.000x | 0.1158 | 0.0000 | n/a | n/a |
| gyro | x-fore-aft | 10 | 0.4053 | 1.000x | 0.2521 | 1.003x | 0.0298 | 0.0073 | n/a | n/a |
| gyro | x-fore-aft | 20 | 0.4055 | 1.001x | 0.2529 | 1.006x | 0.0597 | 0.0146 | n/a | n/a |
| gyro | x-fore-aft | 30 | 0.4056 | 1.001x | 0.2536 | 1.008x | 0.0895 | 0.0219 | n/a | n/a |
| gyro | y-port | 10 | 0.4052 | 1.000x | 0.2501 | 0.995x | 0.0258 | 0.0073 | n/a | n/a |
| gyro | y-port | 20 | 0.4051 | 1.000x | 0.2486 | 0.989x | 0.0516 | 0.0145 | n/a | n/a |
| gyro | y-port | 30 | 0.4051 | 1.000x | 0.2473 | 0.984x | 0.0773 | 0.0218 | n/a | n/a |
| gyro | z-vertical | 10 | 0.4053 | 1.000x | 0.2523 | 1.003x | 0.0386 | 0.0074 | n/a | n/a |
| gyro | z-vertical | 20 | 0.4054 | 1.000x | 0.2531 | 1.007x | 0.0772 | 0.0147 | n/a | n/a |
| gyro | z-vertical | 30 | 0.4055 | 1.001x | 0.2537 | 1.009x | 0.1158 | 0.0221 | n/a | n/a |
| unmodeled | x-fore-aft | 10 | 0.4048 | 0.999x | 0.2501 | 0.994x | 0.0298 | 0.0298 | n/a | n/a |
| unmodeled | x-fore-aft | 20 | 0.4046 | 0.998x | 0.2485 | 0.988x | 0.0597 | 0.0597 | n/a | n/a |
| unmodeled | x-fore-aft | 30 | 0.4047 | 0.999x | 0.2468 | 0.981x | 0.0895 | 0.0895 | n/a | n/a |
| unmodeled | y-port | 10 | 0.4055 | 1.001x | 0.2510 | 0.998x | 0.0258 | 0.0258 | n/a | n/a |
| unmodeled | y-port | 20 | 0.4059 | 1.002x | 0.2501 | 0.994x | 0.0516 | 0.0516 | n/a | n/a |
| unmodeled | y-port | 30 | 0.4066 | 1.003x | 0.2484 | 0.988x | 0.0773 | 0.0773 | n/a | n/a |
| unmodeled | z-vertical | 10 | 0.4024 | 0.993x | 0.2502 | 0.995x | 0.0386 | 0.0386 | n/a | n/a |
| unmodeled | z-vertical | 20 | 0.4001 | 0.987x | 0.2490 | 0.990x | 0.0772 | 0.0772 | n/a | n/a |
| unmodeled | z-vertical | 30 | 0.3983 | 0.983x | 0.2483 | 0.987x | 0.1158 | 0.1158 | n/a | n/a |

## Derivative band of the deployable model

Port arm at 30 cm, pooled over the same eight seas.

| Cutoff [Hz] | 3D disp [m] | 3D / CG | Residual [m/s^2] | Residual / installed |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.4107 | 1.013x | 0.0465 | 0.601 |
| 2 | 0.4060 | 1.002x | 0.0255 | 0.330 |
| 5 | 0.4049 | 0.999x | 0.0114 | 0.147 |
| 10 | 0.4049 | 0.999x | 0.0133 | 0.172 |
| 15 | 0.4051 | 1.000x | 0.0218 | 0.282 |
| 25 | 0.4056 | 1.001x | 0.0416 | 0.538 |
| 50 | 0.4076 | 1.006x | 0.0877 | 1.135 |
| 100 | 0.4154 | 1.025x | 0.1516 | 1.960 |

## Interpretation

The comparison isolates one installation effect: rigid-body rotational
acceleration at the sensor location.  No filter covariance, OU schedule,
pseudo-measurement, vibration guard, startup rule, or quality threshold is
retuned for the off-CG cases.  The exact-model arm returns to the CG baseline
to numerical precision, so the whole unmodeled penalty is deterministic and
recoverable rather than an intrinsic OU-III limit.  The gyro-derived arm shows
how much of that is available to firmware that has only the measured rate, and
the cutoff sweep shows that its single design parameter is two-sided: too
narrow a derivative band and the low-pass phase lag misaligns a correction of
the right size, too wide and differentiated gyro noise dominates.

The injected force is deterministic, but estimator error can increase or decrease
because it combines with existing residuals. Read every sea and axis separately;
a reduction below the CG score is not evidence that the installation is better.

Study matrix: 8 records, 3 axes, 3 offsets, 3 modelling arms.

Figures are written here and mirrored byte-for-byte into `doc/kalman_ou_iii/`
for the article.
