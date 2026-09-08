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

Source commit: `1cb65a24379b74bf59443701c9ecb3be03f0574d`.

## Pooled results

| Arm | Axis | Offset [cm] | 3D disp [m] | 3D / CG | Max roll/pitch RMS [deg] | Tilt / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed | Tilt excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | cg | 0 | 0.4046 | 1.000x | 0.2550 | 1.000x | 0.0000 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 10 | 0.4046 | 1.000x | 0.2550 | 1.000x | 0.0298 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 20 | 0.4046 | 1.000x | 0.2550 | 1.000x | 0.0597 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 30 | 0.4046 | 1.000x | 0.2550 | 1.000x | 0.0895 | 0.0000 | n/a | n/a |
| exact | y-port | 10 | 0.4046 | 1.000x | 0.2550 | 1.000x | 0.0258 | 0.0000 | n/a | n/a |
| exact | y-port | 20 | 0.4046 | 1.000x | 0.2550 | 1.000x | 0.0516 | 0.0000 | n/a | n/a |
| exact | y-port | 30 | 0.4046 | 1.000x | 0.2550 | 1.000x | 0.0773 | 0.0000 | n/a | n/a |
| exact | z-vertical | 10 | 0.4046 | 1.000x | 0.2550 | 1.000x | 0.0386 | 0.0000 | n/a | n/a |
| exact | z-vertical | 20 | 0.4046 | 1.000x | 0.2550 | 1.000x | 0.0772 | 0.0000 | n/a | n/a |
| exact | z-vertical | 30 | 0.4046 | 1.000x | 0.2550 | 1.000x | 0.1158 | 0.0000 | n/a | n/a |
| gyro | x-fore-aft | 10 | 0.4047 | 1.000x | 0.2557 | 1.003x | 0.0298 | 0.0073 | n/a | n/a |
| gyro | x-fore-aft | 20 | 0.4047 | 1.000x | 0.2565 | 1.006x | 0.0597 | 0.0146 | n/a | n/a |
| gyro | x-fore-aft | 30 | 0.4048 | 1.000x | 0.2572 | 1.009x | 0.0895 | 0.0219 | n/a | n/a |
| gyro | y-port | 10 | 0.4046 | 1.000x | 0.2536 | 0.995x | 0.0258 | 0.0073 | n/a | n/a |
| gyro | y-port | 20 | 0.4046 | 1.000x | 0.2521 | 0.989x | 0.0516 | 0.0145 | n/a | n/a |
| gyro | y-port | 30 | 0.4046 | 1.000x | 0.2507 | 0.983x | 0.0773 | 0.0218 | n/a | n/a |
| gyro | z-vertical | 10 | 0.4047 | 1.000x | 0.2558 | 1.003x | 0.0386 | 0.0074 | n/a | n/a |
| gyro | z-vertical | 20 | 0.4048 | 1.000x | 0.2566 | 1.006x | 0.0772 | 0.0147 | n/a | n/a |
| gyro | z-vertical | 30 | 0.4049 | 1.001x | 0.2571 | 1.008x | 0.1158 | 0.0221 | n/a | n/a |
| unmodeled | x-fore-aft | 10 | 0.4043 | 0.999x | 0.2536 | 0.994x | 0.0298 | 0.0298 | n/a | n/a |
| unmodeled | x-fore-aft | 20 | 0.4041 | 0.999x | 0.2519 | 0.988x | 0.0597 | 0.0597 | n/a | n/a |
| unmodeled | x-fore-aft | 30 | 0.4042 | 0.999x | 0.2501 | 0.981x | 0.0895 | 0.0895 | n/a | n/a |
| unmodeled | y-port | 10 | 0.4048 | 1.000x | 0.2546 | 0.998x | 0.0258 | 0.0258 | n/a | n/a |
| unmodeled | y-port | 20 | 0.4053 | 1.002x | 0.2537 | 0.995x | 0.0516 | 0.0516 | n/a | n/a |
| unmodeled | y-port | 30 | 0.4060 | 1.003x | 0.2520 | 0.988x | 0.0773 | 0.0773 | n/a | n/a |
| unmodeled | z-vertical | 10 | 0.4011 | 0.991x | 0.2537 | 0.995x | 0.0386 | 0.0386 | n/a | n/a |
| unmodeled | z-vertical | 20 | 0.3981 | 0.984x | 0.2525 | 0.990x | 0.0772 | 0.0772 | n/a | n/a |
| unmodeled | z-vertical | 30 | 0.3956 | 0.978x | 0.2517 | 0.987x | 0.1158 | 0.1158 | n/a | n/a |

## Derivative band of the deployable model

Port arm at 30 cm, pooled over the same eight seas.

| Cutoff [Hz] | 3D disp [m] | 3D / CG | Residual [m/s^2] | Residual / installed |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.4086 | 1.010x | 0.0465 | 0.601 |
| 2 | 0.4055 | 1.002x | 0.0255 | 0.330 |
| 5 | 0.4045 | 1.000x | 0.0114 | 0.147 |
| 10 | 0.4045 | 1.000x | 0.0133 | 0.172 |
| 15 | 0.4046 | 1.000x | 0.0218 | 0.282 |
| 25 | 0.4050 | 1.001x | 0.0416 | 0.538 |
| 50 | 0.4069 | 1.005x | 0.0877 | 1.135 |
| 100 | 0.4103 | 1.014x | 0.1516 | 1.960 |

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
