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

Source commit: `fee96e2b3be2d3d9d67b654f76ebea92b291df57`.

## Pooled results

| Arm | Axis | Offset [cm] | 3D disp [m] | 3D / CG | Max roll/pitch RMS [deg] | Tilt / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed | Tilt excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | cg | 0 | 0.3065 | 1.000x | 0.1712 | 1.000x | 0.0000 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 10 | 0.3065 | 1.000x | 0.1712 | 1.000x | 0.0386 | 0.0000 | n/a | 100.0% |
| exact | x-fore-aft | 20 | 0.3065 | 1.000x | 0.1712 | 1.000x | 0.0772 | 0.0000 | n/a | 100.0% |
| exact | x-fore-aft | 30 | 0.3065 | 1.000x | 0.1712 | 1.000x | 0.1159 | 0.0000 | 99.9% | 100.0% |
| exact | y-port | 10 | 0.3065 | 1.000x | 0.1712 | 1.000x | 0.0269 | 0.0000 | n/a | n/a |
| exact | y-port | 20 | 0.3065 | 1.000x | 0.1712 | 1.000x | 0.0538 | 0.0000 | n/a | n/a |
| exact | y-port | 30 | 0.3065 | 1.000x | 0.1712 | 1.000x | 0.0807 | 0.0000 | n/a | n/a |
| exact | z-vertical | 10 | 0.3065 | 1.000x | 0.1712 | 1.000x | 0.0462 | 0.0000 | n/a | n/a |
| exact | z-vertical | 20 | 0.3065 | 1.000x | 0.1712 | 1.000x | 0.0924 | 0.0000 | n/a | n/a |
| exact | z-vertical | 30 | 0.3065 | 1.000x | 0.1712 | 1.000x | 0.1386 | 0.0000 | n/a | 100.0% |
| gyro | x-fore-aft | 10 | 0.3062 | 0.999x | 0.1724 | 1.007x | 0.0386 | 0.0075 | n/a | 91.5% |
| gyro | x-fore-aft | 20 | 0.3060 | 0.998x | 0.1736 | 1.014x | 0.0772 | 0.0150 | n/a | 92.3% |
| gyro | x-fore-aft | 30 | 0.3057 | 0.997x | 0.1747 | 1.020x | 0.1159 | 0.0225 | 117.3% | 93.0% |
| gyro | y-port | 10 | 0.3065 | 1.000x | 0.1740 | 1.017x | 0.0269 | 0.0073 | n/a | n/a |
| gyro | y-port | 20 | 0.3065 | 1.000x | 0.1768 | 1.033x | 0.0538 | 0.0146 | n/a | n/a |
| gyro | y-port | 30 | 0.3066 | 1.000x | 0.1795 | 1.049x | 0.0807 | 0.0219 | n/a | n/a |
| gyro | z-vertical | 10 | 0.3067 | 1.000x | 0.1717 | 1.003x | 0.0462 | 0.0076 | n/a | n/a |
| gyro | z-vertical | 20 | 0.3068 | 1.001x | 0.1720 | 1.005x | 0.0924 | 0.0152 | n/a | n/a |
| gyro | z-vertical | 30 | 0.3070 | 1.002x | 0.1720 | 1.005x | 0.1386 | 0.0229 | n/a | 61.9% |
| unmodeled | x-fore-aft | 10 | 0.3077 | 1.004x | 0.1856 | 1.084x | 0.0386 | 0.0386 | n/a | n/a |
| unmodeled | x-fore-aft | 20 | 0.3094 | 1.009x | 0.2026 | 1.183x | 0.0772 | 0.0772 | n/a | n/a |
| unmodeled | x-fore-aft | 30 | 0.3113 | 1.015x | 0.2212 | 1.292x | 0.1159 | 0.1159 | n/a | n/a |
| unmodeled | y-port | 10 | 0.3063 | 0.999x | 0.1713 | 1.001x | 0.0269 | 0.0269 | n/a | n/a |
| unmodeled | y-port | 20 | 0.3063 | 0.999x | 0.1714 | 1.001x | 0.0538 | 0.0538 | n/a | n/a |
| unmodeled | y-port | 30 | 0.3066 | 1.000x | 0.1715 | 1.002x | 0.0807 | 0.0807 | n/a | n/a |
| unmodeled | z-vertical | 10 | 0.3051 | 0.995x | 0.1720 | 1.005x | 0.0462 | 0.0462 | n/a | n/a |
| unmodeled | z-vertical | 20 | 0.3040 | 0.992x | 0.1727 | 1.009x | 0.0924 | 0.0924 | n/a | n/a |
| unmodeled | z-vertical | 30 | 0.3037 | 0.991x | 0.1734 | 1.013x | 0.1386 | 0.1386 | n/a | n/a |

## Derivative band of the deployable model

Port arm at 30 cm, pooled over the same eight seas.

| Cutoff [Hz] | 3D disp [m] | 3D / CG | Residual [m/s^2] | Residual / installed |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.3119 | 1.018x | 0.0531 | 0.657 |
| 2 | 0.3079 | 1.004x | 0.0301 | 0.373 |
| 5 | 0.3065 | 1.000x | 0.0133 | 0.165 |
| 10 | 0.3065 | 1.000x | 0.0137 | 0.170 |
| 15 | 0.3066 | 1.000x | 0.0219 | 0.272 |
| 25 | 0.3066 | 1.000x | 0.0416 | 0.516 |
| 50 | 0.3083 | 1.006x | 0.0877 | 1.087 |
| 100 | 0.3150 | 1.028x | 0.1516 | 1.877 |

## Interpretation

The comparison isolates one installation effect: rigid-body rotational
acceleration at the sensor location.  No filter covariance, OU schedule,
pseudo-measurement, vibration guard, startup rule, or quality threshold is
retuned for the off-CG cases.  The exact-model arm returns to the CG baseline
to numerical precision, so the whole unmodeled penalty is deterministic and
recoverable rather than an intrinsic OU-III limit.  The gyro-derived arm shows
how much of that is available to firmware that has only the measured rate, and
the cutoff sweep measures both phase lag and differentiated gyro noise in
the residual force. A lower downstream 3D score can also reflect cancellation
with other estimator errors; its grid minimum is not a physical optimum.

The injected force is deterministic, but estimator error can increase or decrease
because it combines with existing residuals. Read every sea and axis separately;
a reduction below the CG score is not evidence that the installation is better.

Study matrix: 8 records, 3 axes, 3 offsets, 3 modelling arms.

Figures are written here and mirrored byte-for-byte into `doc/kalman_ou_iii/`
for the article.
