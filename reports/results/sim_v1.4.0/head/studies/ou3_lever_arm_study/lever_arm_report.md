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

Source commit: `37d4ddb7a2d8431f91d2807c526f6adc4e22b1f8`.

## Pooled results

| Arm | Axis | Offset [cm] | 3D disp [m] | 3D / CG | Max roll/pitch RMS [deg] | Tilt / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed | Tilt excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | cg | 0 | 0.3751 | 1.000x | 0.1713 | 1.000x | 0.0000 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 10 | 0.3751 | 1.000x | 0.1713 | 1.000x | 0.0386 | 0.0000 | n/a | 100.0% |
| exact | x-fore-aft | 20 | 0.3751 | 1.000x | 0.1713 | 1.000x | 0.0772 | 0.0000 | n/a | 100.0% |
| exact | x-fore-aft | 30 | 0.3751 | 1.000x | 0.1713 | 1.000x | 0.1159 | 0.0000 | 99.9% | 100.0% |
| exact | y-port | 10 | 0.3751 | 1.000x | 0.1713 | 1.000x | 0.0269 | 0.0000 | n/a | n/a |
| exact | y-port | 20 | 0.3751 | 1.000x | 0.1713 | 1.000x | 0.0538 | 0.0000 | n/a | n/a |
| exact | y-port | 30 | 0.3751 | 1.000x | 0.1713 | 1.000x | 0.0807 | 0.0000 | n/a | n/a |
| exact | z-vertical | 10 | 0.3751 | 1.000x | 0.1713 | 1.000x | 0.0462 | 0.0000 | n/a | n/a |
| exact | z-vertical | 20 | 0.3751 | 1.000x | 0.1713 | 1.000x | 0.0924 | 0.0000 | n/a | n/a |
| exact | z-vertical | 30 | 0.3751 | 1.000x | 0.1713 | 1.000x | 0.1386 | 0.0000 | n/a | 100.1% |
| gyro | x-fore-aft | 10 | 0.3748 | 0.999x | 0.1725 | 1.007x | 0.0386 | 0.0075 | n/a | 91.4% |
| gyro | x-fore-aft | 20 | 0.3745 | 0.998x | 0.1737 | 1.014x | 0.0772 | 0.0150 | n/a | 92.2% |
| gyro | x-fore-aft | 30 | 0.3742 | 0.998x | 0.1748 | 1.021x | 0.1159 | 0.0225 | 120.9% | 92.9% |
| gyro | y-port | 10 | 0.3751 | 1.000x | 0.1741 | 1.017x | 0.0269 | 0.0073 | n/a | n/a |
| gyro | y-port | 20 | 0.3751 | 1.000x | 0.1769 | 1.033x | 0.0538 | 0.0146 | n/a | n/a |
| gyro | y-port | 30 | 0.3752 | 1.000x | 0.1797 | 1.049x | 0.0807 | 0.0219 | n/a | n/a |
| gyro | z-vertical | 10 | 0.3752 | 1.000x | 0.1718 | 1.003x | 0.0462 | 0.0076 | n/a | n/a |
| gyro | z-vertical | 20 | 0.3754 | 1.001x | 0.1721 | 1.005x | 0.0924 | 0.0152 | n/a | n/a |
| gyro | z-vertical | 30 | 0.3755 | 1.001x | 0.1722 | 1.005x | 0.1386 | 0.0229 | n/a | 57.7% |
| unmodeled | x-fore-aft | 10 | 0.3761 | 1.003x | 0.1857 | 1.084x | 0.0386 | 0.0386 | n/a | n/a |
| unmodeled | x-fore-aft | 20 | 0.3777 | 1.007x | 0.2027 | 1.183x | 0.0772 | 0.0772 | n/a | n/a |
| unmodeled | x-fore-aft | 30 | 0.3793 | 1.011x | 0.2213 | 1.292x | 0.1159 | 0.1159 | n/a | n/a |
| unmodeled | y-port | 10 | 0.3748 | 0.999x | 0.1714 | 1.001x | 0.0269 | 0.0269 | n/a | n/a |
| unmodeled | y-port | 20 | 0.3748 | 0.999x | 0.1715 | 1.001x | 0.0538 | 0.0538 | n/a | n/a |
| unmodeled | y-port | 30 | 0.3750 | 1.000x | 0.1716 | 1.002x | 0.0807 | 0.0807 | n/a | n/a |
| unmodeled | z-vertical | 10 | 0.3722 | 0.992x | 0.1721 | 1.005x | 0.0462 | 0.0462 | n/a | n/a |
| unmodeled | z-vertical | 20 | 0.3694 | 0.985x | 0.1728 | 1.009x | 0.0924 | 0.0924 | n/a | n/a |
| unmodeled | z-vertical | 30 | 0.3672 | 0.979x | 0.1735 | 1.013x | 0.1386 | 0.1386 | n/a | n/a |

## Derivative band of the deployable model

Port arm at 30 cm, pooled over the same eight seas.

| Cutoff [Hz] | 3D disp [m] | 3D / CG | Residual [m/s^2] | Residual / installed |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.3803 | 1.014x | 0.0531 | 0.657 |
| 2 | 0.3764 | 1.003x | 0.0301 | 0.373 |
| 5 | 0.3750 | 1.000x | 0.0133 | 0.165 |
| 10 | 0.3751 | 1.000x | 0.0137 | 0.170 |
| 15 | 0.3752 | 1.000x | 0.0219 | 0.272 |
| 25 | 0.3751 | 1.000x | 0.0416 | 0.516 |
| 50 | 0.3760 | 1.003x | 0.0877 | 1.087 |
| 100 | 0.3824 | 1.020x | 0.1516 | 1.877 |

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
