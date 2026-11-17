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
| baseline | cg | 0 | 0.4853 | 1.000x | 0.2027 | 1.000x | 0.0000 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 10 | 0.4853 | 1.000x | 0.2026 | 1.000x | 0.0088 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 20 | 0.4853 | 1.000x | 0.2027 | 1.000x | 0.0175 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 30 | 0.4853 | 1.000x | 0.2025 | 0.999x | 0.0263 | 0.0000 | n/a | n/a |
| exact | y-port | 10 | 0.4853 | 1.000x | 0.2025 | 0.999x | 0.0166 | 0.0000 | n/a | n/a |
| exact | y-port | 20 | 0.4853 | 1.000x | 0.2025 | 0.999x | 0.0331 | 0.0000 | n/a | n/a |
| exact | y-port | 30 | 0.4853 | 1.000x | 0.2025 | 0.999x | 0.0497 | 0.0000 | n/a | n/a |
| exact | z-vertical | 10 | 0.4853 | 1.000x | 0.2024 | 0.999x | 0.0185 | 0.0000 | n/a | n/a |
| exact | z-vertical | 20 | 0.4853 | 1.000x | 0.2025 | 0.999x | 0.0370 | 0.0000 | n/a | n/a |
| exact | z-vertical | 30 | 0.4853 | 1.000x | 0.2025 | 0.999x | 0.0556 | 0.0000 | n/a | n/a |
| gyro | x-fore-aft | 10 | 0.4853 | 1.000x | 0.2026 | 1.000x | 0.0088 | 0.0072 | n/a | n/a |
| gyro | x-fore-aft | 20 | 0.4853 | 1.000x | 0.2028 | 1.000x | 0.0175 | 0.0144 | n/a | n/a |
| gyro | x-fore-aft | 30 | 0.4854 | 1.000x | 0.2029 | 1.001x | 0.0263 | 0.0215 | n/a | n/a |
| gyro | y-port | 10 | 0.4853 | 1.000x | 0.2025 | 0.999x | 0.0166 | 0.0072 | n/a | n/a |
| gyro | y-port | 20 | 0.4854 | 1.000x | 0.2030 | 1.001x | 0.0331 | 0.0144 | n/a | n/a |
| gyro | y-port | 30 | 0.4854 | 1.000x | 0.2028 | 1.000x | 0.0497 | 0.0216 | n/a | n/a |
| gyro | z-vertical | 10 | 0.4853 | 1.000x | 0.2031 | 1.002x | 0.0185 | 0.0072 | n/a | n/a |
| gyro | z-vertical | 20 | 0.4854 | 1.000x | 0.2039 | 1.006x | 0.0370 | 0.0144 | n/a | n/a |
| gyro | z-vertical | 30 | 0.4853 | 1.000x | 0.2040 | 1.006x | 0.0556 | 0.0216 | n/a | n/a |
| unmodeled | x-fore-aft | 10 | 0.4848 | 0.999x | 0.2012 | 0.992x | 0.0088 | 0.0088 | n/a | n/a |
| unmodeled | x-fore-aft | 20 | 0.4845 | 0.998x | 0.1994 | 0.984x | 0.0175 | 0.0175 | n/a | n/a |
| unmodeled | x-fore-aft | 30 | 0.4844 | 0.998x | 0.1977 | 0.976x | 0.0263 | 0.0263 | n/a | n/a |
| unmodeled | y-port | 10 | 0.4855 | 1.000x | 0.1866 | 0.920x | 0.0166 | 0.0166 | n/a | n/a |
| unmodeled | y-port | 20 | 0.4863 | 1.002x | 0.1716 | 0.846x | 0.0331 | 0.0331 | n/a | n/a |
| unmodeled | y-port | 30 | 0.4871 | 1.004x | 0.1584 | 0.782x | 0.0497 | 0.0497 | n/a | n/a |
| unmodeled | z-vertical | 10 | 0.4812 | 0.991x | 0.2036 | 1.005x | 0.0185 | 0.0185 | n/a | n/a |
| unmodeled | z-vertical | 20 | 0.4771 | 0.983x | 0.2040 | 1.007x | 0.0370 | 0.0370 | n/a | n/a |
| unmodeled | z-vertical | 30 | 0.4737 | 0.976x | 0.2047 | 1.010x | 0.0556 | 0.0556 | n/a | n/a |

## Derivative band of the deployable model

Port arm at 30 cm, pooled over the same eight seas.

| Cutoff [Hz] | 3D disp [m] | 3D / CG | Residual [m/s^2] | Residual / installed |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.4859 | 1.001x | 0.0260 | 0.523 |
| 2 | 0.4854 | 1.000x | 0.0139 | 0.280 |
| 5 | 0.4854 | 1.000x | 0.0072 | 0.144 |
| 10 | 0.4855 | 1.000x | 0.0125 | 0.252 |
| 15 | 0.4854 | 1.000x | 0.0216 | 0.435 |
| 25 | 0.4860 | 1.001x | 0.0416 | 0.837 |
| 50 | 0.4987 | 1.028x | 0.0877 | 1.766 |
| 100 | 0.4900 | 1.010x | 0.1516 | 3.050 |

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
