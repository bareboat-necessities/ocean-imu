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

Source commit: `1201f2c15274ee95fc230c77a280c2d8ec8ee47f`.

## Pooled results

| Arm | Axis | Offset [cm] | 3D disp [m] | 3D / CG | Max roll/pitch RMS [deg] | Tilt / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed | Tilt excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | cg | 0 | 0.4938 | 1.000x | 0.2276 | 1.000x | 0.0000 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 10 | 0.4938 | 1.000x | 0.2275 | 0.999x | 0.0087 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 20 | 0.4938 | 1.000x | 0.2274 | 0.999x | 0.0173 | 0.0000 | n/a | n/a |
| exact | x-fore-aft | 30 | 0.4938 | 1.000x | 0.2274 | 0.999x | 0.0260 | 0.0000 | n/a | n/a |
| exact | y-port | 10 | 0.4938 | 1.000x | 0.2273 | 0.999x | 0.0161 | 0.0000 | n/a | n/a |
| exact | y-port | 20 | 0.4938 | 1.000x | 0.2274 | 0.999x | 0.0322 | 0.0000 | n/a | n/a |
| exact | y-port | 30 | 0.4938 | 1.000x | 0.2273 | 0.998x | 0.0483 | 0.0000 | n/a | n/a |
| exact | z-vertical | 10 | 0.4938 | 1.000x | 0.2275 | 0.999x | 0.0181 | 0.0000 | n/a | n/a |
| exact | z-vertical | 20 | 0.4938 | 1.000x | 0.2276 | 1.000x | 0.0362 | 0.0000 | n/a | 101.7% |
| exact | z-vertical | 30 | 0.4938 | 1.000x | 0.2275 | 0.999x | 0.0543 | 0.0000 | n/a | 102.4% |
| gyro | x-fore-aft | 10 | 0.4938 | 1.000x | 0.2280 | 1.002x | 0.0087 | 0.0072 | n/a | n/a |
| gyro | x-fore-aft | 20 | 0.4938 | 1.000x | 0.2283 | 1.003x | 0.0173 | 0.0144 | n/a | n/a |
| gyro | x-fore-aft | 30 | 0.4938 | 1.000x | 0.2290 | 1.006x | 0.0260 | 0.0215 | n/a | n/a |
| gyro | y-port | 10 | 0.4938 | 1.000x | 0.2270 | 0.997x | 0.0161 | 0.0072 | n/a | n/a |
| gyro | y-port | 20 | 0.4939 | 1.000x | 0.2261 | 0.994x | 0.0322 | 0.0144 | n/a | n/a |
| gyro | y-port | 30 | 0.4939 | 1.000x | 0.2264 | 0.995x | 0.0483 | 0.0216 | n/a | n/a |
| gyro | z-vertical | 10 | 0.4939 | 1.000x | 0.2280 | 1.002x | 0.0181 | 0.0072 | n/a | n/a |
| gyro | z-vertical | 20 | 0.4939 | 1.000x | 0.2286 | 1.004x | 0.0362 | 0.0144 | n/a | 71.2% |
| gyro | z-vertical | 30 | 0.4940 | 1.000x | 0.2291 | 1.007x | 0.0543 | 0.0216 | n/a | 73.4% |
| unmodeled | x-fore-aft | 10 | 0.4935 | 0.999x | 0.2248 | 0.988x | 0.0087 | 0.0087 | n/a | n/a |
| unmodeled | x-fore-aft | 20 | 0.4933 | 0.999x | 0.2217 | 0.974x | 0.0173 | 0.0173 | n/a | n/a |
| unmodeled | x-fore-aft | 30 | 0.4933 | 0.999x | 0.2190 | 0.962x | 0.0260 | 0.0260 | n/a | n/a |
| unmodeled | y-port | 10 | 0.4941 | 1.001x | 0.2166 | 0.951x | 0.0161 | 0.0161 | n/a | n/a |
| unmodeled | y-port | 20 | 0.4947 | 1.002x | 0.2072 | 0.910x | 0.0322 | 0.0322 | n/a | n/a |
| unmodeled | y-port | 30 | 0.4941 | 1.001x | 0.1971 | 0.866x | 0.0483 | 0.0483 | n/a | n/a |
| unmodeled | z-vertical | 10 | 0.4899 | 0.992x | 0.2291 | 1.007x | 0.0181 | 0.0181 | n/a | n/a |
| unmodeled | z-vertical | 20 | 0.4864 | 0.985x | 0.2311 | 1.015x | 0.0362 | 0.0362 | n/a | n/a |
| unmodeled | z-vertical | 30 | 0.4830 | 0.978x | 0.2332 | 1.025x | 0.0543 | 0.0543 | n/a | n/a |

## Derivative band of the deployable model

Port arm at 30 cm, pooled over the same eight seas.

| Cutoff [Hz] | 3D disp [m] | 3D / CG | Residual [m/s^2] | Residual / installed |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.4941 | 1.000x | 0.0253 | 0.524 |
| 2 | 0.4939 | 1.000x | 0.0135 | 0.281 |
| 5 | 0.4938 | 1.000x | 0.0070 | 0.146 |
| 10 | 0.4938 | 1.000x | 0.0125 | 0.259 |
| 15 | 0.4939 | 1.000x | 0.0216 | 0.448 |
| 25 | 0.4942 | 1.001x | 0.0416 | 0.861 |
| 50 | 0.4956 | 1.004x | 0.0877 | 1.818 |
| 100 | 0.4935 | 0.999x | 0.1516 | 3.140 |

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
