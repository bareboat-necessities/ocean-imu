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

The canonical body directions are x = athwartships, y = fore-aft, and z = vertical.
Scoring uses the trailing **900 s** of each 1200 s record.

The `exact` arm is an oracle bound on what any lever-arm model can recover.
The `gyro` arm is the deployable one: it sees only the noisy, biased rate and
reconstructs `alpha` through a two-pole low-pass at 15 Hz
followed by a causal second-order difference.

Source commit: `dd30ffb2ea6fa2b41d8bb85dd9e5604c3735ece5`.

## Pooled results

| Arm | Axis | Offset [cm] | 3D disp [m] | 3D / CG | Max roll/pitch RMS [deg] | Tilt / CG | Installed [m/s^2] | Residual [m/s^2] | 3D excess removed | Tilt excess removed |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | cg | 0 | 0.5224 | 1.000x | 0.2362 | 1.000x | 0.0000 | 0.0000 | n/a | n/a |
| exact | x-athwartships | 10 | 0.5224 | 1.000x | 0.2362 | 1.000x | 0.1755 | 0.0000 | n/a | 100.1% |
| exact | x-athwartships | 20 | 0.5224 | 1.000x | 0.2362 | 1.000x | 0.3510 | 0.0000 | n/a | 100.0% |
| exact | x-athwartships | 30 | 0.5224 | 1.000x | 0.2362 | 1.000x | 0.5266 | 0.0000 | n/a | 100.0% |
| exact | y-fore-aft | 10 | 0.5225 | 1.000x | 0.2362 | 1.000x | 0.1312 | 0.0000 | n/a | n/a |
| exact | y-fore-aft | 20 | 0.5224 | 1.000x | 0.2362 | 1.000x | 0.2623 | 0.0000 | n/a | n/a |
| exact | y-fore-aft | 30 | 0.5224 | 1.000x | 0.2362 | 1.000x | 0.3935 | 0.0000 | 100.4% | 100.0% |
| exact | z-vertical | 10 | 0.5225 | 1.000x | 0.2362 | 1.000x | 0.2159 | 0.0000 | n/a | n/a |
| exact | z-vertical | 20 | 0.5224 | 1.000x | 0.2362 | 1.000x | 0.4317 | 0.0000 | 100.1% | 100.3% |
| exact | z-vertical | 30 | 0.5224 | 1.000x | 0.2362 | 1.000x | 0.6476 | 0.0000 | 100.0% | 100.1% |
| gyro | x-athwartships | 10 | 0.5222 | 1.000x | 0.2359 | 0.999x | 0.1755 | 0.0353 | n/a | 101.7% |
| gyro | x-athwartships | 20 | 0.5222 | 0.999x | 0.2356 | 0.998x | 0.3510 | 0.0706 | n/a | 100.4% |
| gyro | x-athwartships | 30 | 0.5222 | 1.000x | 0.2354 | 0.997x | 0.5266 | 0.1059 | n/a | 100.3% |
| gyro | y-fore-aft | 10 | 0.5227 | 1.001x | 0.2359 | 0.999x | 0.1312 | 0.0273 | n/a | n/a |
| gyro | y-fore-aft | 20 | 0.5232 | 1.002x | 0.2357 | 0.998x | 0.2623 | 0.0545 | n/a | n/a |
| gyro | y-fore-aft | 30 | 0.5240 | 1.003x | 0.2355 | 0.997x | 0.3935 | 0.0818 | 76.5% | 101.9% |
| gyro | z-vertical | 10 | 0.5225 | 1.000x | 0.2366 | 1.002x | 0.2159 | 0.0428 | n/a | n/a |
| gyro | z-vertical | 20 | 0.5227 | 1.001x | 0.2372 | 1.004x | 0.4317 | 0.0857 | 95.8% | 60.9% |
| gyro | z-vertical | 30 | 0.5231 | 1.001x | 0.2379 | 1.007x | 0.6476 | 0.1285 | 94.6% | 70.9% |
| unmodeled | x-athwartships | 10 | 0.5226 | 1.000x | 0.2550 | 1.080x | 0.1755 | 0.1755 | n/a | n/a |
| unmodeled | x-athwartships | 20 | 0.5242 | 1.003x | 0.3656 | 1.548x | 0.3510 | 0.3510 | n/a | n/a |
| unmodeled | x-athwartships | 30 | 0.5271 | 1.009x | 0.4691 | 1.986x | 0.5266 | 0.5266 | n/a | n/a |
| unmodeled | y-fore-aft | 10 | 0.5235 | 1.002x | 0.2196 | 0.930x | 0.1312 | 0.1312 | n/a | n/a |
| unmodeled | y-fore-aft | 20 | 0.5257 | 1.006x | 0.2177 | 0.922x | 0.2623 | 0.2623 | n/a | n/a |
| unmodeled | y-fore-aft | 30 | 0.5291 | 1.013x | 0.2725 | 1.154x | 0.3935 | 0.3935 | n/a | n/a |
| unmodeled | z-vertical | 10 | 0.5256 | 1.006x | 0.2368 | 1.003x | 0.2159 | 0.2159 | n/a | n/a |
| unmodeled | z-vertical | 20 | 0.5296 | 1.014x | 0.2387 | 1.011x | 0.4317 | 0.4317 | n/a | n/a |
| unmodeled | z-vertical | 30 | 0.5341 | 1.022x | 0.2420 | 1.025x | 0.6476 | 0.6476 | n/a | n/a |

## Derivative band of the deployable model

Fore-aft arm at 30 cm, pooled over the same eight seas.

| Cutoff [Hz] | 3D disp [m] | 3D / CG | Residual [m/s^2] | Residual / installed |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.5545 | 1.061x | 0.4196 | 1.066 |
| 2 | 0.5495 | 1.052x | 0.3507 | 0.891 |
| 5 | 0.5314 | 1.017x | 0.2054 | 0.522 |
| 10 | 0.5254 | 1.006x | 0.1157 | 0.294 |
| 15 | 0.5240 | 1.003x | 0.0818 | 0.208 |
| 25 | 0.5232 | 1.002x | 0.0636 | 0.162 |
| 50 | 0.5236 | 1.002x | 0.0910 | 0.231 |
| 100 | 0.5363 | 1.027x | 0.1521 | 0.386 |

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

Read the pooled ratios above with the per-sea figure beside them.  An RMS over
all eight seas is dominated by the largest, where the injected term is smallest
relative to the wave signal, so pooling understates a penalty that is severe in
the mildest seas for displacement and in the steepest seas for attitude.

Study matrix: 8 records, 3 axes, 3 offsets, 3 modelling arms.

Figures are written here and mirrored byte-for-byte into `doc/kalman_ou_iii/`
for the article.
