# Noise-free model-mismatch ablation

## Purpose

This ablation estimates the residual error floor of the deployed OU-II, OU-III,
and TFG estimators when the simulator contributes no sensor corruption.  It is
intended to separate errors that survive ideal sensor inputs from errors caused
by the stochastic IMU and magnetometer model used in the normal validation.

The experiment is deliberately a **noise-free estimator/model-mismatch floor**,
not a claim of pure plant-model mismatch.  The deployed filters are not retuned
for ideal sensors: their process and measurement covariance assumptions,
adaptation laws, pseudo-measurements, startup logic, and regularization remain
unchanged.  Consequently, the residual includes model/prior mismatch,
regularization bias, estimator-generated bias states, finite adaptation/startup
residue, attitude/translation coupling, and numerical discretization.

## Protocol

`tools/model_mismatch_ablation.py` replays the eight versioned stationary
JONSWAP and PM-Stokes cases from `oceanography-waves-lib` release `v1.1.3`, at
`H_s = {0.27, 1.5, 4.0, 8.5} m`, through each of OU-II, OU-III, and TFG.  Each
simulator is invoked with `--no-noise`; the harness also requires the simulator
to report `noise=false` before accepting its metrics.

The shared simulation runner therefore bypasses accelerometer and gyro white
noise, initial biases and bias random walks, and the magnetometer white noise,
residual bias, scale/cross-axis perturbation, and misalignment.  Magnetometer
updates remain enabled and receive the ideal simulated field.  Physical wave
motion and vessel attitude motion remain exactly those of the source records.

All RMS values use the trailing 900 s of each 1200 s replay, matching the main
validation scoring window.  Since all eight records contribute the same number
of samples, the pooled value is

`RMS_pool = sqrt(mean(RMS_i^2))`,

which is exactly the RMS obtained by concatenating the eight scored windows.

## Results

| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s^2] | Gyro bias 3D [rad/s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | 0.3328 | 0.2241 | 0.2399 | 0.4675 | 20.319 | 0.0900 | 0.2137 | 0.5560 | 0.038405 | 0.0000250 |
| OU-III | 0.2464 | 0.1736 | 0.1401 | 0.3323 | 11.867 | 0.0701 | 0.0647 | 0.2770 | 0.013235 | 0.0000243 |
| TFG | 0.5425 | 0.1993 | 0.1553 | 0.5984 | 13.152 | 0.0888 | 0.1339 | 0.3026 | 0.025373 | 0.0002315 |

The complete per-record table, raw CSV, pooled CSV, and machine-readable
manifest are committed under `reports/results/model_mismatch_ablation/`.

## Interpretation

OU-III has the lowest residual displacement floor of the three families.  Its
pooled vertical RMS is 0.1401 m and its pooled 3-D RMS is 0.3323 m.  OU-II's
corresponding residuals are 0.2399 m and 0.4675 m.  TFG is close to OU-III on
the vertical channel at 0.1553 m, but its X-axis floor is 0.5425 m, more than
twice OU-III's 0.2464 m, which drives TFG's pooled 3-D floor to 0.5984 m.

The residual vertical error scales approximately with wave amplitude within
each family rather than disappearing with ideal sensor inputs.  Across the
individual records the Z error is about 2.9--3.4% of `H_s` for OU-III, about
5.0--5.45% for OU-II, and about 2.94--3.35% for TFG.  That behavior is
consistent with a substantial deterministic estimator/regularization component
rather than an electronics-noise floor.

The zero-bias input also exposes estimator-generated bias states.  OU-III has
the smallest pooled accelerometer-bias residual, 0.0132 m/s^2, compared with
0.0384 m/s^2 for OU-II and 0.0254 m/s^2 for TFG.  TFG's pooled gyro-bias
residual is 2.31e-4 rad/s, roughly an order of magnitude above the two OU
families (~2.4--2.5e-5 rad/s).  These quantities are not sensor bias errors in
this experiment: the true injected biases are identically zero, so they measure
how much bias state the estimator creates while fitting ideal-but-model-mismatched
wave motion.

These results should not be interpreted as proving that all remaining error is
caused by the OU or TFG stochastic model alone.  In particular, the
pseudo-measurements intentionally distort low-frequency physical motion, the
adaptive operating point has finite bandwidth, attitude errors feed gravity
into horizontal acceleration, and the filter continues to use measurement
covariances calibrated for a real sensor.  The ablation measures the combined
floor of those deployed estimator choices under ideal input measurements.
