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

The same run writes the two published SVG figures alongside the tables, using a
fixed Matplotlib hash salt and no creation timestamp so repeated runs on the
same evidence produce byte-identical files.  `--no-plots` skips them for
environments without Matplotlib.

## Results

| Family | X disp [m] | Y disp [m] | Z disp [m] | 3D disp [m] | Z / ref RMS [%] | Roll [deg] | Pitch [deg] | Yaw [deg] | Acc bias 3D [m/s^2] | Gyro bias 3D [rad/s] |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OU-II | 0.3578 | 0.2478 | 0.2372 | 0.4957 | 20.088 | 0.1147 | 0.2277 | 0.5441 | 0.042484 | 0.0000112 |
| OU-III | 0.2772 | 0.2003 | 0.1392 | 0.3692 | 11.790 | 0.1151 | 0.0733 | 0.2609 | 0.021110 | 0.0000097 |
| TFG | 0.5388 | 0.2043 | 0.1494 | 0.5952 | 12.650 | 0.1179 | 0.1839 | 0.4806 | 0.031267 | 0.0002692 |

The complete per-record table, raw CSV, pooled CSV, machine-readable manifest,
and the two published figures are committed under
`reports/results/model_mismatch_ablation/`.

## Figures

`ou_model_mismatch_floor.svg` shows the pooled floor per family across the
displacement channels, the attitude channels, and the two estimator-generated
bias states.  The gyro-bias panel is logarithmic because TFG exceeds the OU
families by more than an order of magnitude on that channel.

`ou_model_mismatch_scaling.svg` shows the per-record residual against `H_s`
on the four heights, separating JONSWAP (solid, filled markers) from PM-Stokes
(dashed, open markers).  The two logarithmic panels carry a slope-one guide;
the center panel normalizes the vertical residual by `H_s`.

Both files are mirrored byte-for-byte into `doc/kalman_ou_iii/` and included in
the OU-III article, which reproduces the pooled table and both figures in
Sec. "Noise-Free Model-Mismatch Floor".

## Interpretation

OU-III has the lowest residual displacement floor of the three families.  Its
pooled vertical RMS is 0.1392 m and its pooled 3-D RMS is 0.3692 m.  OU-II's
corresponding residuals are 0.2372 m and 0.4957 m.  TFG is close to OU-III on
the vertical channel at 0.1494 m, but its X-axis floor is 0.5388 m, nearly
twice OU-III's 0.2772 m, which drives TFG's pooled 3-D floor to 0.5952 m.

On attitude the three families share essentially one roll floor
(0.115-0.118 deg).  The separation is in pitch and yaw, where OU-III reaches
0.0733 deg and 0.2609 deg against 0.2277 deg and 0.5441 deg for OU-II and
0.1839 deg and 0.4806 deg for TFG.

The residual vertical error scales approximately with wave amplitude within
each family rather than disappearing with ideal sensor inputs.  Across the
individual records the Z error is about 2.83-3.41% of `H_s` for OU-III, about
4.91-5.44% for OU-II, and about 3.02-3.62% for TFG, over a 31:1 range of `H_s`.
A sensor-noise-limited floor would instead be roughly fixed in absolute terms
and would fall as a fraction of `H_s` as the sea grows.  The observed
amplitude-proportional behavior is consistent with a substantial deterministic
estimator/regularization component rather than an electronics-noise floor.

The zero-bias input also exposes estimator-generated bias states.  OU-III has
the smallest pooled accelerometer-bias residual, 0.0211 m/s^2, compared with
0.0425 m/s^2 for OU-II and 0.0313 m/s^2 for TFG.  TFG's pooled gyro-bias
residual is 2.69e-4 rad/s, more than an order of magnitude above the two OU
families (9.70e-6 and 1.12e-5 rad/s).  These quantities are not sensor bias
errors in this experiment: the true injected biases are identically zero, so
they measure how much bias state the estimator creates while fitting
ideal-but-model-mismatched wave motion.

These results should not be interpreted as proving that all remaining error is
caused by the OU or TFG stochastic model alone.  In particular, the
pseudo-measurements intentionally distort low-frequency physical motion, the
adaptive operating point has finite bandwidth, attitude errors feed gravity
into horizontal acceleration, and the filter continues to use measurement
covariances calibrated for a real sensor.  The ablation measures the combined
floor of those deployed estimator choices under ideal input measurements.
