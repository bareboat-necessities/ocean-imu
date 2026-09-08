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
JONSWAP and PM-Stokes cases from `oceanography-waves-lib` release `v1.2.1`, at
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

Current numerical results, all 24 per-record rows, figures and the exact
input hashes are in [the generated report](../reports/results/model_mismatch_ablation/model_mismatch_report.md)
and its adjacent manifest. The article table and figures mirror this run.
The incident-sea heights identify the RAO forcing; the reference displacement
is vessel CG displacement. The largest incident seas remain extrapolated
stress tests for this estimated 28 ft sailboat.

The noise-free residual includes estimator regularization, finite adaptation,
attitude coupling and discretization. It is not a measured hydrodynamic model
error or a uniform theorem bound. A nonzero estimated bias is an estimation
error in this experiment, because true injected bias is zero.
