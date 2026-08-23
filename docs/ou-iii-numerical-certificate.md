# OU-III numerical source-funnel certificate

This workflow turns the stability objects in the OU-III manuscript into numerical diagnostics against the **unchanged deployed filter** and the same eight noisy stationary reference records used by the validation suite.

## Scope and claim levels

The command produces two deliberately different claim levels.

1. **Finite replay certificate.**  The tool observes the exact eight noisy executions, reconstructs the full truth error state, builds source words from the branches that occurred, and checks source coverage, SO(3) chart safety, replay-bounded word inequalities, startup handoff/capture, hybrid transitions, and the ordinary OU-III RMS regression sentinels.
2. **Deployment theorem certificate.**  This remains `NOT_ESTABLISHED` until the continuous source cells are enclosed independently of the simulations.  Dense trajectories, Monte Carlo, or a fitted word map are not replacements for a validated enclosure of every source-reachable word.

The distinction is intentional.  A run may say that all eight noisy simulations stay in the candidate funnel while still refusing to claim the theorem for unobserved points in the declared source envelope.

## Reference matrix

The record inventory is imported from `tools/ou_sweep_common.py`, so it is exactly the publication/validation matrix:

| family | Hs (m) |
|---|---:|
| JONSWAP | 0.27 |
| JONSWAP | 1.50 |
| JONSWAP | 4.00 |
| JONSWAP | 8.50 |
| PM-Stokes | 0.27 |
| PM-Stokes | 1.50 |
| PM-Stokes | 4.00 |
| PM-Stokes | 8.50 |

The simulator calls the common `process_wave_file_for_tracker` path.  Consequently the accelerometer, gyro, magnetometer, bias initialization, random walk, sampling rates, default deterministic seeds, and truth scoring are the same ones used by the existing OU regression harness.

## Estimator invariance

`tests/kalman_ou_iii/ou3-certificate-sim.cpp` is an observer around the normal `SeaStateFusion_OU_III<TrackerType::KALMANF>` configuration.  It does not tune the estimator for proof convenience.  In particular:

- the linear OU-III block remains enabled;
- the online tuner and production clamps remain enabled;
- the deployed periodic `a_w` covariance synchronization policy is retained;
- accelerometer bias is held until Live exactly as in the normal configuration;
- the complete `P(:,S)` block is recorded, so the full `S -> attitude` cross-covariance is visible rather than replaced by a Schmidt/block-diagonal surrogate;
- no fixed-tuning or certificate-only gain is applied.

The extra CSV is read-only telemetry.

## Numerical objects

The trace records Live state, held/active accelerometer-bias mode, magnetic lock/refinement, accepted accelerometer/magnetometer branches, a mirror of the self-similar `S=0` cadence, applied `(tau, sigma_aw, r_S)`, innovation NIS values, the quaternion and physical states, covariance diagonals, and `P(:,S)`.

The analyzer reconstructs the error vector

\[
 e_A=(\delta\theta,b_g,v,p,S,a_w,b_a)\in\mathbb R^{21},
\]

and its held-bias counterpart

\[
 e_H=(\delta\theta,b_g,v,p,S,a_w)\in\mathbb R^{18}.
\]

Attitude is handled on SO(3), using the rotation logarithm and

\[
 V_R=1-\cos\theta.
\]

It does not propagate an Euler-angle radius.

### Source graph and words

Nodes retain the fixed-dimensional mode, magnetic gauge state, and compact cells of the applied OU/regularizer schedule.  A one-second source word also records how many accelerometer, magnetometer, and `S=0` updates occurred.  Held and active words are never multiplied into one fictitious square state transition.

For a node `i`, a positive-definite candidate metric is constructed from the normalized replay state cloud.  Repeated source words fit a map `Phi_w`; the report evaluates

\[
 \lambda_w^{\rm gen}=\lambda_{\max}
 \left(P_i^{-1/2}\Phi_w^T P_j\Phi_wP_i^{-1/2}\right).
\]

It also records replay prefix amplification, fitted-map residual bounds, and the replay decrement

\[
 \mu_{W,w}=\frac{W_i(e)-W_j(e^+)}{\|e\|^2}.
\]

These are **candidate/path diagnostics** until the fitted map and source cell are replaced by validated continuous enclosures.

### Handoff and hybrid report

The first Live sample is the measured startup handoff.  The report records its geodesic attitude error, the bias-release transition, magnetic lock/refinement transitions, and a replay capture count/time into the final inner error envelope.  An event absent from all eight records is reported as unobserved rather than silently certified.

### Stochastic report

The raw Gaussian localization calculation uses per-sensor-sigma-normalized pre-gate increments and the same quadratic-form concentration used in the paper.  Replay `b_W` and `v_W` are also supplied to the Freedman expression.  Because those two constants are extracted from the executed trajectory rather than a deterministic source-family bound, the stochastic result is labelled `DIAGNOSTIC_EMPIRICAL_BW_VW` and is not a theorem-level probability certificate.

## Running

After the versioned simulation records are present under `plots/kalman_ou_ii`:

```sh
python3 tools/ou3_numerical_certificate.py \
  --output-dir reports/results/ou3_numerical_certificate
```

The tool builds `tests/kalman_ou_iii/ou3-certificate-sim`, runs all eight records one process at a time, and writes:

- `certificate.json` — complete machine-readable numerical objects and per-word diagnostics;
- `certificate.md` — compact eight-sea summary;
- `*_certificate_trace.csv` — proof telemetry for each replay;
- `logs/*.log` — original simulator/quality-gate output.

The simulator intentionally requires one input per process so a trace can never contain two sea records with ambiguous provenance.

## Reading failures

A failed candidate is meant to expose the obstruction, not just print `FAIL`.  The JSON contains, per mode and source word:

- worst generalized contraction factor `lambda_gen`;
- source-word sample count and branch pattern;
- fitted-map residual infinity norm;
- replay disturbance allowance `gamma_replay`;
- minimum and 5th-percentile replay decrement;
- worst prefix amplification `B_pre_replay`;
- handoff/capture timing;
- stochastic localization/Freedman diagnostics.

Typical interpretations are therefore quantitative: a word may miss contraction with `lambda_gen=1.04`, a handoff may be outside the observed inner funnel, or the replay may be stable while the theorem status remains open because the continuous source cell has not yet been interval-enclosed.

## What closes the manuscript theorem later

The numerical replay machinery is intentionally reusable by a validated backend.  Promoting `deployment_theorem_certificate` from `NOT_ESTABLISHED` requires all of the following without relying on the eight trajectories as coverage evidence:

- validated continuous-source word families for both 18- and 21-state graphs;
- robust verification of the path LMIs throughout every cell;
- rigorous nodewise large-angle sectors `theta_star < pi` using the exact finite SO(3) correction;
- rigorous lower enclosure of the full nonlinear word margin `mu_W`;
- rigorous startup, bias-release, magnetic-regauge, tilt-reset, and cooldown funnel inequalities;
- non-empirical covariance/localization and martingale fluctuation bounds.

Until those objects exist, the report answers the narrower but still important engineering question: **what does the current implementation do relative to the new certificate geometry on the same eight noisy seas used to validate its performance?**
