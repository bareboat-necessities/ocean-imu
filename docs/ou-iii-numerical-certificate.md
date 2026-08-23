# OU-III numerical source-funnel certificate

This workflow instantiates the stability objects in the OU-III manuscript against the **unchanged adaptive filter** and the same eight noisy stationary reference records used by the validation suite.

The central rule is that a stability transition is never inferred from the observed error trajectory. The certificate executable records the estimator's own closed-loop linearized error maps and the analyzer solves the path-Lyapunov inequalities for those maps.

## Claim levels

The report separates four questions that must not be conflated.

1. **Filter regression.** The unchanged adaptive OU-III implementation must pass its existing RMS/quality gates on all eight noisy records.
2. **Exact executed-word linear certificate.** Every valid ordinary-Live source word executed by those records is composed from the estimator's actual prediction, Kalman correction, pseudo-measurement, and MEKF-reset maps. Source/path Lyapunov matrices are solved from LMIs rather than from empirical error covariance. This level passes only when the worst generalized word factor is strictly below one.
3. **Numerical full source-funnel certificate.** In addition to the linear path result, the exact SO(3) sector, nonlinear word margin, startup/invariant funnel, hybrid jumps, and metric-dependent stochastic constants must close numerically.
4. **Deployment theorem certificate.** This additionally requires validated enclosure of the continuous source cells and nonlinear extrema. Eight trajectories cannot by themselves establish this level.

The distinction is intentional. Stable simulations are evidence about the filter; they are not a substitute for the inequalities in the theorem.

## Reference matrix

The record inventory is imported directly from `tools/ou_sweep_common.py`:

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

The simulator uses `process_wave_file_for_tracker`, so accelerometer, gyro and magnetometer noise, bias/random-walk processes, sampling rates, deterministic seeds and truth scoring remain those of the normal OU validation harness.

## Estimator invariance

`tests/kalman_ou_iii/ou3-certificate-sim.cpp` is a host-only observer around `SeaStateFusion_OU_III<TrackerType::KALMANF>`.

It does **not** introduce a proof-specific estimator configuration:

- the linear OU-III block remains enabled;
- the online adaptive tuner and production clamps remain enabled;
- the deployed `a_w` covariance synchronization policy remains enabled;
- accelerometer bias is held/released by the normal startup logic;
- the complete `S=0` gain is retained, including `S -> attitude` cross-covariance;
- no Schmidt restriction, fixed tuning, certificate-only gain or retuned measurement covariance is used.

The host translation unit exposes the MEKF's internal scratch matrices only so they can be recorded after the exact operations that already occurred. Production source behavior is unchanged.

## Exact closed-loop maps

The MEKF already computes the matrices required by the proof. The certificate executable reads those actual values:

- `F_AA_scratch_` for attitude/gyro-bias prediction;
- `F_LL_scratch_` for the exact `(v,p,S,a_w)` OU chain;
- the active/held accelerometer-bias prediction factor;
- `PCt_scratch_ = P H^T` and `K_scratch_` for accepted measurements;
- the exact periodic `S=0` gain, reconstructed from the same predicted covariance and `R_S` used by the filter;
- the left-error MEKF reset transport following quaternion injection.

For an accepted correction, the local closed-loop map is

\[
 A_k^{\rm corr}=G_k\,(I-K_k H_k),
\]

where `G_k` is the actual left-error reset Jacobian. The measurement Jacobian is recovered from the exact pre-update covariance by solving

\[
 P_k H_k^T=P_kC_k^T,
\]

and the reconstruction residual is recorded as a map-integrity check.

`time_update()` performs the periodic `S=0` correction internally before the accelerometer correction, so the certificate executable mirrors this exact ordering. The pseudo-measurement map uses the full `P(:,S)` gain; its attitude rows are never zeroed.

The maps are emitted in 0.25 s blocks by default. A block containing a fixed-dimensional mode/gauge transition is marked as a hybrid block and is not silently mixed into the ordinary-Live path LMI.

## Source words and path LMIs

Held-bias and active-bias ordinary-Live graphs remain separate:

\[
 e_H=(\delta\theta,b_g,v,p,S,a_w)\in\mathbb R^{18},
\]

\[
 e_A=(\delta\theta,b_g,v,p,S,a_w,b_a)\in\mathbb R^{21}.
\]

Source nodes retain the magnetic gauge state and compact cells of the applied `(tau, sigma_aw, r_S)` schedule. Consecutive exact blocks are composed into words at candidate horizons 0.25, 0.5, 1, 2 and 4 s:

\[
 \Phi_w=A_{k+\ell-1}^{\rm cl}\cdots A_k^{\rm cl}.
\]

The analyzer then solves, jointly over node metrics,

\[
 \Phi_w^T P_j\Phi_w-\rho P_i\prec0,
 \qquad P_i\succ0,
\]

with a strict target `rho < 1`. A cutting-plane loop solves a representative subset and then evaluates **every executed exact word**, adding worst violations until it closes or fails.

The numerical scales used before the LMI are coordinate conditioning only. They do not alter the estimator or create the Lyapunov metric. In particular, the old shortcuts

- fitting `Phi_w` from noisy `(X,Y)` trajectory pairs, and
- setting `P_i` to inverse empirical state covariance

have been removed and are explicitly forbidden by validation tests.

For the solved metrics the report computes

\[
 \lambda_w^{\rm gen}=\lambda_{\max}
 \left(P_i^{-1/2}\Phi_w^T P_j\Phi_wP_i^{-1/2}\right)
\]

for every executed word. The exact-linear replay gate passes only when the maximum is below one.

## SO(3), handoff and hybrid layers

Attitude is evaluated geometrically using

\[
 V_R=1-\cos\theta.
\]

The trace still records truth error, startup handoff angle, bias release, magnetic lock/refinement and hybrid transitions. These diagnostics are retained, but they are **not yet the numerical full funnel certificate**.

After the linear LMI closes, the remaining numerical obligations are:

- largest nodewise `theta_star < pi` satisfying the exact finite SO(3) measurement sector;
- positive exact nonlinear source-word infimum `mu_W`;
- startup handoff levels `c0_i`, invariant levels `b_i`, finite `N_H` and `T_H`;
- held-to-active, magnetic-regauge, tilt-reset and cooldown inequalities;
- metric-dependent, non-empirical stochastic `b_W` and `v_W` bounds.

The report therefore uses `BLOCKED_AFTER_LINEAR` rather than pretending that a linear LMI alone proves the full nonlinear/hybrid theorem.

## Running

After the versioned simulation records are present under `plots/kalman_ou_ii`:

```sh
python3 tools/ou3_numerical_certificate.py \
  --output-dir reports/results/ou3_numerical_certificate
```

CI creates an isolated Python environment containing CVXPY/SCS for the path-LMI solve.

Outputs are:

- `certificate.json` — machine-readable status and per-sea obstruction data;
- `certificate.md` — compact eight-sea summary;
- `path_metrics.npz` — solved held/active node Lyapunov matrices when available;
- `*_exact_maps.bin` — exact composed closed-loop filter maps;
- `*_certificate_trace.csv` — truth/source telemetry;
- `logs/*.log` — original simulation and quality-gate output.

## Reading a failure

The first active obstruction is the result that matters.

- Existing RMS failure: filter regression problem.
- Exact-map reconstruction residual too large: certificate instrumentation problem.
- Exact source/path LMI cannot obtain `lambda_worst < 1`: either the selected source-word horizon/partition is too restrictive or the actual local filter dynamics do not meet the linear certificate; the exact offending word is reported.
- Linear LMI passes but the SO(3)/nonlinear margin fails: nonlinear certificate is the bottleneck.
- Nonlinear words pass but `c0 -> b` capture does not close: startup/handoff is the bottleneck.
- Deterministic funnel passes but stochastic probability is useless: concentration/localization bounds are the bottleneck.

This ordering prevents a stable-looking replay from being mislabeled as certified and prevents a bad numerical identification from being mislabeled as filter instability.

## What ultimately closes the theorem

Promoting `deployment_theorem_certificate` from `NOT_ESTABLISHED` requires all numerical objects above to be replaced or surrounded by validated continuous-source enclosures:

- source-complete reachable word families for both 18- and 21-state graphs;
- robust path LMIs throughout every continuous source cell;
- rigorous nodewise large-angle sectors using exact finite SO(3) corrections;
- rigorous lower enclosure of the full nonlinear word margin;
- rigorous startup and hybrid funnel inequalities;
- non-empirical covariance/localization and martingale fluctuation bounds.

The eight noisy simulations then serve the complementary engineering question: **does the actual adaptive OU-III implementation start inside, remain inside, and enter the invariant portion of the numerically verified source funnel on the same conditions used for its performance validation?**
