# OU-III numerical source-funnel certificate

This workflow instantiates the stability objects in the OU-III manuscript against the **unchanged adaptive filter** and the same eight noisy stationary reference records used by the validation suite.

The central rule is that a stability transition is never inferred from a fitted noisy trajectory model. The certificate executable records the estimator's own closed-loop linearized error maps, the analyzer solves path-Lyapunov inequalities for those maps, and later stages use the same solved group-compatible metric for nonlinear/funnel accounting and theorem promotion.

## Claim levels

The pipeline exposes four separate statuses.

1. **Filter regression.** The unchanged adaptive OU-III implementation must pass its existing RMS/quality gates on all eight noisy records.
2. **Exact executed-word linear certificate.** Every valid ordinary-Live source word executed by those records is composed from the estimator's actual prediction, Kalman correction, pseudo-measurement, and MEKF-reset maps. Source/path Lyapunov matrices are solved from LMIs, not inferred from empirical error covariance. This level passes only when the worst generalized word factor is strictly below one.
3. **Numerical neighborhood/source-funnel certificate.** The exact group metric, nonlinear source-word margin, startup/capture funnel, hybrid jumps, and metric-dependent stochastic bounds must close numerically on a nonzero neighborhood. Executed nominal replay accounting is reported separately and cannot satisfy this level by itself.
4. **Deployment theorem certificate.** The continuous source families must be enclosed with validated arithmetic and every robust path/nonlinear/hybrid/stochastic inequality must have a strict verified margin.

Stable simulations, a solved linear LMI, and a deployment theorem are therefore never represented by the same PASS bit.

## Reference matrix

The record inventory is imported directly from `tools/ou_sweep_common.py` and is exactly:

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

The simulator uses `process_wave_file_for_tracker`, preserving the normal accelerometer/gyro/magnetometer noise, bias processes, sample rates, deterministic seeds and truth scoring.

## Estimator invariance

`tests/kalman_ou_iii/ou3-certificate-sim.cpp` is a host-only observer around `SeaStateFusion_OU_III<TrackerType::KALMANF>`. It does not retune the estimator for proof convenience:

- the linear OU-III block remains enabled;
- the adaptive tuner and production clamps remain enabled;
- the deployed `a_w` covariance synchronization policy remains enabled;
- accelerometer bias is held/released by the production startup logic;
- the complete `S=0` gain is retained, including the `S -> attitude` cross-covariance;
- no Schmidt restriction, fixed tuning, proof-specific gain or proof-specific measurement covariance is introduced.

The host translation unit exposes internal scratch matrices only to observe operations that already occurred. Production estimator semantics are unchanged.

## Exact closed-loop maps

The certificate executable consumes the actual matrices used by the MEKF:

- `F_AA_scratch_` for attitude/gyro-bias prediction;
- `F_LL_scratch_` for the exact `(v,p,S,a_w)` OU chain;
- the active/held accelerometer-bias factor;
- `PCt_scratch_ = P H^T` and `K_scratch_` for accepted updates;
- the periodic `S=0` correction with the full `P(:,S)` gain;
- the left-error covariance/reset transport after quaternion injection.

For an accepted correction the local closed-loop map is

\[
A_k^{\rm corr}=G_k(I-K_kH_k).
\]

`time_update()` performs the periodic `S=0` correction before the accelerometer correction, and the certificate map follows that implementation order. Blocks containing bias-mode or hard gauge transitions are marked hybrid and are not multiplied into ordinary same-mode words.

The analyzer checks a reconstruction residual so a stale or incorrectly ordered map cannot silently become a certificate.

## Source words and group-compatible path metric

Held and active normal-Live coordinates remain separate:

\[
e_H=(\delta\theta,b_g,v,p,S,a_w)\in\mathbb R^{18},
\qquad
e_A=(\delta\theta,b_g,v,p,S,a_w,b_a)\in\mathbb R^{21}.
\]

Source nodes retain magnetic gauge state and compact cells of the applied `(tau, sigma_aw, r_S)` schedule. Exact blocks are composed at candidate horizons 0.25, 0.5, 1, 2 and 4 s:

\[
\Phi_w=A_{k+\ell-1}^{\rm cl}\cdots A_k^{\rm cl}.
\]

The linear and nonlinear theorem now use one metric geometry. The numerical LMI is constrained to

\[
\overline P_i=\operatorname{blkdiag}\!\left(\frac{a_{R,i}}{2}I_3,P_{\xi,i}\right),
\]

which is the local quadratic of

\[
W_i(R_e,\xi)=a_{R,i}(1-\cos\theta)+\xi^TP_{\xi,i}\xi.
\]

Thus the SDP cannot exploit attitude/linear cross terms that disappear when attitude is lifted back to SO(3).

The path solver enforces

\[
\Phi_w^T\overline P_j\Phi_w-\rho\overline P_i\prec0
\]

with a strict target below one. A cutting-plane loop starts from representative words and then evaluates **every executed exact word**, adding the largest violations until the family closes or fails. The old trajectory-fit `Phi_w` and inverse-empirical-covariance metric paths are removed and prohibited by validation tests.

## Completion stage

`tools/ou3_certificate_completion.py` runs after the exact-map/LMI stage. It loads `path_metrics.npz` and evaluates the exact group metric on the eight executed trajectories.

For each executed word it reports

\[
\lambda_w^{\rm gen},\qquad
\gamma_w^{\rm replay}=\max\{0,W^+-\lambda_wW^-\},
\]

and the observed endpoint decrement. Across each fixed-dimensional mode it forms the replay disturbance envelope

\[
b_m^{\rm replay}=\frac{\gamma_m}{1-\lambda_m},\qquad \lambda_m<1,
\]

then evaluates startup handoff, bias release, magnetic lock/refinement and the finite capture recurrence

\[
c_{n+1}=\lambda_m c_n+\gamma_m.
\]

The resulting `N_H`/`T_H` is an **executed-replay funnel accounting result**. It answers whether the actual noisy trajectories are compatible with the solved metric and disturbance allowance. It is not promoted to a neighborhood theorem because the replay does not establish an infimum over nearby nonlinear states.

Outputs from this stage are:

- `completion.json`;
- `completion.md`;
- `enclosure_contract.json`.

## Theorem-promotion gate

`tools/ou3_validate_enclosure.py` is the machine gate from numerical evidence to the deployment theorem. A validated interval/Taylor-model backend must provide outward-rounded source-cell enclosures. The validator does not trust a supplied PASS flag.

For every interval word

\[
\Phi_w=C_w+\Delta_w,\qquad |\Delta_w|\le R_w,
\]

it independently checks a sound robust path-LMI bound. Since

\[
\|\Delta_w\|_2\le\|R_w\|_F=:r_w,
\]

all matrices in the box satisfy

\[
\lambda_{\max}(\Phi_w^TP_j\Phi_w-P_i)
\le
\lambda_{\max}(C_w^TP_jC_w-P_i)
+2\|P_jC_w\|_2r_w+\|P_j\|_2r_w^2.
\]

The bound must be strictly negative. In addition, the validated input must provide strict/finitely bounded values for:

- source completeness in both H and A graphs;
- finite source-prefix gain;
- `theta_star` in `(0, pi)`;
- positive lower `alpha_R` and finite upper `beta_R` for the exact SO(3) sector;
- strictly positive lower `mu_W` for every nonlinear word family;
- startup, held-to-active, magnetic-regauge, tilt-reset and cooldown inward jump margins;
- source-uniform `Sigma_bar`, `b_W`, `v_W`, and a finite-horizon stochastic failure-probability bound.

Promotion also requires provenance stating that the bounds came from validated outward-rounded arithmetic and were generated from source cells rather than trajectory fitting.

Only when all of those checks pass may the validator emit

`deployment_theorem_certificate = PASS`.

## Running

After the versioned simulation records are present under `plots/kalman_ou_ii`:

```sh
python3 tools/ou3_numerical_certificate.py \
  --output-dir reports/results/ou3_numerical_certificate

python3 tools/ou3_certificate_completion.py \
  --certificate-dir reports/results/ou3_numerical_certificate \
  --data-dir plots/kalman_ou_ii
```

The GitHub Actions workflow performs both stages automatically and uploads all reports, exact maps, traces and metrics.

When a validated continuous-source enclosure has been generated, theorem promotion is checked with:

```sh
python3 tools/ou3_validate_enclosure.py \
  --certificate-dir reports/results/ou3_numerical_certificate \
  --enclosure validated_enclosure.json
```

## Reading a failure

The first active obstruction is the scientific result:

- RMS/quality failure -> filter regression problem;
- map reconstruction residual failure -> certificate instrumentation problem;
- exact path LMI cannot reach `lambda_worst < 1` -> source partition/horizon or real local filter-dynamics problem;
- path LMI passes but neighborhood `mu_W`/`theta_star` fails -> nonlinear certificate bottleneck;
- nonlinear words pass but capture/jump margin fails -> startup/hybrid bottleneck;
- deterministic funnel passes but stochastic bound is ineffective -> concentration/localization bottleneck;
- replay/numerical layers pass but validated source-cell enclosure fails -> theorem conservatism or deployment-envelope problem.

This ordering prevents a stable-looking replay from being mislabeled as certified and prevents a bad certificate implementation from being mislabeled as filter instability.

## Ultimate question

The pipeline is designed to answer, without ambiguity:

**Does the current adaptive OU-III filter pass its ordinary performance tests, enter a numerically verified stability neighborhood on all eight noisy reference seas, and satisfy the complete source-reachable deployment theorem?**

Those are reported as separate gates. A final answer can therefore be `yes`, `no`, or `certificate inconclusive` at the exact layer where the mathematics stops, with the offending word/node and numerical margin retained in the artifact.
