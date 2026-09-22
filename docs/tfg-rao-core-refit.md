# TFG core-coefficient refit on the pinned vessel RAO records

## Decision

Retain the shipping control

\[
(c_\tau,c_\sigma,k_{S,x},k_{S,y})=(1.0,0.8,1.15,1.15),
\qquad Q_{b_g}=3.75\times10^{-10}.
\]

The completed study found real metric tradeoffs and an asymmetric physical-performance basin, but **no tested replacement passed the unchanged promotion and deterministic regression requirements**. This is not a proof that the legacy constants are globally optimal. Historical pre-RAO comments are not evidence that these values were fitted to the current vessel dataset.

The best non-baseline arm in the final paired, gate-aware selection pool was frozen at `(1.0, 0.9, 1.05, 1.2)` before the sealed evaluation. It was already ineligible for shipping because of deterministic regression failures. Its sealed comparison is retained as negative promotion evidence, not presented as a successful retune. No second candidate was chosen after observing the holdout.

OU-III, the common wave-period front end, the process-noise improvement from PR #572, and every quality limit remain unchanged. Neither the completed generic tuning campaign nor the attitude/bias tuning campaign was repeated.

## Histories, objective and provenance

Every coefficient arm uses all eight pinned 28-ft vessel-RAO records: JONSWAP and PM-Stokes at significant wave heights 0.27, 1.5, 4.0 and 8.5 m. Each replay is scored over the final 900 seconds, or 180,000 samples. The dataset loader verifies the archive and individual input hashes.

Training seeds are `default, 11, 23`; independent refinement seeds are `31, 47, 59`. The sealed seeds are `104729, 130363, 155921, 196613`. A nondefault seed sets both the sensor-noise and initialization streams. The final paired audit verified the exact logged noise/initialization seed lines, input hashes and scoring windows for all 32 holdout histories across TFG baseline, frozen challenger and shipping OU-III.

The selection score is the equal-history average of weighted log error ratios against the paired baseline. Weights are 2 for heave percentage of Hs and 3-D displacement RMS; 0.5 for roll, pitch and gyro-bias RMS; 0.25 for yaw; and 0.75 for accelerometer-bias RMS. Each additional violation on a paired history adds 0.025. Net violation totals cannot hide a new failure on a formerly better history. Missing or nonfinite metrics make an arm ineligible while retaining its evidence.

Pooled RMS means `sqrt(mean(case_RMS^2))`, not the arithmetic mean of case RMS values. Both are reported, together with equal-case geometric changes. Before holdout, the protocol explicitly clarified this distinction and conservatively applied the original promotion guards to both aggregations. The objective and numerical thresholds were not relaxed.

| Executed phase | Configurations, including control where applicable | New replays |
|---|---:|---:|
| Broad axes and joint screen | 87 | 2,088 |
| Coupled boundary extension | 22 | 528 |
| Independent four-dimensional refinement | 82 | 1,968 |
| Finer asymmetric extension | 55 | 1,296; 24 control cases reused |
| Low/high-regime cross-check | 2 challengers | 48 |
| Interior, sigma bracket, symmetry and reversed-axis controls | 6 challengers | 144 |
| Interior training confirmation | 1 challenger | 24 |
| Nearby training feasibility checks | 4 challengers | 96 |
| Final bounded near-control checks on all six selection seeds | 4 challengers | 192 |
| Sealed TFG control/challenger and OU-III comparison | 3 arms | 96 |
| Subsequent passive allocation diagnostic | 3 arms | 24 plus one trace-on/off check |

The coefficient-selection phases comprise 6,384 executed replays and 253 distinct coefficient tuples. The finer search reused its existing control and executable rather than restarting a completed campaign. Operational compilation, portability and deterministic validation checks are separate from these counts.

GitHub Actions runs `35669711189`, `35672004933` and `35673733797` contain the executed broad, extended, four-dimensional and finer response surfaces. Their merged case records, including all failed arms, are retained under `reports/results/tfg_rao_core_refit/{broad,extend,refine,fine}`. Local continuation used the same portable TFG executable produced by run `35672004933`, SHA-256 `214b746ed8eb89ff1c5b9713a834e6cb243eae0f6a261e35b876224db989b8e9`. Its source commit is `5167aa6d218dbb19c0baf69c1f75cb71eb9bfa5a`; the local checkout provenance is recorded separately rather than misattributed to that build commit.

The challenger was frozen in PR #571 commit `d325b26fe66937346ef57bf98bac7edc2f43093c`, before any sealed replay. The coefficient hash is `206cab6cef27cd3f5c09d7dc181fc3013eac205507f6666f7ebb112478783496`. The unchanged portable OU-III comparison executable has SHA-256 `28733e2037ac999fe35f3466656756a7c8828db9f8cfb71ccb89e22ad5286cbe`.

## Surface, boundaries and interactions

The explored envelope reaches `c_tau = 0.45–1.85`, `c_sigma = 0.35–3.5` and horizontal factors `0.30–2.70`. The broad screen includes 64 joint tau/sigma/symmetric-XY combinations, not only one-dimensional sweeps. Later refinement explicitly separates X and Y.

The broad gate-aware score favors the control. Its best unpenalized improvement is only about 0.039%, with additional violations and worse pooled displacement behavior. Larger sigma values can improve bias metrics while worsening pitch or heave and producing more failures. All 21 non-control boundary-extension arms are worse in the balanced physical score even before gate penalties.

Tau and horizontal regularization interact strongly. Away from cadence and regularizer clamps, the SpectralMSE expression supplies the search coordinate `k_xy proportional to tau^(-41/14)`. This is not an invariance of the nonlinear filter: process covariance, vertical regularization, cadence limits and the lower rS bound still change. The exploratory broad and four-dimensional quadratic fits have R-squared about 0.991 and 0.994, respectively; they summarize interactions, not certify an optimum.

The first independent 81-point four-dimensional grid contains 54 asymmetric and 27 symmetric tuples. Its best asymmetric arm `(1.025, 1.05, 0.85, 1.45)` improves the balanced refinement score by about 0.88%, versus about 0.38% for the best symmetric arm. Because X and Y touch opposite grid faces, a further 54-point joint grid extends X down to 0.60 and Y up to 1.70 while refining tau and sigma.

The extended basin moves into the X/Y interior, near X=0.75 and Y=1.3–1.4. A separately evaluated rounded interior point `(1.03, 1.10, 0.75, 1.30)` improves the balanced refinement score by about 1.01%. Sigma brackets 1.0, 1.1, 1.2 and 1.4 turn over near 1.1 at that point. An equal-product symmetric control and a reversed-X/Y control are both worse, demonstrating an orientation effect rather than only a geometric-mean regularization effect.

However, that interior point fails four of the eight deterministic records, including low-wave heave and large-wave displacement/bias bars. Its seven additional paired selection violations reject promotion. Four nearby training checks also fail unchanged deterministic bars. The final bounded near-control probe does not find a feasible replacement. The control remains the constrained selection winner.

## Sealed comparison

Negative changes mean lower error. Record wins compare each record after pooling its four heldout seeds. All values below are physical errors, not the selection score.

| Metric | Retained TFG pooled RMS | Frozen challenger pooled RMS | Change | Record wins | Paired wins |
|---|---:|---:|---:|---:|---:|
| Heave RMS, m | 0.177888059 | 0.177906902 | +0.0106% | 2/8 | 8/32 |
| Heave RMS, %Hs | 3.916381569 | 3.918334380 | +0.0499% | 2/8 | 8/32 |
| 3-D displacement RMS, m | 0.410512283 | 0.409099005 | -0.3443% | 8/8 | 28/32 |
| Roll RMS, degrees | 0.514092593 | 0.510139215 | -0.7690% | 4/8 | 17/32 |
| Pitch RMS, degrees | 0.413952161 | 0.402497984 | -2.7670% | 8/8 | 24/32 |
| Yaw RMS, degrees | 3.175740916 | 3.151566099 | -0.7612% | 8/8 | 24/32 |
| Accelerometer-bias RMS, m/s^2 | 0.111928911 | 0.110190810 | -1.5529% | 8/8 | 27/32 |
| Gyro-bias RMS, rad/s | 0.000105476170 | 0.000105399930 | -0.0723% | 8/8 | 26/32 |

The balanced composite improves by 0.5784%, below the predeclared 1% promotion threshold. Heave worsens under both pooled and mean-case aggregation and wins only two records, below the five-record requirement. These are the sealed rejection reasons. There are no additional paired holdout violations, but this does not erase the earlier deterministic failures.

All 32 randomized holdout histories fail at least one existing regression bar for each of the TFG control, challenger and OU-III. These failures are retained, not excluded from averages or relabeled as passes. The deterministic regression bars are not a probabilistic performance certificate over arbitrary initialization/noise seeds.

| Pooled metric | Shipping OU-III | Retained TFG change versus OU-III | Challenger change versus OU-III |
|---|---:|---:|---:|
| Heave RMS, %Hs | 3.858931473 | +1.489% | +1.539% |
| 3-D displacement RMS, m | 0.399279209 | +2.813% | +2.459% |
| Roll RMS, degrees | 0.340549277 | +50.960% | +49.799% |
| Pitch RMS, degrees | 0.381783483 | +8.426% | +5.426% |
| Yaw RMS, degrees | 2.982075622 | +6.494% | +5.684% |
| Accelerometer-bias RMS, m/s^2 | 0.087151184 | +28.431% | +26.436% |
| Gyro-bias RMS, rad/s | 0.000122309481 | -13.763% | -13.825% |

TFG does not uniformly outperform OU-III on these histories. Its gyro-bias advantage persists, while displacement, attitude and accelerometer-bias tradeoffs remain. These are eight fixed RAO motion histories with independent sensor/initialization seeds, not independent vessel or wave-realization validation.

## Sea-state regimes and the existing adaptation law

The low-sea broad preference for tighter horizontal factors and the high-sea preference for a different tau/sigma combination persist as opposing tradeoffs on fresh refinement seeds. The low-sea proposal improves small-wave 3-D displacement but worsens larger-wave displacement by roughly 10–15%. The high-sea proposal improves some larger-wave records while worsening low-wave heave and displacement. A single new constant is not justified by that conflict.

The passive continuation additionally finds the lower rS bound active for about 20.75% and 25.80% of the two 0.27-m control records at diagnostic seed 31, but not for the larger records. This identifies a meaningful regime boundary in the existing adaptation law. A future modification should be expressed through causal onboard quantities such as estimated motion amplitude and regularizer saturation, not true Hs or dataset record names. No record-specific oracle or unvalidated adaptive law is shipped here.

Across all 87 broad coefficient arms, the reported wave period is identical for each paired history. This study identifies no reason to retune the shared `WavePeriodEstimator`.

## Attitude / bias / wave-acceleration allocation

After closing the coefficient selection and sealed decision, 24 passive diagnostic replays compare retained TFG, the sealed challenger and the already-rejected interior point on all eight records at seed 31. The test-only reader uses the actual retained accelerometer `H`, `K` and residual `r`:

\[
G_j=H_jK_j,\qquad u_j=G_jr,
\quad j\in\{\theta,a_w,b_a\}.
\]

All component responses therefore share accelerometer residual units. Individual squared response energies and signed cross terms are retained; their sum is checked against the net response energy. Signed trace fractions are not probabilities or fractions of Fisher information.

For retained TFG, the a_w block accounts for 99.9289–99.9894% of the signed operator trace across the eight records. The attitude share is 0.0103–0.0695%, and the accelerometer-bias share is 0.000282–0.001582%. In the component-energy decomposition, the a_w share is above 99.9999%; the separately reported cross terms prevent a misleading additive-energy interpretation. The H K block-decomposition closure is exact in the recorded double-precision analysis.

This establishes a strong instantaneous allocation toward the modeled wave acceleration, which remains after retuning attempts. It does not establish temporal unobservability or a covariance-reset defect. Magnetic and integral updates, model dynamics and covariance transport can still carry information not represented by this accelerometer-only readout. No estimator change is inferred solely from these fractions.

The reader samples at 20 Hz over the same final 900-second interval. One additional same-binary trace-on/off replay produced identical parsed physical metrics and identical quality-gate results. The diagnostic executable hash is `bd4faf48c672e3f09b701064e9950212b1832fd81b5a3ec8e972c8944e101e2d`.

## Reproduction and evidence retention

`tools/tfg_rao_refit_protocol.json` contains the frozen selection/promotion rules. `tools/tfg_rao_refit_analysis.py` produces the full metric surfaces, boundary contacts, regime summaries and paired scores. `tools/tfg_rao_refit_report.py` evaluates the sealed result without selecting another candidate. `tools/tfg_information_allocation.py` reads the actual passive diagnostic, including its trace-on/off check.

Existing completed phases should be consumed from their retained case records rather than restarted. The replay cache verifies executable, producer, input, seed, effective override, payload and log hashes. Every attempted case, including fatal execution errors or gate failures, has a retained result. The supplemental case archive preserves all requested physical endpoints and exact violation messages for the locally continued phases; the full downloadable evidence bundle additionally retains their original parsed results and logs.
