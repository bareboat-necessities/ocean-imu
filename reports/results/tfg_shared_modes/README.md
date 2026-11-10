# TFG shared-mode CARMA: measured results

The 27-state shared-mode TFG filter is the most accurate of the four
estimators in the condition-balanced position comparison. Most of its gain
comes from the shared stochastic prior; the TFG geometry adds a small further
improvement at about twice the host cost of shared QMEKF. It remains an
experimental estimator: the shipping OU-III default is unchanged.

| Shared TFG versus | Heave RMS change | Horizontal RMS change |
|---|---:|---:|
| QMEKF shared | -0.04% [-0.06, -0.02] | -0.28% [-0.40, -0.16] |
| TFG independent | -5.93% [-6.19, -5.67] | -2.81% [-2.93, -2.68] |
| QMEKF independent | -5.95% [-6.20, -5.69] | -3.09% [-3.23, -2.95] |

Changes are condition-balanced geometric means over 24 independent paired
seeds, with paired 95% intervals. Each of the 27 conditions has equal
log-ratio weight. The primary study runs 1200 seconds at 200 Hz and scores
from 300 seconds. There are 6,052 estimator trajectories in all, including
registered ablations, reset diagnostics, pinned replays, transition
diagnostics and the two exploratory reduced-heave studies
(`experiment_counts.json`).

Per condition, shared TFG has the lowest mean heave RMS of the four filters in
22/27 conditions and is at most 0.031% above the lowest in the others. It has
the lowest horizontal RMS in 17/27, within 0.08% of the lowest in all but
one. The exception is colored sensor noise, where its horizontal error is
0.48% above shared QMEKF. No primary practical improvement over shared QMEKF
reaches the predeclared 5% margin: geometry alone is a small effect.

## Reduced heave and the excitation cap

The shared prior maps the measured band cross-spectrum to a Hermitian
excitation `[[H, b], [b^H, sigma_z^2]]`. It keeps the measured
horizontal-on-vertical regression in `b`, with the horizontal prior
`H = sigma_h^2 I`. Its applied coherence `q = b^H H^-1 b / sigma_z^2` grows as
`1/k^2` when the heave RAO is scaled by `k`.

- **Legacy cap.** Above `q = 0.999` the legacy mapping scaled `b` down. That
  understates the regression and leaves a nearly singular prior: only 0.1% of
  the along-wave horizontal variance stays independent of heave. Under reduced
  heave (`k = 0.15`) this costs both shared filters their heading. In the
  12-seed diagnosis shared QMEKF/TFG reach 21.7/25.7 mm heave, 90.1/56.3
  degrees mean heading RMS and 12/12 and 7/12 heading failures above 30
  degrees.
- **Default policy.** Keep `b` and, when `q > 0.8`, enlarge `H` only along the
  coupling: `H <- H + ((lambda - 1) / (q sigma_z^2)) b b^H` with
  `lambda = q / 0.8`. The applied coherence is then 0.8 and the cross-wave
  prior is unchanged. In the primary reduced-heave condition shared TFG has
  4.26/20.92 mm heave/horizontal error, against 5.39/23.39 mm for
  independent TFG, with no heading failure in 24 runs.
- **Controls (12 seeds).** Cap targets 0.7 and 0.9 give 4.21 and 5.10 mm
  shared-TFG heave, against 4.38 mm at 0.8. Without confidence weighting,
  estimated coupling reaches 3.71/19.78 mm and known coupling 3.79/19.74 mm.
  The estimator therefore recovers the known-coupling benefit, and confidence
  weighting gives part of it up here. Real-only coupling gives 6.98/23.64 mm,
  worse than independent axes, so the quadrature term is required.
  Omitting reset transport is immaterial in this condition.
- **Design choice.** Enlarging the whole horizontal block by `lambda`, instead
  of only the coupling direction, also inflates the cross-wave prior. In
  development runs of the `k = 0.25-0.4` dose-response cells (seeds
  1001-1012) it made shared horizontal error 3-5% worse than independent TFG;
  the table below shows the directional form. The target 0.8 was not retuned
  for it.
- **Development overlap.** The target 0.8 was chosen on seeds 1001-1003 of the
  `k = 1, 0.5, 0.3, 0.15` cells. Those seeds are also primary seeds, and
  `k = 1` and `k = 0.15` are the primary `jonswap` and `weak_heave`
  conditions. The cap does not bind in any other primary condition: every
  other primary and ablation row is identical with and without it.

### Dose response

`tools/tfg_shared_modes_heave_sweep.py` scales the heave RAO by `k` on the
first 12 primary seeds; horizontal motion and sensor streams are unchanged.
The repeated `k=1` and `k=0.15` cells reproduce all 132 overlapping
published runs exactly.

| k | Legacy-cap coherence | Shared TFG, legacy cap (heave / horizontal mm, heading deg) | Shared TFG, default policy | Independent TFG |
|---:|---:|---:|---:|---:|
| 1 | 0.334 | 13.8 / 21.3, 1.8 | 13.8 / 21.3, 1.8 | 14.6 / 21.7, 1.8 |
| 0.6 | 0.936 | 8.8 / 18.4, 1.8 | 8.8 / 18.2, 1.8 | 10.6 / 18.6, 1.8 |
| 0.5 | 0.999 (cap) | 15.2 / 27.4, 2.5 | 7.7 / 17.8, 1.8 | 9.6 / 18.0, 1.8 |
| 0.4 | 0.999 (cap) | 22.2 / 35.8, 3.9 | 6.6 / 17.5, 1.8 | 8.4 / 17.8, 1.8 |
| 0.3 | 0.999 (cap) | 29.7 / 57.3, 15.6 | 5.4 / 17.4, 1.8 | 7.3 / 18.2, 1.8 |
| 0.25 | 0.999 (cap) | 29.0 / 85.8, 33.7 | 4.8 / 17.7, 1.8 | 6.6 / 18.9, 1.8 |
| 0.15 | 0.999 (cap) | 25.7 / 111.9, 56.3 | 4.4 / 21.2, 1.8 | 5.5 / 23.5, 1.8 |

Means over 12 seeds; heading is mean RMS.

- **Legacy cap.** The known prior's coherence follows 0.334/k^2 (0.523, 0.685,
  0.936 at k = 0.8, 0.7, 0.6) and is pinned at the cap from k = 0.5. The
  shared-QMEKF heave change relative to real-only coupling goes from -5.5% at
  k = 1 and -16.5% at k = 0.6 to +62.4% [+55.8, +68.9] at k = 0.5, where the
  cap binds. Heading fails from k = 0.4 (QMEKF) and k = 0.3 (TFG). Raising
  only `S_factor` to 2 (`/sf2` rows) moves the observed boundary to between
  k = 0.3 and 0.25, as the 1/k^2 law predicts (k = 0.29).
- **Default policy.** In all 12 cells shared TFG has the lowest mean heave of
  the four filters. It has the lowest horizontal error in 9 and is within
  0.02% of the lowest in the other 3. Relative to independent TFG its heave
  changes by -27.9% to -5.7% and its horizontal error by -9.8% to -1.4% over
  the default-spreading factors. No shared run exceeds 30 degrees heading.

## Other conditions and limits

- All trajectories remain below the inherited position/nonfinite divergence
  flag, and no operational shared run exceeds 30 degrees RMS heading.
- Secondary endpoints are condition dependent. At 800 rpm, shared TFG
  lowers tilt RMS from 0.45 to 0.31 degrees, heading from 3.02 to 2.10 degrees
  and accelerometer-bias error from 7.84 to 5.40 mg relative to shared QMEKF.
  During maneuvers it raises heading from 1.75 to 2.28 degrees and gyro-bias
  RMS by 10.0% [6.9, 13.2]. These are exploratory; only the primary position
  tests use Holm correction.
- Nominal joint physical NEES is 39.9 for shared TFG and 40.5 for
  shared QMEKF, against a 9-state reference mean of 9, despite conservative
  marginal position coverage. Geometry does not fix calibration. One-second
  operational covariance audits remain symmetric, finite and PSD; the float
  stress test records a small rounding eigenvalue.
- The pinned-record heave geometry changes are below 0.013%, and transition
  recovery differences are practically small.

## Computation and hardware limits

| Measurement | QMEKF shared | TFG shared |
|---|---:|---:|
| Host core mean, median of 3 repeats | 3.981 us | 7.987 us |
| ESP32-S3 core object | 12,504 bytes | 14,472 bytes |
| Matched full AtomS3R application flash | 1,094,467 bytes | 1,029,687 bytes |
| Full application globals | 31,660 bytes | 31,660 bytes |

A separate four-repetition counterbalanced check rotates every core through
each execution position (`host_order_*.csv`). Its shared-QMEKF/TFG medians are
3.946/8.720 us at fixed intervals (ratio 2.21), and 6.733/11.660 us with
interval jitter (ratio 1.73). Host variation remains substantial.

These timings, the Xtensa probe and the firmware builds were measured before
the cap policy. The policy adds one float to the coupling configuration and a
few operations to the 0.1 s model commit; the measurements have not been
repeated. Core objects are allocated at construction, separately from firmware
globals. The Xtensa stack figures are compiler call-chain lower bounds, **not**
measured peak stack. The on-device 200 Hz probe compiles, but no connected
board was available: actual deadline, minimum heap and peak-stack measurements
remain unverified. Host timing cannot establish those properties.

## Files and reproduction

- `protocol.json`: primary/ablation seeds, cases, policy, windows and hashes.
- `runs.csv.gz`: compact full-rate per-run scores; no raw sample streams.
- `paired.csv`, `aggregate.csv`, `interactions.csv`, `ablations.csv`: paired
  statistics including negative and insignificant changes.
- `dataset_*`, `reset_*`, `transient_*`, `weak_*`, `heave_sweep_*`: diagnostic
  protocols and results. Pinned replay intervals condition on fixed wave
  records. `/cap0` rows are the legacy cap; `/cap0.7` and `/cap0.9` bracket
  the target.
- `execution_provenance.json`: the source revision, runner and bench binary of
  this execution.
- `regression_gates.json`, `review_validation.json`: `make all`, the
  validation suite, lint and timing-order checks with original thresholds.
- `shipping_faithfulness.json`: fresh OU-III finite diagnostic reproduces both
  original raw byte hashes; only its source binding and rounding-level summary
  values change. No theorem claim or shipping behavior is changed.
- `mathematical_verification.txt`, `refactor_regression.json`: independent
  27-state checks and exact original/post-refactor four-wrapper traces.
- `figures/` and `traces/`: bounded SVG/PGF figures and their compact drawing
  inputs. Scores and transient recovery use full-rate data before decimation.
- `environment.json`, `host_*.csv`, `esp32s3_size_stack.csv`,
  `operation_counts.json`, `firmware_builds.json`: computational evidence.

From the repository root, reanalyse the committed scores without simulations:

```sh
python3 tools/tfg_shared_modes_study.py analyse
python3 tools/tfg_shared_modes_heave_sweep.py analyse
python3 tools/tfg_shared_modes_report.py
(cd doc/kalman_tfg_shared_modes && latexmk -lualatex -halt-on-error kalman_tfg_shared_modes.tex)
```

Fresh experiments must use a fresh results directory to avoid mixing compiler,
source or binary hashes. Full commands, APIs and hardware instructions are in
[`docs/tfg-carma-shared-modes.md`](../../../docs/tfg-carma-shared-modes.md).
The standalone article is
[`kalman_tfg_shared_modes.tex`](../../../doc/kalman_tfg_shared_modes/kalman_tfg_shared_modes.tex).
