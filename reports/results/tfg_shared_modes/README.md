# TFG shared-mode CARMA: measured results

The reusable 27-state TFG implementation works, but the evidence does **not**
support a general replacement of shared QMEKF. Keep it experimental and keep
the shipping OU-III default. The main practical problem affects both shared
geometries: under reduced heave the shared excitation mapping caps the
quadrature coupling into a singular prior (see the dose response below).

| TFG shared versus QMEKF shared | Condition-balanced change | Paired 95% interval |
|---|---:|---:|
| Heave RMS | +0.499% | [+0.305, +0.694]% |
| Horizontal RMS | -1.433% | [-2.025, -0.838]% |

Each condition has equal log-ratio weight within each independent seed. The
primary study uses 24 seeds, 27 conditions, 1200 seconds at 200 Hz, scored from
300 seconds. There are 5,920 estimator trajectories including registered
ablations, reset diagnostics, pinned replays, transition diagnostics and the
two exploratory reduced-heave studies. See `experiment_counts.json` for the accounting.
The registered simulation source and original protocols remain recoverable at
the revision in `execution_provenance.json`. Subsequent constructor lint repairs
retain the same expressions; repeated equivalence checks and counterbalanced
host timing are recorded in `review_validation.json`. Runner updates make replay
row ordering deterministic and derive report claims from the measured tables.

## Failure and condition dependence

- Reduced heave: shared QMEKF averages 21.8 mm heave / 163.2 mm horizontal and
  95.2 degrees heading RMS; shared TFG averages 25.0 / 131.4 mm and 70.9 degrees.
  Independent TFG averages 5.4 / 23.4 mm and 1.9 degrees. RMS heading exceeds
  30 degrees in 24/24 shared-QMEKF and 17/24 shared-TFG primary runs.
- All trajectories remain below the inherited position/nonfinite divergence
  flag. That flag **does not** imply useful estimation: the heading losses and
  enormous joint NEES above are practical failures.
- In the separate 12-seed exploratory diagnosis, real-only and uncoupled
  controls remove the large heading failure. Known quadrature coupling still
  fails; confidence and its PSD fallback only partly contain the problem.
  This implicates prior/attitude/bias interaction rather than establishing a
  unique cause or validating a replacement policy.
- The exploratory heave-RAO sweep (below) traces this failure to the shared
  excitation mapping, not to TFG geometry.
- Ordinary position improvements are mostly tiny. At 800 rpm the geometry
  improves tilt, heading and accelerometer-bias tracking, but maneuvering
  worsens heading (1.75 to 2.28 degrees) and gyro-bias tracking. These are
  exploratory secondary endpoints; only primary position tests use Holm
  correction. See all metrics and paired intervals in the CSV files.
- Nominal joint physical NEES is about 40 against a 9-state reference mean of
  9, despite conservative marginal position coverage. Geometry does not fix
  calibration. One-second operational covariance audits remain symmetric,
  finite and PSD; the float stress test records a small rounding eigenvalue.
- The pinned-record heave geometry changes are below 0.013%; transition
  recovery differences are practically small. Reset transport is retained as
  the verified default even though omitting it changes the reduced-heave
  failure mode.

## Reduced-heave dose response and mechanism

`tools/tfg_shared_modes_heave_sweep.py` scales the heave RAO by `k` on the
first 12 primary seeds; horizontal motion and sensor streams are unchanged.
The repeated `k=1` and `k=0.15` cells reproduce all 132 overlapping published
runs exactly, and a rebuilt binary reproduces every first-run row.

| k | Known-prior coherence | QMEKF shared heave vs real-only | Heading > 30 deg (Q/T) |
|---:|---:|---:|---:|
| 1 | 0.334 | -5.5% [-5.8, -5.1] | 0/0 of 12 |
| 0.6 | 0.936 | -16.5% [-17.6, -15.5] | 0/0 |
| 0.5 | 0.999 (cap) | +62.5% [+55.9, +69.0] | 0/0 |
| 0.3 | 0.999 (cap) | +257% [+234, +280] | 11/0 |
| 0.25 | 0.999 (cap) | +255% [+229, +280] | 11/7 |

- The wrapper's horizontal prior is `S_factor` times the vertical one
  (default 1). The excitation mapping scales the measured horizontal-on-
  vertical regression by the vertical prior variance, so the applied multiple
  coherence grows as 0.334/k^2 (measured 0.523, 0.685, 0.936 at k=0.8, 0.7,
  0.6). From k=0.5 it is pinned at the 0.999 positive-semidefinite cap. The
  capped prior is singular: no horizontal acceleration remains independent
  of heave.
- Up to the cap the quadrature prior helps more as coherence rises. At the cap
  it turns into heave harm (+62% at k=0.5), while mean coupling confidence is
  still 0.96. Heading fails later: k=0.4 for QMEKF and k=0.3 for TFG.
- Intervention: raising only `S_factor` to 2 (`CSM/sf2`, `TCSM/sf2`) divides
  the coherence by four and predicts a boundary at k=0.29. Observed:
  -22.0% [-23.3, -20.8] at k=0.3 and +50.1% [+40.6, +59.7] at k=0.25, with no
  heading failure in any intervention cell. At k=1 the oversized prior cuts
  the shared gain to -1.5%.
- TFG geometry delays but does not prevent the collapse. At k=0.3 it changes
  heading by -48.5 [-65.5, -31.4] degrees and horizontal error by -70 mm
  relative to shared QMEKF, but heave by +5.3 mm. Both stay far worse than the
  independent-axis filters, whose heading is 1.8 degrees at every factor.
- At k=0.15, broad spreading reduces and long-crested seas increase the
  shared-QMEKF heading error. The estimated-coupling TFG intervals include
  zero.

Implication: the reduced-heave failure belongs to the shared excitation
mapping used by both geometries. Further shared-mode work should first
replace cross-term shrinkage at the cap by a horizontal prior consistent with
the measured amplitude ratio, and make cap binding visible to the confidence.
It should then re-register the complete study, including calm water, where
horizontal band power is dominated by leakage. This sweep retunes nothing.

## Computation and hardware limits

| Measurement | QMEKF shared | TFG shared |
|---|---:|---:|
| Host core mean, median of 3 repeats | 3.981 us | 7.987 us |
| ESP32-S3 core object | 12,504 bytes | 14,472 bytes |
| Matched full AtomS3R application flash | 1,094,467 bytes | 1,029,687 bytes |
| Full application globals | 31,660 bytes | 31,660 bytes |

A separate four-repetition counterbalanced check rotates every core through
each execution position (`host_order_*.csv`). Its shared-QMEKF/TFG medians are
3.946/8.720 us at fixed intervals (ratio 2.21), and 6.733/11.660 us with interval
jitter (ratio 1.73). Host variation remains substantial; this check retains the
original publication timings and supports the same additional-cost conclusion.

Core objects are allocated at construction, separately from firmware globals.
The Xtensa stack figures are compiler call-chain lower bounds, **not** measured
peak stack. The on-device 200 Hz probe compiles, but no connected board was
available: actual deadline, minimum heap and peak-stack measurements remain
unverified. Host timing cannot establish those properties. Exact structured
kernel operation counts, interval-jitter timings, stack frames and linked
sizes are recorded separately.

## Files and reproduction

- `protocol.json`: primary/ablation seeds, cases, policy, windows and hashes.
- `runs.csv.gz`: compact full-rate per-run scores; no raw sample streams.
- `paired.csv`, `aggregate.csv`, `interactions.csv`, `ablations.csv`: paired
  statistics including negative and insignificant changes.
- `dataset_*`, `reset_*`, `transient_*`, `weak_*`, `heave_sweep_*`: diagnostic
  protocols and results. Pinned replay intervals condition on fixed wave
  records.
- `regression_gates.json`: uninterrupted `make all` exits 0 on a pristine PR
  head and on the follow-up branch (all 19 suite scripts, 1,152 validation
  tests, original thresholds); earlier workspace failures are recorded.
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
python3 tools/tfg_shared_modes_report.py
(cd doc/kalman_tfg_shared_modes && pdflatex -halt-on-error kalman_tfg_shared_modes.tex)
(cd doc/kalman_tfg_shared_modes && pdflatex -halt-on-error kalman_tfg_shared_modes.tex)
```

Fresh experiments must use a fresh results directory to avoid mixing compiler,
source or binary hashes. Full commands, APIs and hardware instructions are in
[`docs/tfg-carma-shared-modes.md`](../../../docs/tfg-carma-shared-modes.md).
The standalone article is
[`kalman_tfg_shared_modes.tex`](../../../doc/kalman_tfg_shared_modes/kalman_tfg_shared_modes.tex).
