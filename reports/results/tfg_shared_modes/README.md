# TFG shared-mode CARMA: measured results

The reusable 27-state TFG implementation works, but the evidence does **not**
support a general replacement of shared QMEKF. Keep it experimental and keep
the shipping OU-III default. The main practical problem is quadrature coupling
under reduced heave, affecting both shared geometries.

| TFG shared versus QMEKF shared | Condition-balanced change | Paired 95% interval |
|---|---:|---:|
| Heave RMS | +0.499% | [+0.305, +0.694]% |
| Horizontal RMS | -1.433% | [-2.025, -0.838]% |

Each condition has equal log-ratio weight within each independent seed. The
primary study uses 24 seeds, 27 conditions, 1200 seconds at 200 Hz, scored from
300 seconds. There are 4,204 estimator trajectories including registered
ablations, reset diagnostics, pinned replays, transition diagnostics and the
exploratory reduced-heave study. See `experiment_counts.json` for the accounting.
The source is frozen for the simulations; subsequent runner changes only add
output routing and reanalysis from committed scores. Original source is
recoverable at the revision in `execution_provenance.json`.

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

## Computation and hardware limits

| Measurement | QMEKF shared | TFG shared |
|---|---:|---:|
| Host core mean, median of 3 repeats | 3.981 us | 7.987 us |
| ESP32-S3 core object | 12,504 bytes | 14,472 bytes |
| Matched full AtomS3R application flash | 1,094,467 bytes | 1,029,687 bytes |
| Full application globals | 31,660 bytes | 31,660 bytes |

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
- `dataset_*`, `reset_*`, `transient_*`, `weak_*`: diagnostic protocols and
  results. Pinned replay intervals condition on fixed wave records.
- `regression_gates.json`: all 19 original suite scripts completed, including
  1,152 validation tests, with original thresholds; execution repairs and
  resumed command groups are recorded explicitly.
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
