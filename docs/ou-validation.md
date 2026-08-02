# Paired OU-II/OU-III validation

`tools/ou_validation.py` implements the statistically powered validation path for
the two OU filter families. It does not replace the executable regression gates:
the simulators still calculate and enforce their historical trailing 60-second
thresholds. The full experiment separately scores the trailing 900 seconds of each
20-minute realization.

The repository versions the completed ten-seed study used by the manuscript in
`reports/results/ou_validation/`. It contains 300 simulator rows: five scenarios,
two filter families, three tuning modes, and ten paired seed triplets.

## Experiment design

- Each repetition has three independent seeds: wave realization, IMU noise/random walk,
  and sensor-error initialization.
- A repetition uses the identical seed triplet for OU-II and OU-III and for every
  adaptation ablation, enabling paired inference.
- Wave realizations use one random Fourier phase per retained frequency,
  shared by world velocity and Euler channels over the JONSWAP model's
  0.02--1.60 Hz fundamental-plus-bound-harmonic band. Displacement and
  acceleration are analytically derived from the same randomized velocity
  spectrum, eliminating finite-record boundary leakage while preserving the
  retained primitive auto- and cross-spectra. Quaternion, body specific force,
  and body angular rate are then reconstructed.
- Full mode defaults to ten predeclared seed triplets. Seeds can be supplied as
  one value (broadcast) or equally sized comma-separated lists.
- The non-stationary case uses a C2 quintic transition from
  $H_s=1.5$ m, $T_p=5.7$ s to $H_s=4.0$ m, $T_p=11.4$ s between 420 and 780 seconds.
  Velocity and acceleration include the exact first- and second-derivative
  terms introduced by the time-varying blend.

The adaptation ablations are:

- `Adaptive`: the online tuner remains enabled.
- `FixedNominal`: parameters are held at the noise-free full-trace operating point
  for the nominal $H_s=1.5$ m, $T_p=5.7$ s sea.
- `FixedOracle`: parameters are held at the known stationary-sea operating point;
  for the transition it uses the known final-sea operating point. This is a
  clairvoyant fixed endpoint baseline, not a deployable online estimator.

## Running

Fetch the versioned simulation data and run a short integration check:

```bash
make fetch-sim-data
python3 tools/ou_validation.py --mode smoke
```

Run the preregistered full defaults:

```bash
python3 tools/ou_validation.py --mode full
```

Configure independent seed streams explicitly:

```bash
python3 tools/ou_validation.py --mode full \
  --wave-seeds 11,29,47,71 \
  --imu-seeds 101,211,307,401 \
  --initialization-seeds 1009,1103,1201,1301
```

Useful controls include `--stationary-input` (repeatable),
`--skip-nonstationary`, `--duration-sec`, `--window-sec`, and
`--bootstrap-resamples`. `--mode smoke` is only an integration test: its single
short realization is not statistical evidence.

## Outputs and interpretation

The output directory contains:

- `ou_validation_raw.csv`: one row per paired filter/mode run, including the
  historical gate result;
- `ou_validation_summary.csv`: $n$, mean, sample standard deviation, normal 95%
  interval, and bootstrap 95% interval;
- `ou_validation_paired_effects.csv`: paired mean differences, bootstrap intervals,
  Cohen's $d_z$, and small-sample-corrected Hedges' $g_z$;
- `ou_validation.json`: protocol, calibration points, raw observations, and all
  statistics;
- `ou_validation_manifest.json`: command, versions, Git state, source-file hashes,
  result-file hashes, seed triplets, and the stationary normalized aggregate;
- `ou_validation_table.tex`, `ou_validation_publication.tex`, and three SVG
  figures for publication workflows.

Differences are always `left - right`. Thus negative displacement or attitude
error differences favor the left-hand estimator. A confidence interval spanning
zero is inconclusive at that interval level. Effect sizes should be interpreted
with their paired sample count and interval, not in isolation.

Across the four stationary JONSWAP seas, the within-seed mean normalized
vertical RMS is 7.22% of $H_s$ for OU-III and 8.31% for OU-II. The paired
difference is -1.086 percentage points with a bootstrap 95% interval of
[-1.208, -0.970]. OU-III is lower vertically in every stationary sea, but is
higher during the transition and has higher 3D displacement RMS in four of the
five scenarios; the nominal-sea 3D interval spans zero. Conclusions should
therefore remain conditional on the evaluated JONSWAP family rather than be
read as uniform estimator dominance.

The follow-on OU-III parameter-sensitivity and degradation-case protocol is
documented separately in [`ou-robustness.md`](ou-robustness.md). It reuses the
same paired seed triplets but does not alter the confirmatory comparison or
select a new reported operating point from reference errors.

The historical final-60-second pass/fail thresholds are intentionally retained
as deterministic regression sentinels. They were calibrated to the original
realization: 61 of the 300 kinematically projected surrogate runs satisfy all
of them. They
must not be interpreted as ensemble acceptance criteria; the statistical study
uses the raw long-window metrics and paired intervals instead.

## Direct simulator controls

Both OU executables accept repeated `--input PATH` arguments. Their validation
environment variables are:

- `W3D_IMU_SEED`, `W3D_INIT_SEED`, or combined `W3D_SEED`;
- `W3D_VALIDATION_WINDOW_SEC` for the machine-readable `VALIDATION_METRICS` line;
- `W3D_WRITE_TIMESERIES=0` to suppress large diagnostic CSVs;
- `W3D_TUNING_MODE=adaptive`, `fixed_nominal`, or `fixed_oracle`, plus the
  family-specific `W3D_FIXED_*` values.

If no seed variable is present, the original deterministic 1234/5678/9012
realization and its historical random-number draw order are retained.
