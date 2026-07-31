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

- Each repetition has three independent seeds: wave phase, IMU noise/random walk,
  and sensor-error initialization.
- A repetition uses the identical seed triplet for OU-II and OU-III and for every
  adaptation ablation, enabling paired inference.
- Wave realizations use one random Fourier phase per frequency, shared by all
  world-motion and Euler channels. This preserves auto-spectra and cross-spectra.
  Quaternion, body specific force, and body angular rate are then reconstructed.
- Full mode defaults to ten predeclared seed triplets. Seeds can be supplied as
  one value (broadcast) or equally sized comma-separated lists.
- The non-stationary case transitions smoothly from
  $H_s=1.5$ m, $T_p=5.7$ s to $H_s=4.0$ m, $T_p=11.4$ s between 420 and 780 seconds.

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

The completed study does not show uniform dominance. Averaged within each seed
across the four stationary JONSWAP seas, normalized vertical RMS is 8.73% of
$H_s$ for OU-III and 8.98% for OU-II, but the scenario-level difference changes
sign and OU-III has higher 3D displacement RMS in every evaluated scenario.
Adaptive tuning improves on FixedNominal at the smallest and largest stationary
seas, but it is worse than FixedNominal during the controlled transition. The
manuscript tables report the paired bootstrap intervals and effect sizes.

The historical final-60-second pass/fail thresholds are intentionally retained
as deterministic regression sentinels. They were calibrated to the original
realization: only 53 of the 300 phase-surrogate runs satisfy all of them. They
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
