# Conventional heave baselines: HPDI and FDDI

Two conventional comparators for OU–III's vertical channel. Both consume the
**same measurement-only inputs** OU–III's front end produces: the Mahony-proxy
vertical acceleration `a_up` (paper Eq. 73) and, for period-scaled cutoffs, the
canonical zero-crossing period `Tz_hat` (Eq. 80). Reference motion is used only
to score results and to select frozen parameters on a development split.

| Method | Causal | Where | Transfer from acceleration to heave |
|---|---|---|---|
| HPDI-classic | yes | `src/hpdi/HeaveHPDI.h` | `s^(n-2) / D(s)`, i.e. Butterworth HP of order n then 1/s² (m = n) |
| HPDI-compensated | yes | same, m < n | `[1 − R(s)/D(s)] / s²` |
| FDDI | no (offline) | `tools/heave_baselines.py` | `−W(f) / (2πf)²`, raised-cosine high-pass taper |

## Design summary

`D(s)` is the order-n Butterworth polynomial at cutoff ω_c, and `R(s)` is its
lowest m terms. The displacement gain relative to ideal double integration is
`G = 1 − R/D`.

- **Classic (m = n).** `G = s^n/D`. The error `|G − 1|` decays only as
  `c_{n−1} ω_c/ω`. This is the persistent phase lead of a causal high-pass, and
  is why textbook heave filters distort the wave band.
- **Compensated (m < n).** `|G − 1| ~ c_{m−1} (ω_c/ω)^(n−m+1)`. The cost is gain
  peaking near the corner (|G| up to ~2 for n = 4, m = 3; ~6 for n = 6, m = 3).
- **m ≥ 3.** Zero steady-state heave for a constant accelerometer bias, and
  bounded variance under a random-walk bias.
- **m = 2.** Leaves an offset of `c_2 b / ω_c²`. The member n = 3, m = 2 with
  ω_c = ω_R reproduces the reduced SpectralMSE response, paper Eq. (86), so the
  HPDI family contains OU–III's own scalar reduction.

**Discretization.** The trapezoidal rule (bilinear transform) is applied in delta
form on an observer-canonical realization with frequency-normalized states. This
stays float-accurate at ω_c·h ≈ 10⁻³. `selfcheck` verifies agreement with
`scipy.signal.bilinear` of H(s) to about 1e-8.

**Embedded footprint.** `HeaveHPDI<float, 6>`:
- 472 B of state, no heap.
- About 3 KB of code at `-Os` with `-fno-exceptions -fno-rtti`.
- Per sample: n² + 2n MACs plus one `expm1f`.
- Coefficient refresh about every 0.1 s (sample-and-hold, like OU–III).

## Build and test

```sh
make -C tests/hpdi            # unit test (67 checks), bilinear and scorer self-checks
```

`tests/hpdi` is part of `make all` and needs no simulation data. CMake builds
`test_heave_hpdi` and `hpdi_replay`, and the CI `cmake` job runs both checks.
The header compiles cleanly from C++11 through C++20 under `-Wconversion
-Wdouble-promotion`.

## Scoring

Every method is scored with the simulation harness's own scorer
(`util/W3dSimCommon.cpp`), so baseline and estimator numbers are directly
comparable:

- window: the last `floor(float(900) / float(1/200))` = 180000 samples; a
  segment `[t0, t1)` covers samples `floor(float(t0)/float(dt))` to
  `floor(float(t1)/float(dt))`;
- error: `float(estimate) - float(z_ref)`;
- RMS: squares accumulated in order into a float with a fused multiply-add
  (the harness compiles to `vfmadd231ss` under `-march=native`), then
  `sqrt(sum / N)` in float;
- `%Hs = rms * (100 / Hs)` in float, with `Hs` the record's nominal height.

`hpdi_replay --harness-window 900` and `heave_baselines.py` both implement it;
`selfcheck` tests the Python FMA emulation against exact rationals and against
`hpdi_replay`.

## Harness export

`tests/hpdi/heave_export-<family>` (OU-III, OU-II, TFG, PII, TVG-NLO) replays
one record through the shipped estimator and writes its per-sample vertical
error. Each build compiles that family's simulator source unchanged and runs its
own `main()`; two link-time wraps observe it (see `heave_export.cpp`). No
simulator, estimator or evidence-closure file changes, and the OU-III, OU-II
and TFG builds print `VALIDATION_METRICS` lines byte-identical to the shipped
simulators'. The OU-III build also exports, from that same filter's front end
after each sample:

| Column | Content |
|---|---|
| `a_up` | Mahony-proxy up-positive vertical acceleration, gravity removed (Eq. 73), after the vibration guard OU–III sees |
| `Tz_hat` | canonical log-period estimate (Eq. 80); NaN before its startup-usable gate |
| `z_ref` | the record's `disp_z`, which the harness scores as up-positive |

The study driver merges the five exports into one CSV per record:
`t, a_up, Tz_hat, z_ref, err_<F>, zhat_<F>` for each family `F`. `err_<F>` is
the harness's float error; `zhat_<F>` is `err_<F> + z_ref` in double.

## Study

```sh
python3 tools/heave_baseline_study.py run --work runs/heave_baselines
python3 tools/heave_baseline_study.py analyze --work runs/heave_baselines
```

`run` builds the binaries, replays every record, verifies the export, tunes,
freezes and evaluates; `analyze` rebuilds the statistics, `summary.md` and the
plot from the work directory. Results go to `reports/results/heave_baselines/`.

| Split | Records | Draws |
|---|---|---|
| eval | 8 stationary RAO records, controlled crossfade, low–high–low; Hs = 0.05 m and 30 s ramp (secondary) | the ten predeclared seed triplets of `tools/ou_validation.py`, phase-randomized as there |
| eval | 8 pinned records: default draw, noise-free (`--no-noise`), nominal-cruise engine vibration | simulator default |
| dev | 8 pinned records | comparator-retuning draws `default, 11, 23, 61001, 62003` (`W3D_IMU_SEED = W3D_INIT_SEED`) |

The dev wave histories (the pinned records, not phase-randomized) and sensor
draws are disjoint from eval's; they share the hull and the incident spectra,
as the comparator retuning did. The default draw is also a dev draw, so the
baselines' default-draw (Table X layout) entries are in-sample.

The driver checks, and fails otherwise, that:
- the Python rescoring of every exported error series equals the harness's
  printed %Hs in float32, for the whole window and every segment;
- the export builds' metric lines equal the shipped simulators' byte for byte;
- the CLI `evaluate` path and the driver score the frozen baselines identically;
- no tuned optimum lies on a grid edge.

It also records how the replayed metrics compare with the committed evidence
bundles, the sign of `z_ref` against `a_up`, and where the baseline code reads
reference motion.

## Standalone use

```sh
python3 tools/heave_baselines.py tune --manifest runs/manifest.csv --out runs/hpdi_fddi_frozen.json
python3 tools/heave_baselines.py evaluate --manifest runs/manifest.csv \
    --params runs/hpdi_fddi_frozen.json --out runs/hpdi_fddi_results.csv \
    --compare-csv runs/ou3_results.csv --compare-method OU-III
python3 tools/heave_baselines.py rescore --manifest runs/manifest.csv --out runs/rescored.csv
```

The manifest lists `path,family,hs,split,seed` (extra columns are ignored);
`split=dev` records select the frozen parameters and `split=eval` records score
them. `tune` grid-searches FDDI (fixed or period-scaled cutoff, crossed with the
taper ratio) and HPDI (`3 <= m <= n <= 6`, fixed or period-scaled cutoff), and
freezes FDDI, the best `m = n` (classic) and the best `m < n` (compensated)
HPDI; it flags an optimum on a grid edge in the frozen file. `evaluate` writes
long-format per-record rows and, with `--compare-csv` (columns
`method,family,hs,seed,pct_hs`), a seed-level paired contrast with a 10,000
resample percentile bootstrap (PCG64 seed 20260317) and an exact sign-flip test.
`--float32` replays HPDI in single precision, as deployed on the ESP32-S3.
`--flip-z` negates a NED-down reference.

## Reporting notes

- FDDI is zero-phase and uses future samples. Report it as an offline reference
  bound, not a real-time competitor.
- Both baselines estimate heave only. Compare them on the vertical endpoint, and
  keep OU–III's 3-D and attitude outputs as a separate claim.
- `synth` writes idealized vertical-only records: no RAO, no tilt, no proxy, and
  the true T_z. Use them for plumbing tests only; they are not evidence.

## Results

`reports/results/heave_baselines/summary.md` holds every table, contrast and
check; `heave_baselines_long.csv` holds one row per method, record and scoring
interval; `frozen.json` the selected parameters; `heave_baselines_vertical.svg`
the article's Fig. 6 layout with the baselines added. On the primary-style
stationary JONSWAP endpoint the causal HPDI forms are worse than OU-III and the
offline FDDI reference is better; see the summary for the intervals.
