# Paired OU-II/OU-III validation

`tools/ou_validation.py` implements the paired inferential validation path for
the two OU filter families. It does not replace the executable regression gates:
the simulators still calculate and enforce their own thresholds over the
trailing 900 seconds of each 20-minute realization, which is the same window
the full experiment scores.

The repository versions the completed ten-seed study used by the manuscript in
`reports/results/ou_validation/`. It contains 840 simulator rows over ten
paired seed triplets: nine scenarios (four JONSWAP seas, four PM-Stokes seas,
one controlled transition) and two filter families. The three primary modes run
on all nine scenarios; the two covariance-policy controls and the two
OU-III-only channel-freeze modes run on the five JONSWAP-plus-transition
scenarios.

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
- That transition is a **crossfade between two independently phase-randomized
  stationary records**, not a continuously evolving spectrum. During the blend
  both spectra coexist, so there is no single intermediate $T_p$ and the sea can
  be bimodal. Because the records are independent their variances add, so the
  midpoint effective height is $\sqrt{0.5^2(1.5)^2+0.5^2(4.0)^2}\approx 2.14$ m
  rather than the 2.75 m of a linear $H_s$ ramp.
- The 900 s score therefore covers 300--420 s of the pure start sea, 420--780 s
  of blend, and 780--1200 s of the pure endpoint sea. 47% of the scored record
  is the endpoint sea alone, which favours the endpoint-calibrated fixed
  reference. `transition_window_composition_sec` in the manifest records the
  split for whatever window is configured.

The adaptation ablations are:

- `Adaptive`: the online tuner remains enabled.
- `FixedNominal`: parameters are held at the noise-free full-trace operating point
  for the nominal $H_s=1.5$ m, $T_p=5.7$ s sea, in every scenario.
- `FixedOracle`: parameters are held at the operating point calibrated from the
  stationary sea being scored; for the transition it uses the known final-sea
  point. This is a **scenario-calibrated fixed reference**, not a deployable
  online estimator and not an optimum.
- `AdaptiveHeldCovariance` and `FixedNominalHeldCovariance`: matched controls
  that repeat their partner mode with the periodic covariance re-alignment
  switched off (see below).

There is no separate "oracle" solver. All three of the first modes run the same
filter; the fixed modes simply freeze `(tau, sigma_aw, r_S)` after the normal
startup/Live transition. Each fixed triple is obtained by running the adaptive
filter once on a noise-free, unrandomized 1200 s record, reading its final
`tau` and `sigma_aw`, and computing
`r_S = clip(0.35 * sigma_aw * tau**3, 0.4, 400)` (OU-III), or, for OU-II,
`r_p0 = clip(0.6 * sigma_aw * tau**2, 0.05, 150)` together with
`r_v0 = clip(1.1 * sigma_aw * tau, 0.01, 40)`. Both families derive `tau` from
the same wave-band zero-crossing period, `tau = T_z / 2`; neither reads the
acceleration-band frequency tracker for tuning.
No fixed point is optimized against displacement error. The exact
frozen values for the committed study are in `fixed_tuning_points` in the
manifest and are typeset by `ou_validation_tuning_points.tex`.

Those values are the vertical/base parameters. The filter derives the applied
anisotropic values internally: OU-III uses
`(1.87*sigma_aw, 1.87*sigma_aw, sigma_aw)` for the stationary acceleration
standard deviation and `diag(r_S, r_S, r_S)**2` for the integral
pseudo-measurement covariance; OU-II uses `1.5*sigma_aw` horizontally and, since
the operating point moved to the wave band, an isotropic `r_p0`.

### Covariance policy and its control

The filter re-aligns the posterior marginal `P_aw_aw` with the stationary OU
covariance once per adaptation period, keeping the cross-covariances it has
learned. It stops the marginal settling far below the level the process model
considers stationary, which keeps the accelerometer gain responsive when the sea
state changes. The inflation is not small -- typically a factor of 40 to 100 --
and the operation is not a consistent posterior update: keeping the raw
cross-covariances while replacing the marginal rescales the implied correlation
coefficients by the square root of the marginal change.

`W3D_AW_COV_SYNC=congruent` performs the same re-alignment as a congruence,
which reaches the same marginal, leaves the whitened cross-covariance untouched
and stays positive semi-definite by construction. It is measurably *worse* than
the deployed overwrite, because consistency propagates the inflation into the
cross-covariances and the filter cannot absorb it. The conclusion is that the
periodic re-alignment should be retired or bounded rather than made
self-consistent; the `*HeldCovariance` modes are what price that.

Earlier revisions applied the re-alignment inside the adaptation path only.
That made it run in `Adaptive` and never in a fixed mode that stops re-tuning,
so comparing the two confounded *whether parameters adapt* with *whether part
of the covariance is periodically re-aligned*. The re-alignment is now driven
independently of the tuner, so all three primary modes apply it identically.
The `*HeldCovariance` modes switch it off at matched tuning settings so the
policy can be priced on its own. Use `W3D_AW_COV_SYNC=reconfigure` to run that
policy directly; `periodic` is the default.

## Running

Fetch the versioned simulation data and run a short integration check:

```bash
make fetch-sim-data
python3 tools/ou_validation.py --mode smoke
```

Run the declared full defaults:

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
`--skip-nonstationary`, `--duration-sec`, `--window-sec`, `--jobs`, and
`--bootstrap-resamples`. `--mode smoke` is only an integration test: its single
short realization is not statistical evidence.

`--mode full` is roughly 840 twenty-minute simulator replays and takes hours,
so the committed bundle is regenerated by the dispatch-only `regenerate` job in
`.github/workflows/ou-validation.yml` rather than on a workstation or in the
pull-request gate. That job runs the validation and robustness bundles as
separate matrix legs, because the two together do not fit in one runner
lifetime. It uploads the result as an artifact rather than committing it: a
regenerated bundle has to be read against the manuscript before it replaces the
committed one.

### Restating a bundle without re-simulating

The replays are the expensive half of this study and are fully determined by
their seed triplets, so a change to how the rows are *summarized* does not
require re-running them:

```bash
python3 tools/ou_validation.py \
  --restat-from reports/results/ou_validation/ou_validation.json \
  --output-dir reports/results/ou_validation
```

This reads `raw_runs` back out of the bundle and rewrites every derived file --
summaries, paired effects, LaTeX tables, generated macros, bundle, and manifest
-- with the statistics of the current source. Restating the committed bundle
unchanged reproduces its derived files byte for byte, which
`tests/validation/test_ou_validation.py` asserts; that is what makes it safe to
use for adding or changing an interval construction. It cannot invent rows:
whatever ensemble the source bundle scored is the ensemble the restated bundle
reports, so widening the seed set still requires a `--mode full` run.

Protocol fields that describe what was *run* (seeds, durations, transition
bounds) are carried over untouched. Fields that describe how the rows are
summarized are taken from the current source, so a restated bundle cannot
disagree with the manuscript generated beside it.

## Outputs and interpretation

The output directory contains:

- `ou_validation_raw.csv`: one row per paired filter/mode run, including the
  historical gate result;
- `ou_validation_summary.csv`: $n$, mean, sample standard deviation, normal 95%
  interval, and bootstrap 95% interval;
- `ou_validation_paired_effects.csv`: paired mean differences, bootstrap intervals,
  Cohen's $d_z$, and small-sample-corrected Hedges' $g_z$;
- the `stationary_normalized_aggregate` block of the bundle and manifest, which
  carries the declared primary endpoint under four constructions on the same
  paired differences: the percentile bootstrap, a Student-t interval, an exact
  sign test, and an exact paired randomization (sign-flip) test enumerated over
  all $2^n$ sign patterns. Descriptive contrasts keep the bootstrap alone; four
  p-values on every one of them would enlarge the family of tests without
  adding evidence. At $n=10$ both exact tests bottom out at $2^{-9}=0.002$
  regardless of the data, so that floor is a property of the design;
- `ou_validation.json`: protocol, calibration points, raw observations, and all
  statistics;
- `ou_validation_manifest.json`: command, versions, Git state, source-file hashes,
  result-file hashes, seed triplets, and the stationary normalized aggregate;
- `ou_validation_tuning_points.tex`: the exact frozen operating points used by
  the fixed modes;
- `ou_validation_transition.csv` and `ou_validation_transition.svg`: a decimated
  time series of one transition realization -- blend weight, reference rolling
  $H_s$, reference and estimated vertical displacement, error, and the applied
  `tau`/`sigma_aw`/`r_S` against the two fixed levels;
- `ou_validation_table.tex`, `ou_validation_publication.tex`, and the SVG
  figures for publication workflows.

Differences are always `left - right`. Thus negative displacement or attitude
error differences favor the left-hand estimator. A confidence interval spanning
zero is inconclusive at that interval level. Effect sizes should be interpreted
with their paired sample count and interval, not in isolation.

Across the four stationary JONSWAP seas, the within-seed mean normalized
vertical RMS is 4.98% of $H_s$ for OU-III and 6.68% for OU-II. The paired
difference is -1.697 percentage points with a bootstrap 95% interval of
[-1.781, -1.623], a Student-t interval of [-1.794, -1.601], and all ten
seed-level differences negative (exact sign and paired sign-flip tests both at
$p=0.002$, the smallest two-sided p-value ten pairs can produce). OU-III is
lower vertically in every stationary sea and during
the transition, and has higher 3D displacement RMS in four of the five
scenarios. In the fifth, the $H_s=8.5$ m sea, it is lower with an interval
excluding zero; that is a scenario-specific exception, not a general
three-dimensional advantage. Which band and which axis produce that pattern is
attributed in [`ou-3d-error-attribution.md`](ou-3d-error-attribution.md).

Both families now take the operating point from the same wave-band
zero-crossing period. Earlier bundles tuned only OU-III that way, which made
the reported vertical difference (-3.233 points) roughly twice what the
architecture alone accounts for, and made the 3D count look like 1 of 5 rather
than 4 of 5. Conclusions should remain conditional on the evaluated JONSWAP
family rather than be read as uniform estimator dominance.

The follow-on OU-III parameter-sensitivity and degradation-case protocol is
documented separately in [`ou-robustness.md`](ou-robustness.md). It reuses the
same paired seed triplets but does not alter the confirmatory comparison or
select a new reported operating point from reference errors.

The simulator pass/fail thresholds are intentionally retained as
deterministic regression sentinels. They are calibrated to the deterministic
realization: only a minority of the kinematically projected surrogate runs
satisfy all of them, and the exact count is regenerated into
`\OUValidationGatePasses`. They
must not be interpreted as ensemble acceptance criteria; the statistical study
uses the raw long-window metrics and paired intervals instead.

## Direct simulator controls

Both OU executables accept repeated `--input PATH` arguments. Their validation
environment variables are:

- `W3D_IMU_SEED`, `W3D_INIT_SEED`, or combined `W3D_SEED`;
- `W3D_VALIDATION_WINDOW_SEC` for the machine-readable `VALIDATION_METRICS` line;
- `W3D_WRITE_TIMESERIES=0` to suppress large diagnostic CSVs;
- `W3D_TUNING_MODE=adaptive`, `fixed_nominal`, or `fixed_oracle`, plus the
  family-specific `W3D_FIXED_*` values;
- `W3D_AW_COV_SYNC=periodic` (default) or `reconfigure` for the a_w
  covariance-alignment policy.

If no seed variable is present, the original deterministic 1234/5678/9012
realization and its historical random-number draw order are retained.
