# Paired OU-II/OU-III validation

`tools/ou_validation.py` implements the paired inferential validation path for
the OU-II and OU-III filter families. The committed publication evidence is in
`reports/results/ou_validation/`.

The current full bundle was regenerated for the shortened controlled transition
and the refitted `r_S` averaging horizon. Its immutable replay source commit is
`9276aee2c05804d61079d383d8657daa894c16b3`. It was produced outside GitHub
Actions, so the manifest's `workflow_run` is null and its `build_environment`
records the toolchain that did produce it; the replay hashes, not the runner,
are what the contract checks. The manifest separates
that replay provenance from any later statistical restatement. Replay provenance
pins the simulator/filter dependency closure, build files, versioned inputs, and
normalized raw-row hash; a later restatement records its analysis context in a
separate block and cannot replace those replay pins.

The completed study contains 840 simulator rows over ten paired seed triplets:
nine scenarios (four JONSWAP seas, four PM--Stokes seas, and one controlled
transition) and two filter families. The three primary modes run on all nine
scenarios; the covariance-policy controls and OU-III-only channel-freeze modes
run on the five JONSWAP-plus-transition scenarios.

## Experiment design

- Each repetition has three independent seeds: wave realization, IMU noise/random
  walk, and sensor-error initialization.
- A repetition uses the identical seed triplet for OU-II and OU-III and for every
  adaptation ablation, enabling paired inference.
- Wave realizations use one random Fourier phase per retained frequency, shared by
  world velocity and Euler channels. Displacement and acceleration are derived
  analytically from the same randomized velocity spectrum before quaternion,
  body specific force, and body angular rate are rebuilt.
- Full mode uses ten predeclared seed triplets. A single supplied seed is
  broadcast; equally sized lists remain paired element by element.
- The non-stationary case uses a C2 quintic crossfade from the $H_s=1.5$ m,
  $T_p=5.7$ s sea to the $H_s=4.0$ m, $T_p=11.4$ s sea between 540 and 660 s.
  The two endpoint records are independently phase-randomized, so the blend is
  deliberately bimodal rather than a continuously evolving single spectrum.
  The crossfade used to span 420--780 s. A 360 s ramp is quasi-static against
  the adaptation horizons it is meant to exercise -- the slowest is the
  two-stage $\sigma_a$ EWMA at roughly $2K T_z$, about 34 s on the largest
  reference sea -- so the score could not separate tracking lag from
  steady-state accuracy. At 120 s the crossfade is about three times that
  memory: still a sea state changing rather than a step, but fast enough for a
  lagging schedule to show up as an error.
- The final 900 s are scored. For the transition that window contains 300--540 s
  of the start sea, 540--660 s of blend, 660--780 s of endpoint-sea run-on, and
  780--1200 s of settled endpoint sea. The run-on interval is scored apart from
  the settled one because a schedule that averages too long does not only lag
  during the blend, it carries the old sea into the new one, and that cost
  lands after the crossfade; pooled into a single endpoint interval the two
  partly cancel. Shortening the crossfade from the front leaves the settled
  interval where it was, so settled numbers remain comparable with those the
  360 s protocol published.

The principal adaptation modes are:

- `Adaptive`: online sea-state tuning is enabled.
- `FixedNominal`: tuning is frozen at the noise-free full-trace operating point
  of the nominal $H_s=1.5$ m sea.
- `FixedOracle`: tuning is frozen at the noise-free operating point of the sea
  being scored; for the transition it uses the known endpoint. This is a
  scenario-calibrated reference, not a deployable online estimator and not an
  error-minimizing oracle.
- `AdaptiveHeldCovariance` and `FixedNominalHeldCovariance`: matched controls with
  periodic $P_{a_wa_w}$ re-alignment disabled.
- `AdaptiveRSOnly` and `AdaptiveOUOnly`: OU-III-only controls that freeze one of
  the two adaptation channels.

Both OU families derive their operating time scale from the same wave-band
zero-crossing period. The exact frozen points used by the committed study are in
`fixed_tuning_points` in the manifest and in
`ou_validation_tuning_points.tex`.

For current OU-III the translational OU acceleration prior and the integral
pseudo-measurement are isotropic: the stationary acceleration standard deviation
is `(sigma_aw, sigma_aw, sigma_aw)` and the integral pseudo-measurement covariance
uses the same `r_S` standard-deviation scale on X, Y, and Z. Thus the committed
Monte Carlo rows and the current estimator configuration now refer to the same
implementation.

## Statistical evidence versus simulator regression gates

The Monte Carlo evidence is a continuous-metric study, not a pass/fail study.
Every completed replay with a valid machine-readable metrics record is included
in the paired analysis. The statistical CSV/JSON bundle intentionally does not
export `quality_gate_pass`, `simulator_return_code`, or pass/fail counts.

The simulators may still retain deterministic threshold checks as executable
regression diagnostics for their canonical deterministic realization. Those
thresholds are **not** cohort acceptance criteria for randomized Monte Carlo
replays. Keeping that distinction in the schema prevents a reader from
mistaking a legacy regression threshold miss for a failed statistical
experiment.

The bundle protocol records this explicitly as:

- `simulator_regression_gates_exported: false`;
- `replay_inclusion_rule`: all completed replays with machine-readable metrics;
  deterministic simulator regression thresholds are not statistical acceptance
  criteria.

## Current full-study result

For the declared primary endpoint, the four stationary JONSWAP seas are first
averaged within each seed triplet and the paired seed-level values are then
analyzed. The current isotropic bundle gives:

- OU-II Adaptive: **6.464% of $H_s$** mean normalized vertical RMS;
- OU-III Adaptive: **4.729% of $H_s$**;
- OU-III minus OU-II: **-1.735 percentage points**;
- percentile-bootstrap 95% interval: **[-1.811, -1.668]** points;
- Student-t 95% interval: **[-1.823, -1.647]** points;
- all ten paired seed-level differences are negative; the exact sign test and
  exact paired sign-flip test both give $p=0.001953125$.

The current isotropic rerun also changes the three-dimensional interpretation of
the JONSWAP-plus-transition comparison. OU-III has lower paired mean 3-D
 displacement RMS than OU-II in **all four stationary JONSWAP seas and the
transition**, with all ten paired seeds improving in each of those five
scenarios. This supersedes prose attached to the earlier operating point; the
generated tables/macros are the source of truth for numerical publication
statements.

PM--Stokes remains a separate declared ensemble and is not pooled into the
confirmatory JONSWAP endpoint.

Differences throughout the generated paired-effect tables are `left - right`.
For displacement/attitude errors, a negative value therefore favors the
left-hand estimator. Descriptive intervals should be interpreted with their
paired sample count; they do not create additional independent experiments.

## Provenance contract

`ou_validation_manifest.json` has two provenance layers. `replay_provenance` is
immutable for a set of simulator rows: it records the replay source commit, the
transitive repository-local simulator/filter dependency closure, the simulator
build Makefiles, input hashes, and a SHA-256 of the normalized raw replay CSV.
The build environment is also recorded when a new full replay is generated
(compiler identity/version, Eigen identity, Python, NumPy, and platform), but
those environment fields are informational rather than cross-platform validity
keys.

A later `--restat-from` may recompute summaries, paired effects, tables, or
publication text from those same rows. Before writing anything it compares the
current replay dependency closure with the immutable replay closure and verifies
the raw-row hash. If a filter header, simulator translation unit, shared replay
dependency, or recorded build file changed, restatement fails and requires a
full simulator replay. If only approved analysis/presentation code changed, the
restatement records a separate `restatement` block with its analysis hashes and
source-bundle hash; it never changes the replay commit or replay dependency
hashes.

`tools/ou_evidence_contract.py --check` enforces the replay pins and result
inventory. In a Git checkout it can additionally validate the recorded replay
commit; in a GitHub source ZIP, release archive, Zenodo artifact, or copied tree
without `.git`, it falls back to read-only manifest/file-hash verification and
retains the same dependency-hash checks. `--auto` cannot promote a legacy
restatement or old rows to a new estimator revision.

The existing isotropic bundle was migrated to this schema without replay only
after Git history proved that its replay-producing dependency closure had not
changed since the recorded full replay and that the current normalized raw CSV
was byte-identical to the archived full-replay CSV. No numerical replay metric
or paired statistic was changed by that provenance migration.

## Running

Fetch the versioned simulation data and run the short integration profile:

```bash
make fetch-sim-data
python3 tools/ou_validation.py --mode smoke
```

Run the full declared profile directly:

```bash
python3 tools/ou_validation.py --mode full
```

The production CI path is `.github/workflows/ou-validation.yml`. Full mode is
sharded across three validation and three robustness jobs, combines each study,
checks the manuscript/evidence contract, and then commits the regenerated
bundles and mirrored publication inputs together. The repository `build`
workflow calls that full regeneration on `main`; a manual full dispatch is also
available.

Useful controls include `--stationary-input` (repeatable),
`--skip-nonstationary`, `--duration-sec`, `--window-sec`, `--jobs`, and
`--bootstrap-resamples`. Smoke mode is an integration check only and is not
inferential evidence.

### Restating a bundle without re-simulating

A pure change to summary/statistical presentation can reuse the archived raw
rows:

```bash
python3 tools/ou_validation.py \
  --restat-from reports/results/ou_validation/ou_validation.json \
  --output-dir reports/results/ou_validation
```

Restatement cannot substitute for a replay when any file in the recorded
implementation closure has changed. The evidence contract will report the
mismatch and require a new `--mode full` generation. This prevents a statistical
or formatting restat from laundering rows produced by an older estimator into a
current-code publication bundle.

## Outputs

The output directory contains:

- `ou_validation_raw.csv`: one row per completed filter/mode replay with
  continuous metrics and seed identifiers; no simulator regression-gate result;
- `ou_validation_summary.csv`: sample size, mean, sample standard deviation and
  confidence intervals;
- `ou_validation_paired_effects.csv`: paired mean differences, bootstrap
  intervals, Cohen's $d_z$, and Hedges' $g_z$;
- `ou_validation.json`: protocol, calibration points, raw observations, and all
  statistics;
- `ou_validation_manifest.json`: provenance and self-hashes, including
  implementation and analysis-pipeline hashes;
- `ou_validation_tuning_points.tex`: exact frozen operating points;
- transition CSV/SVG diagnostics and generated LaTeX/SVG publication inputs.

The follow-on OU-III sensitivity/degradation protocol is documented in
[`ou-robustness.md`](ou-robustness.md). It reuses the same paired seed triplets
but does not alter the confirmatory comparison or select a new reported operating
point from reference errors.

## Direct simulator controls

Both OU executables accept repeated `--input PATH` arguments. Validation-related
environment variables include:

- `W3D_IMU_SEED`, `W3D_INIT_SEED`, or combined `W3D_SEED`;
- `W3D_VALIDATION_WINDOW_SEC` for the machine-readable metrics line;
- `W3D_WRITE_TIMESERIES=0` to suppress large diagnostic CSVs;
- `W3D_TUNING_MODE=adaptive`, `fixed_nominal`, or `fixed_oracle`, plus the
  family-specific `W3D_FIXED_*` values;
- `W3D_AW_COV_SYNC=periodic` or `reconfigure` for the $a_w$ covariance policy.

If no seed variable is present, the canonical deterministic realization and its
historical random-number draw order are retained for simulator regression use.
