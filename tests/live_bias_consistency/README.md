# Live accelerometer-bias consistency experiments

Baseline: `9018fa4b57d9bd10983ae4e77ff95ce509551e17` (main when this task began).
The executable reuses `CarmaBenchTruth`, `W3dSimCommon`, `Bench3::fill_config`,
`SharedModesFusion`, and native QMEKF/TFG correction, Joseph and reset kernels.
The test-only adapter is the earlier follow-up probe, extended for independent
core/frontend timing changes, persistent-state adaptation freezes and recovery.
It does not replace native prediction or measurement algebra.

## Build and checks

Set `EIGEN_DIR` to Eigen 3.4.0 (vendored Eigen is preferred; the Makefiles retain
the system `/usr/include/eigen3` fallback). From the repository root:

```sh
make -C tests/live_bias_consistency EIGEN_DIR="$EIGEN_DIR"
(cd tests/live_bias_consistency && bash run_tests.sh)
EIGEN_INCLUDE_DIR="$EIGEN_DIR" make all EIGEN_DIR="$EIGEN_DIR"
```

The additional tests cover changing multi-axis rates, second-order convergence,
variable intervals and missed samples, exact rotation increments, first-sample
priming, gyro-bias sensitivity, native covariance propagation, correlated-noise
bounds, endpoint lever-arm corrections, rejected IMU packets, native sensor and
oracle Jacobians, and samplewise parity with an uninstrumented fusion wrapper.
Configuration-time selection and the immediate setter path must also produce
identical state and full covariance through startup, handoff and bias activation.
Dense Schmidt/Joseph tests additionally hold the bias mean, preserve its
cross-covariance and check both force and magnetic kernels. Existing quality
gates and thresholds are untouched.

Build the *same* harness against untouched main for a genuine binary baseline:

```sh
git worktree add --detach /tmp/ocean-imu-main 9018fa4b57d9bd10983ae4e77ff95ce509551e17
make -C tests/live_bias_consistency EIGEN_DIR="$EIGEN_DIR" \
  BASELINE_SRC=/tmp/ocean-imu-main/src live-bias-baseline
python tools/live_bias_consistency.py --stage baseline --jobs 2 \
  --binary tests/live_bias_consistency/live-bias-baseline
python tools/live_bias_consistency.py --stage baseline --label compatibility \
  --jobs 2 --binary tests/live_bias_consistency/live-bias-followup
```

`BIAS_UNCHANGED_BASELINE` excludes calls to new APIs; production headers come
from the untouched worktree. No historical research branch or downloaded
research artifact is required. Use a fresh baseline object build if changing
`BASELINE_SRC` (remove `baseline-*.o` and `baseline-*.d`).

## Reproduce results

Python dependencies: NumPy, pandas, SciPy and Matplotlib. No simulated dataset
download is needed for these synthetic experiments. `make all` separately uses
the repository's pinned oceanography-waves-lib v1.3.1 moving-vessel archive.

```sh
python tools/live_bias_consistency.py --stage screen --jobs 4
python tools/live_bias_consistency.py --stage long --jobs 3
python tools/live_bias_consistency.py --stage noise --jobs 3
python tools/live_bias_consistency.py --stage recovery --jobs 2
python tools/live_bias_consistency.py --stage late_recovery --noise schmidt \
  --label late_recovery --jobs 6
python tools/live_bias_consistency.py --stage calibration --jobs 4
python tools/live_bias_consistency.py --stage combined --jobs 3
python tools/live_bias_consistency.py --stage production --jobs 3
python tools/live_bias_consistency.py --stage production_recovery --jobs 4
python tools/main_gyro_rollout.py --jobs 3
python tools/live_bias_consistency.py --stage startup --jobs 2
python tools/live_bias_consistency.py --stage latency --jobs 2
OPENBLAS_NUM_THREADS=1 python tools/live_bias_window_batch.py --jobs 2
python tools/live_bias_csv_timing.py
python tools/live_bias_csv_replay.py --jobs 2
python tools/live_bias_csv_replay.py --all-records --localize --jobs 4
python tools/live_bias_consistency_analysis.py
```

All stages use seeds 1001–1006, 200 Hz IMU and about 25 Hz magnetometer. Each
seed/motion produces one continuous JONSWAP/vessel-RAO history with constant
geographic yaw, or independent 8-degree/0.13 Hz roll and 5-degree/0.19 Hz pitch
with zero translation. Sensor noise and constant turn-on biases are the shared
simulator models; bias random walk is disabled in physical truth. All scores
use every 200 Hz sample, independently of plotting decimation.

`screen`: 1200 s, independent core/frontend changes and sea-state/coupling
freezes at 600 s. `long`: 7200 s wave cases and 3600 s rocking; late scores are
5400–7200 s and 2400–3600 s. `recovery`: 3600 s with additive estimated-bias
`[.04,-.03,.02]` m/s² or 0.01-rad tilt injections at 600 s, retaining covariance;
a separate arm corrects the magnetic reference at that time. `late_recovery`
runs the gyro correction from startup and Schmidt holding with the same native
prior, injects at 3600 s, and continues through 7200 s; unperturbed paired arms
separate recovery from an existing floor. Its corrected-reference arm starts
with the ordinary acquired reference and replaces it at 3600 s. `startup` enables
the input convention from the first sample. `latency` generates magnetic
measurements 0, 5 or 20 ms old but applies them at arrival; it diagnoses an
uncompensated delay, without asserting that these delays describe hardware.

`noise` changes only declared white-noise parameters, using the known synthetic
sensor noise: gyro density `GYR_SIGMA*sqrt(.005)`, magnetic SD `MAG_SIGMA_UT`,
or all three (accelerometer SD `ACC_SIGMA`). It is an ablation, not a retuned
shipping default. Both geometries receive the same parameters in every case. `calibration`
separates a matched initial bias prior, a covariance-only late change, all
matched white noise plus prior, Schmidt holding, and Schmidt holding with the
matched prior. Both exact-reference and ordinary learned-reference cases run
3600 s. `combined` runs the timing correction from the first sample for 7200 s
with legacy holding, Schmidt holding, or Schmidt plus the matched prior. Its
magnetic-reference arms distinguish conditional diagnostics from ordinary
estimation.

`endpoint`, `core`, `frontend`, and `both` localize gyro preprocessing. The
original raw sensor checksum is unchanged across these arms; the production
cores/frontend perform preprocessing. The optional `interval` diagnostic
instead generates `Log(R_prev^T R_current)/dt` as a simulated sensor measurement.
Do not describe interval and instantaneous inputs as identical effective rates.

`learned_fixed` is the ordinary arm: acquired/refined magnetic reference, fixed
correct zero hard iron, and no truth supplied to estimation. `exact_fixed`
provides an explicitly diagnostic exact magnetic reference. `oracle` additionally
observes true world acceleration with SD .001 m/s² through native update kernels
before each accelerometer correction, beginning at 600 s. It is artificial and
is not proposed as a deployed sensor. `correct_reference` is another labelled
truth intervention. Truth attitude and IMU biases are otherwise used only for
sensor generation, scoring and labelled diagnostic injections.

Every run retains state/covariance and physical history across interventions.
Applied magnetic corrections are instrumented at the native kernel; enabled
bias gates are read from the native core. The wrapper's magnetic update count
is not interpreted as directional information. Reference/offset changes,
innovations, physical covariance, coverage and nuisance conditioning are logged.

The offline window experiment spans 60 s / about 17 wave cycles and uses the
labelled exact-reference diagnostic history; its field is not an acquired
calibration. It retains all 18 uncertain CARMA motion states, per-node attitudes and a constant gyro
bias, uses actual magnetic packet times at about 2.5 Hz, and integrates the
original 200 Hz gyro between nodes. It performs sparse Gauss–Newton
relinearization with line search. It uses the captured CARMA dynamics/excitation
and a stationary motion prior, **not** an acceleration oracle or a second copy
of the filter's acceleration estimate. Its independent zero-mean bias prior
has SD `8e-3*g/sqrt(3)`, the synthetic turn-on uniform distribution's SD;
sensitivity arms multiply motion covariance by four or bias SD by three.
These are declared priors. Their information is separated from sensor factors
when eliminating nuisance variables. Target scales are .01 rad and .05 m/s². The five-dimensional tilt/bias
Schur analysis expresses initial attitude perturbations about world X/Y and
also eliminates world yaw; the full six-dimensional attitude/bias analysis is
retained separately.
Numerical rank regularizers and PSD-floor diagnostics are retained in JSON.

Window estimates never feed back into the live filter, so they cannot count
any of its samples twice. They are an offline reference diagnostic, with a
coarser motion grid and constant-window biases, not an equivalent shipping
smoother. Seed-level Student intervals use six independent seeds. Serial
coverage fractions and NEES are diagnostics, not independent binomial samples
or a claim of calibrated uncertainty. No stability theorem or shipping
counterexample is asserted.

The runner always executes fresh runs. Binary hashes, commands and paired
sensor checks are recorded per stage. Full histories are regenerable under the
ignored `runs/` directory. Committed seed-level tables and bounded SVG/CSV
vertices are the review evidence; figures are never used to calculate RMS.

The existing pinned-record replay has an explicit corrected input contract:

```sh
W3D_CSM_GYRO_INPUT=instantaneous W3D_WRITE_TIMESERIES=0 \
  tests/kalman_carma_shared_modes/kalman_carma_shared_modes-sim \
  --input /path/to/pinned/wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv
```

The same environment variable works for the dedicated TFG shared-mode replay.
This is labelled diagnostic (`GATED=0`); `--regression` requires the unchanged
baseline contract. Rounded CSV timestamps are audited by
`live_bias_csv_timing.py`; fixed index cadence is authoritative for the pinned
records. An irregular/dropped real recording needs actual acquisition timing.


## Live defaults and broader rollout

OU-II, OU-III, TFG, CARMA/shared-mode and NLO live paths default to
instantaneous sampling. Their Config exposes `gyro_sample_mode`; select
`IntervalRate` explicitly for interval averages/preintegrated inputs. Native
Kalman `time_update` mathematical kernels keep their interval-rate compatibility
contract; normal live wrappers set the corrected mode automatically. The marine frontend retains its endpoint integration by default, with a separate
`frontend_gyro_sample_mode` diagnostic selector. Held-bias policies retain their
existing defaults; CARMA Schmidt handling is an explicit audited option.

The extended focused tests run independent multi-axis rotation/covariance
checks through native OU-II/OU-III/TFG. Whole-wrapper interval references cover
constant-rate default startup and the optional centered frontend under changing
rotation; samplewise Config/setter parity checks the selected default path.
`production_recovery` repeats the one-hour injections/reference repair with
these actual defaults, keeping covariance and physical histories continuous.
NLO gets the same default reference-path check; PII retains endpoint feedback by default after a heave replay regression, and its optional centered path gets a reference check. NLO's handoff test requires the first
live attitude to equal the current bootstrap seed, without a second prediction.

`production` uses the unmodified sampling/held-bias Config defaults with common benchmark tuning (core centering, endpoint
frontend, legacy holding), for six seeds of wave/rocking and labelled oracle controls. `main_gyro_rollout`
runs OU-II, OU-III and TFG on identical continuous constant-yaw wave histories;
each run steps both an explicit interval-rate control and the actual default
wrapper on the same sensor samples. It uses the shared `CarmaBenchTruth` and
sensor-noise infrastructure, keeps the acquired magnetic reference, fixes the
correct zero hard iron, and never supplies acceleration or attitude truth to
estimation. Runs last 7200 s and score 5400–7200 s. Family priors are retained;
this is a within-family timing comparison, not retuning between motion models.

Bare shared-mode replay now uses and gates production defaults. `--regression`
selects the unchanged historical baseline explicitly. The CSV comparison
runs both plus the isolated timing-only arm; no baseline threshold is weakened.

Independent CARMA, OU-III, TFG and TFG-CARMA also expose `--regression` for
the locked historical integrator. Their test scripts run that fixture with
every original limit and the corrected-default fixture. The report lists the
distinct yaw references and retains common-limit comparisons. Reproduce with
`W3D_COLLECT_ALL_GATES=1 W3D_WRITE_TIMESERIES=0`; run each binary first with
`--regression`, then without it. Legacy replay rejects `--no-noise`, `--nomag`
and non-default seeds.

After rendering figures, retain the ignored full-rate `runs/` directory outside
`reports/results` before checking the published evidence fingerprint: the gate
intentionally fingerprints every protected results file, including ignored
files. Restore the directory temporarily for history re-analysis.
