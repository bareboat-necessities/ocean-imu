# Wave startup and live handover

This host-only experiment reads the current v1.4.0 following-sea and head-sea
records. `StartupLearner` has no access to truth, filenames, wave height, vessel
speed, inclination, declination or field magnitude. It consumes calibrated
body-axis accelerometer, instantaneous gyroscope and timestamp-matched
magnetometer samples. This full-history solver serves as the experimental reference. Production uses the fixed-memory implementation in `src/kalman_common/WaveStartup.h`.

`ShippingBaseline` runs the shared `ProxyStartupFusion` with the previous startup policy through the
OU-III wrapper and stops at its **first live handover, before a live magnetic
correction**. The test-only private-access shim reads the acquired magnetic
reference; it does not modify any estimator equation or API. Bias and
quaternion seeds are read from the actual live state.


The OU-II, OU-III and CARMA wrappers and the TFG-OU proxy bootstrap use the
same streaming solver (880 bytes on the measured x86-64 host). It accumulates gyro-transported field and force moments,
fits a constant body gyro bias, then constructs magnetic north and down at the
handover sample. The fitted bias and field enter the live core along with the
attitude. The gyro interval history, configured covariance, accelerometer-bias
release policy and tuner clocks continue through handover. The ordinary proxy
handover remains available when the fit lacks a fresh field or fails its
residual checks, and for explicit startup hard-iron estimation. Set
`wave_startup_enabled=false` for a paired control. The host driver's
`--production-pair` switch scores both policies on the same sensor stream
and records `joint_seed_applied` for each live handover.

```sh
make -C tests/startup_learning
tests/startup_learning/startup-test
python3 tools/startup_learning.py --seeds 1,2,3 --offsets 0,7
```

The runner verifies both archive/CSV pins using `tools/sim_dataset.py`, makes
one immutable noisy power-on history per trial, and feeds all methods that
same prefix. Only sensor generation and scoring access the CSV truth. The
baseline is scored at its own handover epoch; candidate `matched` rows use
that exact epoch. Separate `fixed` rows use 10, 15, 20, 25, 30 and 35 seconds.
No estimator receives samples beyond the horizon being scored.

## What is being estimated

| Quantity | Startup treatment |
|---|---|
| Gravity magnitude | Known, 9.80665 m/s² for these simulations; never fitted |
| Gravity direction | Unit vector, two physical orientation degrees of freedom |
| Attitude | Current-epoch body-to-magnetic-NED quaternion, all three axes |
| Local magnetic field | Unconstrained 3-D vector in the initial gyro frame; transformed to magnetic NED |
| Gyroscope bias | Regularized finite-window fit, retaining a prior in weak directions |
| Accelerometer bias | Calibration-prior seed; **not claimed learned** in a short moving window |
| Angular rate | Latest measured body rate minus gyro-bias seed, in rad/s |
| Next prediction history | Last raw gyro sample and its acquisition timestamp |
| Position, velocity, wave state | Outside this experiment; cannot be initialized from these attitude statistics alone |
| True north | Not learned; external declination is needed to convert magnetic north to geographic north |

Existing hard/soft-iron and axis calibration remain prerequisites. A coherent
magnetic offset can look like a different local field on a short, weakly
excited record. This experiment explicitly measures that degradation.

## Methods

All candidates propagate a relative body-to-initial-frame quaternion from the
first sample, with no fictitious prediction before that sample. Instantaneous
gyro endpoints use trapezoidal integration. The terminal attitude is therefore
at the **last acquisition epoch**, not the mean heading in the window.

* `uniform`: average gyro-transported force and magnetic observations using
  the gyro calibration prior.
* `tapered`: replace the force mean by the article's normalized
  `30 x² (1-x)²` quadrature. Its endpoint value and derivative vanish, which
  suppresses bounded-displacement contamination. Only the direction of this
  intermediate force statistic is used; gravity magnitude is fixed.
* `magfit`: estimate gyro bias by minimizing the scatter of the transported
  magnetic samples, eliminating the unknown constant field analytically.
  Every bias trial repropagates the complete gyro history.
* `joint`: additionally penalize the tapered linear slope of transported
  force. This informs drift about the magnetic axis. Its 35-second model
  standard deviation is 0.005 m/s³, enlarged as `(35/T)^3` for shorter windows.
  It is a finite-window wave-contamination assumption, not thousands of
  independent gravity measurements or a proof that wave acceleration averages
  exactly to zero.
* `joint_cal`: the same joint fit also estimates a constant residual hard-iron
  vector, with a 1.2 µT calibration prior. The field is still eliminated
  analytically. This is a separate ablation, because fitting an unnecessary
  calibration parameter costs accuracy when the magnetometer is already clean.

The gravity direction and horizontal magnetic projection define a right-handed
magnetic-NED frame. That frame is composed with the **terminal** relative
rotation. Magnetic magnitude and inclination are outputs, not inputs.

## Measurement contract and limits

Sensor axes are body NED; quaternion maps body NED to magnetic world NED.
Specific force is negative down at rest. Magnetic values are microtesla, rates
rad/s and acceleration m/s². Magnetometer samples run at 25 Hz, IMU at 200 Hz.
An arriving magnetic sample must first be associated with its acquisition
epoch; arbitrary delivery-time reuse is not supported. Invalid samples,
duplicate/backward timestamps, gaps above 0.1 s and windows beyond 35 s are
rejected without mutating the history.

Noise profiles separate ideal sensors, noisy/biased IMU with calibrated
magnetometer, the repository's residual magnetic calibration model, and a
10 µT body-frame field step at 15 s. All three sensor streams have separate
seed streams; residual calibration and measurement noise are also separate.

The development set uses seeds 1–3 at offsets 0 and 7 s. The recorded validation
set uses seeds 11–15 at 67, 307 and 607 s, after freezing both fits. The focused
CI workflow repeats both directions and checks paired yaw improvement, tilt,
gyro-bias error and complete fixed-horizon outputs. Separate 35-second gates
include every startup, including shipping deadline misses: worst calibrated
yaw ≤2°, worst tilt ≤1° in both profiles, residual-calibration yaw RMS ≤2.5°,
and gyro-bias vector RMS ≤0.04°/s. Run it locally with:

```sh
python3 tools/startup_learning.py --seeds 11,12,13,14,15 \
  --offsets 67,307,607 --profiles calibrated,residual --check \
  --output reports/results/startup_learning/validation
python3 tools/report_startup_learning.py reports/results/startup_learning/validation
```

The optional report renderer needs Matplotlib; the experiment and checks use
only a C++20 compiler, Eigen and the Python standard library. `metrics.csv.gz`
contains six checkpoints plus exact matched-time comparisons, including the
actual quaternion, field, bias/rate seeds and conditional gyro-bias covariance.

Readiness is a diagnostic heuristic based on elapsed time, horizontal field,
magnetic scatter, force norm, weighting sensitivity, fresh magnetic data and
four-sigma consistency with the supplied calibration priors. The report counts
false readiness against an explicit target of yaw error ≤2° and tilt ≤1°.
A small residual is not an accuracy certificate. In particular, an unavailable
heading, a deadline fallback, and a verified accurate seed must remain distinct
when comparing estimator policies.

The field-step control revealed that a constant-field optimizer can absorb
interference into a very large gyro/hard-iron correction. Such fits remain in
the evidence but cannot be marked ready. The fit equations and accuracy
settings remain fixed; the prior-consistency guard is an additional model
rejection, not a truth-based choice between candidate attitudes.

The experiment reports the fitted gyro covariance conditional on the calibration
and wave model. Attitude, reference, calibration and wave uncertainties are
correlated. The deployed filters keep their configured initial covariance and
continue the gyro clock from the last observed sample.

## Reproducing the force-slope comparison

Run both weighting models with the same seeds, offsets and sensor profiles.
The comparison checks input hashes, epoch pairing and shipping control values.
CI runs both cohorts and uploads the full raw metrics. Example for the newer
480-start cohort:

```sh
for model in fixed scaled; do
  flag=""
  if [ "$model" = scaled ]; then flag=--scale-force-slope; fi
  python3 tools/startup_learning.py --seeds 31,32,33,34,35 --offsets 41,191,421 \
    --profiles calibrated,residual --check $flag --output "/tmp/startup-$model"
done
python3 tools/compare_startup_slope.py /tmp/startup-fixed /tmp/startup-scaled /tmp/confirmation.csv
```

## Reproducing the handover component comparison

Use the pinned v1.4.0 input directories produced by the startup runner. Both
commands build the real OU-III replay. The study changes only which startup
values enter the core. The control checks compare the old and new endpoints
against separate builds using the intact production handover. Each sea
direction contributes its 4 m JONSWAP record to these checks.

```sh
python3 tools/startup_handover_components.py --verify-controls --jobs 2 --seeds 7 41 \
  --data tests/startup_learning/data --output /tmp/startup-components-early
python3 tools/startup_handover_components.py --long --verify-controls --jobs 2 --seeds 7 41 101 \
  --data tests/startup_learning/data --output /tmp/startup-components-steady
```

## Tracing the live startup gain

The retention driver reuses the component driver's compiler setup and the
normal OU-III simulator. Source-anchored temporary headers select experimental
handover or refinement interventions. Every policy consumes the same full
1,200 s input. The receipt includes input, executable and source hashes.
The optional `--seed` argument selects an independent sensor realization.

The archived bias-prior experiment used estimator commit `c3481ce`. Its analysis
driver was published in `abfd578`. To reproduce the archived scores, prepare a
separate worktree with those two recorded inputs, then run the policy set
from the bias-prior receipt.
Run the following from the repository root and set `replay` to a new directory.

```sh
replay=/tmp/ocean-imu-bias-prior-replay
git worktree add --detach "$replay" c3481cef635faee2d77efb5e4beba85a903e37e1
git -C "$replay" restore --source=abfd57856fc008556bf64e004d82f4b99f5e90f7 --worktree -- \
  tools/startup_live_retention.py tools/startup_handover_components.py tools/startup_source_snapshot.py
cd "$replay"
python3 tools/sim_dataset.py --case all --dest tests/startup_learning/data
python3 tools/startup_live_retention.py --study bias_prior --data tests/startup_learning/data \
  --output /tmp/startup-bias-prior-scores \
  --policies bias_prior bias_prior_early bias_prior_tilt bias_prior_tilt_early deployed
```

The receipt records the SHA-256 values for the estimator files and the analysis
driver. The component runner writes matched production logs under
`production-controls/` and instrumented logs alongside `records.csv`.

The shared source snapshot is captured before compilation and checked after
compilation and before writing the results; changed source or HEAD rejects the
receipt. Each completed replay log is preserved for audit.

```sh
for study in timeline reference refinement retention bias_prior frontend; do
  python3 tools/startup_live_retention.py --study "$study" --jobs 2 \
    --data tests/startup_learning/data --output "/tmp/startup-live-$study"
done
```

Complete the four additional long-run component interventions with:

```sh
python3 tools/startup_handover_components.py --long --jobs 2 --seeds 7 41 101 \
  --arms previous_tilt previous_yaw previous_field previous_gyro_bias \
  --data tests/startup_learning/data --output /tmp/startup-components-restored
```

The deterministic magnetic-step probe prints accepted and rejected cases.
The noisy probe reproduces the 1,440 matched clean/step controls. The archived
previous results retain the acceptance defect before the continuity guard:

```sh
c++ -O1 -std=c++20 -I src -isystem third_party/eigen \
  tests/startup_learning/field-step-probe.cpp -o /tmp/startup-field-step
/tmp/startup-field-step
```

```sh
c++ -O1 -std=c++20 -I src -isystem third_party/eigen \
  tests/startup_learning/field-noise-probe.cpp -o /tmp/startup-field-noise
/tmp/startup-field-noise
```

The guarded production startup qualification uses the same two cohorts as
`production-handover.json`. Its receipt records the base commit and exact
modified header, binary, inputs and outputs; it is separate from the full
live-navigation matrix. Both controls and wave-start measurements are in
`reports/results/startup_learning/field_steps`.

## Tracing gravity through the live handover

Build test-only versions of the OU-III and TFG wave replays and run each of the
16 pinned records through the previous and deployed startup, one-time exact
tilt and field controls, and a sensor-only wave-acceleration seed. The
accelerometer and gyro measurements seen by each paired run are identical.
Run from the repository root with Eigen installed and the replay compiler
selected.

```sh
python3 tools/sim_dataset.py --case all --dest tests/startup_learning/data
python3 tools/startup_gravity_oracle.py --data tests/startup_learning/data \
  --families ou3 tfg --seeds default 7 41 --jobs 4 \
  --policies previous deployed previous_exact_tilt exact_tilt exact_field exact_both measured_wave \
  --output gravity-handover
```

The command above exercises the checked-out production configuration.
To reproduce the 2° TFG baseline from the gravity-retention study on its
original guarded source, create a detached worktree and copy the two
test-only files into it before running the same commands there.

```sh
git worktree add --detach ../ocean-imu-gravity-baseline abfd57856fc008556bf64e004d82f4b99f5e90f7
cp tools/startup_gravity_oracle.py ../ocean-imu-gravity-baseline/tools/
cp tests/startup_learning/GravityHandoverOracle.h ../ocean-imu-gravity-baseline/tests/startup_learning/
```

Each compressed CSV contains either the handover, the first live updates, or
the scored time windows. The receipt records source, executable, input and
output hashes. The [study](../../docs/startup-gravity-retention.md) gives the
meaning of each intervention and the observed accuracy.

## Replaying the TFG tilt-uncertainty comparison

The previous guarded TFG build is `abfd578` and the accepted-fit uncertainty
change is `68f0978`. Build the same executable in a separate checkout at
each revision using the pinned compiler. For each of the 16 records, use
the identical input path and sensor seed with both binaries. This example
runs the 8.5 m following JONSWAP case on the default sensor draw.

```sh
make -C tests/kalman_tfg kalman_tfg-sim
record=tests/startup_learning/data/following/wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv
env -u W3D_SEED OCEAN_IMU_SIM_CASE=following W3D_WRITE_TIMESERIES=0 \
  W3D_COLLECT_ALL_GATES=1 W3D_VALIDATION_WINDOW_SEC=900 \
  W3D_VALIDATION_SEGMENTS=early:40:60,late:300:1200 \
  tests/kalman_tfg/kalman_tfg-sim --input "$record"
```

For the other sensor draws set `W3D_SEED` to 7, 41 or 101. Set
`OCEAN_IMU_SIM_CASE=head` for the head-sea records. The per-run
`VALIDATION_SEGMENT` lines contain the paired scores. The archived receipt
retains the exit status against the limits present at the tested revision.
