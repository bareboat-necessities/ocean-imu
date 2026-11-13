# Tangent-bias shared-mode research

Baseline: `7712a0e046bb653839a4d582ec893d7671ac1e5f` (current main at study
start, including the measured-endpoint integration and sample-path optimization).

`KalmanTangentCarmaSharedModes` implements the 27-state **reduced tangent-bias
hybrid**: the velocity/bias subgroup is `T(SE(3))`, with the same Euclidean
CARMA motion model. It does not implement the full 30-state `T(SE_2(3))`
INS with a virtual position-input bias. The article distinguishes the two
and derives why the standard autonomous INS-error theorem cannot simply
be substituted for the CARMA dynamics.

The benchmark uses one TFG adapter instantiation for TFG and tangent modes;
its geometry switch is exactly the setting enabled by the new estimator
constructor. This avoids duplicating the truth, sensors, startup, adaptation,
correction, covariance write or benchmark code.

## Build and native tests

```sh
make -C tests/kalman_tangent_shared_modes build
make -C tests/kalman_tangent_shared_modes test
make all
```

Use serial compilation on small hosts. Eigen is vendored when available,
otherwise `/usr/include/eigen3` is used. The root CMake build includes all
three executables. `geometry-test` checks group/adjoint algebra, finite
reset derivatives, structured versus dense covariance reset, PSD, the
fixed-heading rocking information control and the Joseph penalty for a
nonoptimal gain. The smoke envelope is fixed before holdout evaluation;
it rejects gross numerical/accuracy failures and is not an adoption claim.

## Controlled cohort

```sh
python3 tools/tangent_shared_modes_study.py --profile development \
  --raw /tmp/tangent-development --output /tmp/tangent-development-summary
python3 tools/tangent_shared_modes_study.py --profile acceptance \
  --raw /tmp/tangent-holdout --workers 2
python3 tools/tangent_shared_modes_report.py --raw /tmp/tangent-holdout
```

Development seeds are 7301–7302; untouched holdout seeds are 4101–4106.
The acceptance profile uses 24 physical conditions, six independent seeds,
and two additional wave noise/prior controls. Normal rows run 900 s;
`long` runs 7200 s. Every variant in a row consumes the same measured
history. Sensor hashes are required to agree, and binary/source hashes and
commands are retained. Resume is allowed only with the same binary hash
and command. A nonzero experiment exit is retained as a failure, never
silently dropped from the study.

Normal filters receive only IMU, temperature and magnetometer samples.
The truth is the existing independent directional spectral/RAO synthesis.
The ordinary learned-reference arm does not supply the actual field to the
estimator. `exact_reference` and `reference_error` explicitly set a fixed
field as diagnostic controls. No acceleration oracle is used.

Primary scoring is 20 Hz from full-rate live updates. Covariances/traces are
sampled at 1 Hz. Scores are not recomputed from publication traces. The
reported position and drift are wave-relative quantities; no sensor in
this configuration observes absolute geographic position or steady drift.
Metrics and units are defined in the standalone article.

The optional bias veto changes the bias gain and uses the complete Joseph
update. It is a conditional fading-memory excitation monitor, **not** a
full marine observability certificate. Defaults remain disabled. `-G`
enables it before startup; `-LG` enables it at 300 s through a dedicated
configuration method that leaves coupling/adaptation/covariance intact.
Native held-bias semantics remain in force. TAN uses the existing TFG
retraction while a bias is held, avoiding a second-order shift of that mean.

Early perturbation arms apply a fixed estimated-state error at first Live;
late arms do so at 300 s. Thus startup is shared and live convergence is
separable. A convergence time is the first 30 s continuously below 1 degree
attitude or 0.01 m/s² three-dimensional bias error; missing convergence is
reported as missing, not averaged as a negative time.

## Untouched-main comparison

```sh
git worktree add --detach /tmp/ocean-imu-tangent-baseline \
  7712a0e046bb653839a4d582ec893d7671ac1e5f
make -C tests/kalman_tangent_shared_modes tangent-baseline \
  BASELINE_SRC=/tmp/ocean-imu-tangent-baseline/src
tests/kalman_tangent_shared_modes/tangent-baseline --scenario wave \
  --duration 900 --seed 7301 --filters QSM,TFG,CARMA --output /tmp/old
tests/kalman_tangent_shared_modes/tangent-study --scenario wave \
  --duration 900 --seed 7301 --filters QSM,TFG,CARMA --output /tmp/new
```

Require equal trace bytes and sensor hashes. Compare scores excluding runtime
and host object size, which are expected to differ. All baseline factory
units compile with `BASELINE_SRC` before the working source include path;
`TANGENT_BASELINE` removes access to new optional APIs. No baseline header is
copied into the working source tree.

## Pinned vessel-response replay

```sh
make ensure-sim-data
cd tests/kalman_tangent_shared_modes
W3D_WRITE_TIMESERIES=0 W3D_COLLECT_ALL_GATES=1 ./tangent-replay
```

This uses the exact existing shared replay and the **unchanged TFG bars**
from `tests/common/TfgSharedRegressionBars.h`. A failed candidate replay
remains an experimental negative result; its bars are not widened and it
does not replace the accepted default replay. Original QMEKF, TFG,
independent CARMA and unrelated root quality gates remain mandatory.

## Window observability

```sh
python3 tools/tangent_shared_modes_observability.py \
  /tmp/tangent-holdout/wave_4101_native --field-sd 0 \
  --output /tmp/information-known.json
python3 tools/tangent_shared_modes_observability.py \
  /tmp/tangent-holdout/wave_4101_native --field-sd 2 \
  --output /tmp/information-uncertain.json
```

`--field-sd inf` leaves the reference free. `--prior-scale 4` changes only
the declared motion prior in this offline diagnostic. All poses, gyro bias,
motion states and reference nuisance parameters are eliminated by sparse
least-squares projection. The script never reads truth columns. It records
solver convergence and numerical covariance floors. One-second gyro
quadrature and frozen exported CARMA parameters are explicit diagnostic
approximations. Sensor-only freeing of acceleration is a relaxation, not
a theorem about bounded marine motion.

## Publication and feasibility

The standalone article is
`doc/kalman_tangent_shared_modes/kalman_tangent_shared_modes.tex`.
`tools/tangent_shared_modes_report.py` generates tables, confidence
comparisons and bounded SVG/PGF charts using the repository's shared
`plots/vector_export.py`. Raw traces stay outside git. A focused CI workflow
builds native tests and the article from committed, source-bound evidence.
The full cohort is separately reproducible; routine CI runs a fixed smoke.

Runtime is instrumented **process CPU time** per filter update, not elapsed
wall time from parallel runs. Object storage includes the heap-owned core,
excludes allocator metadata, and uses the host ABI. There is no claim of
measured ESP32-S3 runtime, memory high-water mark or target build acceptance.
