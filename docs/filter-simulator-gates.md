# Filter simulator quality gates

The regression baseline is oceanography-waves-lib **v1.4.0**, 28 ft RAO at
5 kn, with ±3° heading oscillations every 20 s. Following and head seas have
independent limits. `docs/moving-vessel-replays.md` describes the inputs and commands.

`tools/filter_sim_gates.py` covers OU-II, OU-III, TFG, legacy TFG, independent
CARMA, QMEKF/TFG/tangent shared-mode CARMA, TFG-CARMA, AHRS, NLO, PII and all
heave-only methods. Current and interval-rate profiles are both scored where
supported. The heave-only methods additionally retain Mahony, proxy and truth
frontends where applicable. There are **544 independently gated cells**:
272 filter/profile/sea combinations in each direction.

All complete records contain 240,000 samples at 200 Hz. Scores use the trailing
900 s, except the AHRS attitude regression, which retains its trailing 60 s.
Each three-axis INS replay gates heave, wave Y displacement, 3-D displacement, roll, pitch, heading,
accelerometer Z/3-D bias and gyroscope 3-D bias. AHRS gates its three attitude
components, PII gates heave/heading, and vertical-only methods gate heave.
Metrics that a method does not estimate are not invented.

`reports/results/filter_sim_gates/per_record_limits.json` records each case,
family, profile, CSV, measured value and limit. The limit is the measured RMS
plus 0.5%, rounded upward to 1e-5. This is a deterministic regression allowance,
not a statistical confidence bound or an accuracy requirement for arbitrary
sensor noise. Input hashes are distinct even though the two archives share
filenames. `RecordReplayLimits.inc` is generated from this reviewed evidence.
`CurrentReplayLimits.h` holds the overall ceilings; the per-record limits
prevent one direction or sea from borrowing another's larger allowance.
The original qualifying receipts are retained as `baseline_following.json`
and `baseline_head.json` beside the gate table and checked against its hashes.
Current replay reports may be refreshed without changing those locked caps.
The Y gate is in metres RMS, independent of the 3-D percentage gate. It is
enabled for all 272 horizontal default/interval cells, including legacy TFG;
vertical-only and explicitly ungated diagnostic configurations do not acquire
an invented horizontal requirement. The steering-separated baseline and all
cap increases are documented in
[`heading_limit_update.md`](../reports/results/filter_sim_gates/heading_limit_update.md).

The fixed-draw host acceptance build is Linux x86-64, Debian GCC
**14.2.0-19** and Eigen 3.4. The compiler package build is part of the protocol:
Ubuntu GCC 14.2.0-4ubuntu2 produces a different startup transition on the
following-sea 4 m JONSWAP record, despite reporting the same upstream version.
The shared CI action extracts five official, SHA256-pinned Debian compiler and
header packages into a private directory, checks the relocated frontend and
compiles/runs a probe; it does not replace the runner's system packages.
CI selects that compiler explicitly and uses the shared AVX2 host flags instead
of selecting AVX-512 or other reductions from the runner CPU. The tangent
replay retains its existing scalar/SSE build. Compiler/SIMD changes can alter
adaptive mode decisions over a long replay; their results need separate
qualification rather than being treated as identical fixed-draw arithmetic.
These host build settings do not change firmware or estimator code.

The v1.4.0 baseline is initialized explicitly from complete measurements.
The old v1.3.1 heading period was 120 s, so its thresholds are not carried over
as if the input bytes had remained unchanged. `tools/lock_simulation_baseline.py
--write-reviewed-baseline --reason "reason for the baseline change"` is a review action, never a CI step. Ordinary runs
never modify acceptance limits. `--measure` retains failures for investigating
a new dataset and is not used by regression CI.

`make sim-regression` and `simulation-matrix.yml` run both directions with all
profiles. Missing/duplicate/short records, nonfinite metrics, invalid case
names and failed numeric gates are errors. Complete scores and failure logs
are retained before returning a failing status. The case-aware study matrix
also rejects partial dataset invocations rather than silently appending them.
Startup/engine studies additionally enable the regular Gerstner, cnoidal and
Fenton controls with `W3D_ALL_WAVE_TYPES`. Those controls retain their diagnostic
global bars; enabling them does not bypass the eight irregular-wave record caps.

Current JSON/CSV metrics and SVG charts are in
`reports/results/sim_v1.4.0/{following,head}`. Their receipts record the input
and executable hashes and actual exit status. Older speed-control bundles are
historical comparisons, explicitly bound to their original v1.3.1 input pins.
