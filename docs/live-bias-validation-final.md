# Live-bias validation receipt

Validated source revision:
`60bbf96bd4a9d3e24f3df37896d5759648a113ac`.
Main is synchronized through
`b915b0470911a49de6ddcf2ff6dad79cb966f295`.

## Local validation

The complete, uninterrupted serial `make all` finished with exit code 0 on
13 November 2026. All 21 suite scripts ran, including the native legacy/default
replays, focused timing and covariance checks, NLO, PII, wave direction and
simulation-convention tests. The validation suite ran 1152 tests in 173.028 s:
OK, with one existing skip.

Eigen 3.4.0 was selected for both C++ Makefiles and the Python-driven device
compass checks. The command was:

```sh
EIGEN_INCLUDE_DIR="$EIGEN_DIR" make all EIGEN_DIR="$EIGEN_DIR"
```

Here `EIGEN_DIR` names the installed Eigen include directory. The completed
local log has SHA-256
`a084178099e360088a59b4c38d40c5ab9b62e51055094041b478330ccc804231`.

## CI validation

On the validated revision, the evidence-contract workflow independently passed
all 1152 tests (one skip). The OU validation gate, TFG validation, heave-only,
wave-direction, both device workflows, both lever-arm studies, source-bound
OU3 proof checks and stability-study rendering passed. All quality gates passed:
Python, cppcheck, strict warnings, coverage, ASan, UBSan and Clang-tidy.

The fresh full-study workflow is
[43370428471](https://github.com/bareboat-necessities/ocean-imu/actions/runs/43370428471).
All six replay shards, both combined bundles and the TFG companion passed.
The final publisher passed manuscript/evidence admission, generated both
fingerprints and committed the bundle as
`6c27f1a074d93e15422a7f53dac2e3989a247968`. The local replay/results fingerprint
check, OU evidence contract and CI-qualified TFG companion check all pass.

The published replay fingerprint is
`ad3b07fcf6de6637ced1d615de38ce1d7f444ab1b0f5875e819c75aa4f536797`;
the results fingerprint is
`07e18bb45a2a21bef51c121f16aa72ee63c83907c8331c49dcf162606e45be9b`.
The pinned simulation ZIP SHA-256 is
`87e3cee9505439ec3c0a3ea8c7516789a8f48e3412314c1dcd483894908d433c`.
The full study was genuinely replayed on the repaired source revision; no old
raw rows were rebound to changed code.

Both Linux and Windows Arduino jobs passed, linking every shipped sketch with
the Arduino IDE default flags on the validated revision.

The estimator, firmware and replay C++ sources are byte-identical to
`81d16672fa4ff398f9d84247ccb2201f18f29e62`, on which all quality gates and Arduino
builds passed. The archive-fixture repair copies the shared replay headers into
source archives; changed estimator sources and compiler flags still fail closed.
The complete new replay obeys the conservative fingerprint rule for that repair.

## Scope

The ordinary acquired-reference six-seed, two-hour comparison gives horizontal
accelerometer-bias RMS 72.80 to 68.52 mm/s² for shared QMEKF and 55.22 to
56.88 mm/s² for shared TFG. Both paired 95% intervals include zero. Reliable
live unknown-motion bias recovery is not established by these results.

Firmware conversion delay and hardware timing remain unqualified. Passing the
finite simulations and source-bound diagnostics does not close the end-to-end
OU3 stability theorem. The original study's environment and source bindings
remain in `reports/results/live_bias_consistency/validation.json`; this receipt
records the completed serial validation and subsequent publication checks.
