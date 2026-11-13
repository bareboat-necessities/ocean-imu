# PR 746 validation

The estimator and simulator sources are those validated at
81d16672fa4ff398f9d84247ccb2201f18f29e62. The subsequent archive-fixture repair
adds the missing shared test headers to copied source archives; it does not
change a filter or a replay. Main is synchronized through
b915b0470911a49de6ddcf2ff6dad79cb966f295.

## Completed verification

- The serial build/test inventory was resumed after execution-session
  interruptions. All C++ suites completed, including the 64 paired native
  legacy/default replay gates. The final remaining-suite command exited zero.
  This is not described as a single uninterrupted `make all` success.
- All 1152 validation tests are checked with Eigen 3.4.0 explicitly selected.
  The initial local attempt without `EIGEN_INCLUDE_DIR` failed its device
  compass dependency lookup; the six source-archive/compass tests subsequently
  passed with that environment variable.
- CI on 81d16672 passed quality-gates (strict warnings, Clang-tidy, cppcheck,
  ASan, UBSan, coverage and Python), Arduino, heave-only, wave-direction, both
  device checks, both lever-arm studies, TFG validation and the non-promoting
  OU3 source-bound proof workflow.
- Full OU validation/robustness run 43361798603 genuinely regenerated all six
  shards, both combined bundles and the TFG companion. Publication failed on
  two archive-fixture checks because the copied archive omitted tests/common.
  That omission is repaired; its four archive tests pass without Git metadata,
  including rejection of altered estimator sources and compile flags.

## Remaining publication checks

The broad fingerprint intentionally includes every test. The archive repair
therefore triggers another full replay on its exact revision. Final publication
must use those newly generated bundles and verify both the replay-source and
results fingerprints. The pending CI status is not reported as successful.

## Interpretation

The ordinary acquired-reference six-seed, two-hour comparison gives horizontal
accelerometer-bias RMS 72.80 to 68.52 mm/s² for shared QMEKF and 55.22 to
56.88 mm/s² for shared TFG. Both paired 95% intervals include zero. The physical
integration correction is supported, but reliable live unknown-motion bias
recovery is not established. Firmware conversion delay and hardware timing
remain unqualified; the end-to-end stability theorem remains open.
