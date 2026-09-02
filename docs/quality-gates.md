# Quality gates

Sanitizer, static-analysis and coverage layer for the C++ and Python in this
repository. Every gate runs the same way locally as it does in CI:

```
tools/quality_gates.sh warnings     # -Wall -Wextra -Wpedantic -Werror
tools/quality_gates.sh sanitizers   # ASan + UBSan over the unit-test suites
tools/quality_gates.sh clang-tidy   # bugprone / cert / performance / portability
tools/quality_gates.sh cppcheck     # warning / performance / portability
tools/quality_gates.sh python       # ruff + py_compile
tools/quality_gates.sh coverage     # gcovr line/branch report
tools/quality_gates.sh all
```

CI runs them in `.github/workflows/quality-gates.yml`. A red gate on a PR names
the subcommand that reproduces it.

Eigen is found via `$EIGEN_DIR`, then `third_party/eigen`, then
`/usr/include/eigen3` — the same order the test Makefiles use.

## What each gate found when it was added

These are the numbers the gates were tuned against, so the next person does not
have to re-derive them.

| gate | baseline | now |
|---|---|---|
| warnings | 1 real finding (`-Wreorder` in `KalmanQMEKF.h`) | 52 TUs clean |
| sanitizers | never run | 6 suites, 0 ASan / 0 UBSan / 0 leak findings |
| cppcheck | 5 findings — 3 false positives, 2 in frozen files | clean |
| clang-tidy | 18 findings, 16 of them noise | 2 blockers, see below |
| python | 96 ruff findings across 122 files | clean |
| coverage | no tooling | 65.6% lines from the six sanitizer dirs |

The codebase was in much better shape than the absence of tooling suggested.
The warnings gate cost exactly one real fix to turn on.

## The one thing that is not green

`clang-tidy` is **non-blocking** (`continue-on-error`) for one specific reason:

```
src/kalman_ou_ii/SeaStateFusionFilter_OU_II.h:1375: error: constexpr variable
  'C_HS' must be initialized by a constant expression
```

`std::sqrt` in a `constexpr` initializer is a GCC extension. Clang rejects it,
and clang-tidy uses the clang frontend, so any translation unit that includes
`SeaStateFusionFilter_OU_II.h` or `_OU_III.h` fails to parse. That is also why
this library does not currently build with clang at all.

The fix is two lines per file — `std::sqrt(2.0f)` → `std::numbers::sqrt2_v<float>`,
which is bitwise identical (`0x3FB504F3`) and so changes no numerical result:

```diff
-        constexpr float K = std::sqrt(2.0f) / std::numbers::pi_v<float>;
+        constexpr float K = std::numbers::sqrt2_v<float> / std::numbers::pi_v<float>;
```

It is not applied here because both files are in the frozen replay dependency
closure (below). Once it lands with an evidence regeneration, drop the
`continue-on-error` from the `clang-tidy` matrix entry.

## The frozen replay closure

`tests/validation` pins the committed OU evidence to an immutable list of 34
implementation files. Any byte change to one of them fails the `validate` job:

```
validation: replay dependency differs from replay provenance: <file>
```

Only a genuine full simulator regeneration may create new replay provenance, so
**these gates never rewrite those files, and neither should an autofix.** The
list lives in `reports/results/*/[study]_manifest.json` under
`replay_provenance.implementation_files`; it covers most of `src/tuner/`,
`src/freq/`, `src/kalman_ou_*/`, `src/wave_dir/`, `src/util/W3dSimCommon.*`, and
the two OU test Makefiles.

Practical consequences:

- `tools/quality_gates.sh` has no `--fix` mode. Adding one would need a frozen
  path filter first.
- Three suppressions in `tools/cppcheck-suppressions.txt` and the clang-tidy
  blocker above exist because the finding is inside the closure.
- Changing compiler flags in `tests/kalman_ou_ii/Makefile` or
  `tests/kalman_ou_iii/Makefile` invalidates the evidence too. That is why the
  warnings gate compiles its own translation units instead of adding flags to
  the existing Makefiles.

## Scope choices, and why

**Sanitizers cover six unit-test suites**, not all fourteen:

```
spike_filter spectrum detrend imu_calibrate wave_dir wave_sim
```

The OU/TFG/AHRS/NLO simulation suites are excluded. Under
`-fsanitize=address,undefined` they run past the 15-minute mark each — a first
attempt timed out three of them at 900 s — which is too slow for a per-PR gate.
They are still compiled and run by the existing `build` matrix in `build.yml`.
If they ever need sanitizer coverage, run them on a schedule rather than per PR.

**Sanitizer builds are memory-hungry.** Instrumenting the Eigen template
instantiations OOM-killed `cc1plus` at `-j4` on a 16 GB machine. `SAN_JOBS`
defaults to 2 for that reason; raise it only if you know the runner is bigger.

**The warnings gate skips `src/AtomS3R/`.** Those files include `<Wire.h>` and
the M5Unified headers, which exist only in the Arduino core; they are covered by
the `build-MCU` matrix in `build.yml`. Watch for this — before the exclusion was
explicit, `ImuCalWizardRunner.cpp` was failing to compile and contributing zero
warnings to the count, which looked like success.

**No formatter.** `black --check` wants to reformat 120 of 122 Python files, and
there is no `.clang-format`. A formatting diff that large buries real findings,
so the Python gate checks for defects only (`F`, `E9`, `B`) and leaves layout
alone.

**`I` and `l` are allowed as names** (`E741`/`E743` are off). `I` is the
identity matrix. Renaming it would make the source harder to check against the
papers it implements.

## Coverage

`tools/quality_gates.sh coverage` writes `reports/coverage/index.html` and
`summary.txt`, and CI uploads them as an artifact. It reports 65.6% lines /
38.2% branches over `src/`, measured from the six sanitizer suites only — it is
a floorless report today, not a gate.

To turn it into a gate, set `COVERAGE_MIN_LINES`:

```
COVERAGE_MIN_LINES=65 tools/quality_gates.sh coverage
```

Do not set a floor from the current number without first widening the suite
list: the figure reflects six of fourteen directories, so it understates the
estimators and would ratchet on the wrong baseline.
