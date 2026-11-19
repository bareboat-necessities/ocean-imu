#!/usr/bin/env bash
#
# Quality gates for ocean-imu. Each subcommand is one CI job, and runs the same
# way locally as it does in .github/workflows/quality-gates.yml — so a gate that
# fails on a PR can be reproduced with a single command.
#
#   tools/quality_gates.sh warnings     compile everything with -Wall -Wextra -Wpedantic -Werror
#   tools/quality_gates.sh sanitizers   run filter replays and unit suites under ASan + UBSan
#   tools/quality_gates.sh clang-tidy   static analysis (.clang-tidy)
#   tools/quality_gates.sh cppcheck     static analysis (tools/cppcheck-suppressions.txt)
#   tools/quality_gates.sh python       ruff over tools/, tests/, plots/ (ruff.toml)
#   tools/quality_gates.sh coverage     gcovr line/branch report over src/
#   tools/quality_gates.sh all          every gate above, in order
#
# Eigen is located the same way the test Makefiles locate it: $EIGEN_DIR, then
# third_party/eigen, then /usr/include/eigen3.
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

: "${CXX:=g++}"
: "${JOBS:=$(nproc 2>/dev/null || echo 2)}"

# Sanitizer builds of the Eigen-heavy sims are memory-hungry: instrumenting
# those template instantiations OOM-killed cc1plus at -j4 on a 16 GB runner.
# Keep the sanitizer build serial-ish regardless of how wide $JOBS is.
: "${SAN_JOBS:=2}"

if [[ -n "${EIGEN_DIR:-}" && -f "$EIGEN_DIR/Eigen/Dense" ]]; then
  EIGEN_INC="-isystem $EIGEN_DIR"
elif [[ -f "$REPO_ROOT/third_party/eigen/Eigen/Dense" ]]; then
  EIGEN_INC="-isystem $REPO_ROOT/third_party/eigen"
else
  EIGEN_INC="-isystem /usr/include/eigen3"
fi

# Spectrum has no host executable. Coverage uses the same utility suites.
UNIT_DIRS=(spike_filter detrend imu_calibrate wave_dir wave_sim)
FUSION_CASES=(ou2 ou3 tfg carma tfg_carma shared tfg_shared tangent)

STATUS=0
step()  { printf '\n\033[1m== %s\033[0m\n' "$*"; }
fail()  { printf '\033[31mFAIL\033[0m %s\n' "$*"; STATUS=1; }
ok()    { printf '\033[32mok\033[0m   %s\n' "$*"; }

# Translation units compiled by the warnings gate: every test driver plus the
# library .cpp files that build off-target. Compiling a TU pulls in the headers
# it uses, so this covers src/ without needing a separate per-header pass.
#
# src/AtomS3R/ is excluded because it includes <Wire.h> and the M5Unified
# headers, which exist only in the Arduino core. Those files are covered by the
# build-MCU matrix in build.yml, which compiles them for the real target.
tu_list() {
  git ls-files 'tests/*.cpp' 'src/*.cpp' | grep -v '^src/AtomS3R/'
}

# CI partitions this ordered inventory. With no shard settings the local
# command runs every translation unit, using the same flags as CI.
gate_host_audit() {
  local gate="$1" shard="${SHARD_NUMBER:-1}" count="${SHARD_COUNT:-1}"
  if (( shard < 1 || shard > count )); then
    fail "$gate: invalid shard $shard/$count"; return
  fi
  local -a units command selected=()
  local -a includes=(-I src -isystem "${EIGEN_INC#-isystem }" -isystem "$REPO_ROOT/third_party/Lie-plusplus/include")
  if [[ "$gate" == warnings ]]; then
    mapfile -t units < <(tu_list)
    command=("$CXX" -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror "${includes[@]}" -c TU -o /dev/null)
  else
    mapfile -t units < <(git ls-files 'tests/*.cpp' | grep -Ev '(-sim|certificate-sim|neighborhood-sim)\.cpp$')
    command=(clang-tidy --quiet TU -- -std=c++20 "${includes[@]}")
  fi
  local i
  for ((i=shard-1; i<${#units[@]}; i+=count)); do selected+=("${units[i]}"); done
  if (( ${#selected[@]} == 0 )); then
    fail "$gate: empty shard $shard/$count"; return
  fi
  step "$gate: shard $shard/$count, ${#selected[@]} of ${#units[@]} translation units"
  local log rc=0
  log=$(mktemp)
  printf '%s\0' "${selected[@]}" | xargs -0 -P "$JOBS" -I TU "${command[@]}" > "$log" 2>&1 || rc=$?
  if (( rc != 0 )) || grep -qE 'warning:|error:' "$log"; then
    cat "$log"
    fail "$gate: shard $shard/$count failed (xargs exit $rc)"
  else
    ok "$gate clean"
  fi
  rm -f "$log"
}

gate_unit_sanitizers() {
  step "sanitizers: AddressSanitizer + UndefinedBehaviorSanitizer"
  # -O1 keeps the instrumented build inside the runner's memory budget while
  # staying fast enough to be worth running.
  local sanflags="-O1 -g1 -fno-omit-frame-pointer -fsanitize=address,undefined -fno-sanitize-recover=all"
  export UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=1"
  export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=1}:halt_on_error=1"

  for d in "${UNIT_DIRS[@]}"; do
    [[ -f "tests/$d/Makefile" ]] || { fail "sanitizers: tests/$d has no Makefile"; continue; }
    local blog rlog
    blog="$(mktemp)"; rlog="$(mktemp)"
    if ! make -B -j "$SAN_JOBS" -C "tests/$d" build CPPFLAGS="$sanflags" > "$blog" 2>&1; then
      tail -30 "$blog"; fail "sanitizers: $d failed to build"
      rm -f "$blog" "$rlog"; continue
    fi
    ( cd "tests/$d" && bash -e ./run_tests.sh ) > "$rlog" 2>&1
    local rc=$?
    local ub as lk
    ub=$(grep -c 'runtime error:' "$rlog")
    as=$(grep -c 'ERROR: AddressSanitizer' "$rlog")
    lk=$(grep -c 'ERROR: LeakSanitizer' "$rlog")
    if (( ub || as || lk )); then
      grep -B2 -A20 'runtime error:\|ERROR: AddressSanitizer\|ERROR: LeakSanitizer' "$rlog" | head -80
      fail "sanitizers: $d — ubsan=$ub asan=$as leak=$lk"
    elif (( rc != 0 )); then
      tail -30 "$rlog"; fail "sanitizers: $d test suite exited $rc"
    else
      ok "sanitizers: $d clean"
    fi
    rm -f "$blog" "$rlog"
  done
}

# Build into a fresh temporary directory so an ordinary or coverage build
# can never satisfy a filter sanitizer target.
filter_sanitizer_case() (
  set -eu
  local suite="$1" tmp source
  tmp=$(mktemp -d)
  trap 'rm -rf "$tmp"' EXIT
  local flags=(-std=c++20 -O1 -g1 -fno-omit-frame-pointer -fsanitize=address,undefined
               -fno-sanitize-recover=all -DEIGEN_UNROLLING_LIMIT=0 -I src
               -isystem "${EIGEN_INC#-isystem }")
  export ASAN_OPTIONS="${ASAN_OPTIONS:-detect_leaks=1}:halt_on_error=1"
  export UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1
  if [[ "$suite" == observers ]]; then
    local tests=(tests/ahrs/qmekf_geometry-test.cpp tests/ahrs/mahony-frame-test.cpp
                 tests/freq/pll-coarse-test.cpp tests/heave_only/heave_only-test.cpp
                 tests/nlo/device_heave-test.cpp tests/pii_observer/local_gravity-test.cpp)
    for source in "${tests[@]}"; do
      step "$source"
      "${SAN_CXX:-clang++}" "${flags[@]}" "$source" -o "$tmp/test" || exit 1
      "$tmp/test" || exit 1
    done
  else
    step "sanitizers: $suite fusion replay"
    "${SAN_CXX:-clang++}" "${flags[@]}" "-DFILTER_${suite^^}" tools/quality/fusion_smoke.cpp -o "$tmp/test" || exit 1
    "$tmp/test"
  fi
)

gate_sanitizers() {
  local requested="${1:-all}" suite
  if [[ "$requested" == all || "$requested" == unit ]]; then
    gate_unit_sanitizers
    [[ "$requested" == unit ]] && return
  fi
  for suite in "${FUSION_CASES[@]}" observers; do
    if [[ "$requested" == all || "$requested" == "$suite" ]]; then
      filter_sanitizer_case "$suite" || fail "sanitizers: $suite"
      [[ "$requested" == "$suite" ]] && return
    fi
  done
  [[ "$requested" == all ]] || fail "unknown sanitizer suite: $requested"
}

gate_cppcheck() {
  step "cppcheck"
  command -v cppcheck >/dev/null || { fail "cppcheck not installed"; return; }
  local log; log="$(mktemp)"
  cppcheck --enable=warning,performance,portability \
           --std=c++20 --language=c++ --quiet --inline-suppr \
           --suppressions-list=tools/cppcheck-suppressions.txt \
           --error-exitcode=0 \
           --template='{severity}: {id}: {file}:{line}: {message}' \
           -I src src/ > "$log" 2>&1

  if [[ -s "$log" ]]; then
    cat "$log"
    fail "cppcheck: $(wc -l < "$log") finding(s)"
  else
    ok "cppcheck clean"
  fi
  rm -f "$log"
}

gate_python() {
  step "python: ruff"
  command -v ruff >/dev/null || { fail "ruff not installed (pip install ruff)"; return; }
  if ruff check .; then ok "ruff clean"; else fail "ruff findings above"; fi

  # py_compile catches syntax errors in files ruff's rule set does not reach.
  if python3 -m py_compile tools/*.py tests/validation/test_*.py 2>/dev/null; then
    ok "py_compile clean"
  else
    python3 -m py_compile tools/*.py tests/validation/test_*.py
    fail "py_compile"
  fi
}

gate_coverage() {
  step "coverage: gcovr over src/"
  command -v gcovr >/dev/null || { fail "gcovr not installed (pip install gcovr)"; return; }
  find tests -name '*.gcda' -o -name '*.gcno' | xargs -r rm -f

  local covflags="-O0 -g --coverage -fprofile-abs-path"
  for d in "${UNIT_DIRS[@]}"; do
    [[ -f "tests/$d/Makefile" ]] || continue
    make -B -j "$JOBS" -C "tests/$d" build CPPFLAGS="$covflags" LDFLAGS="--coverage" >/dev/null 2>&1 \
      || { fail "coverage: $d failed to build"; continue; }
    ( cd "tests/$d" && bash -e ./run_tests.sh ) >/dev/null 2>&1 || true
  done

  mkdir -p reports/coverage
  gcovr --root . --filter 'src/' --exclude 'third_party/' \
        --html-details reports/coverage/index.html \
        --txt reports/coverage/summary.txt \
        --print-summary || { fail "gcovr"; return; }
  ok "coverage report: reports/coverage/index.html"

  # Informational for now. Set COVERAGE_MIN_LINES to turn it into a floor once
  # the sanitizer/coverage dir list covers the estimators.
  if [[ -n "${COVERAGE_MIN_LINES:-}" ]]; then
    local pct
    pct=$(gcovr --root . --filter 'src/' --exclude 'third_party/' --print-summary 2>/dev/null \
          | sed -n 's/^lines: \([0-9.]*\)%.*/\1/p' | head -1)
    if awk "BEGIN{exit !($pct < $COVERAGE_MIN_LINES)}"; then
      fail "coverage: lines ${pct}% below floor ${COVERAGE_MIN_LINES}%"
    else
      ok "coverage: lines ${pct}% >= floor ${COVERAGE_MIN_LINES}%"
    fi
  fi
}

case "${1:-all}" in
  warnings)   gate_host_audit warnings ;;
  sanitizers) gate_sanitizers "${2:-all}" ;;
  clang-tidy) gate_host_audit clang-tidy ;;
  cppcheck)   gate_cppcheck ;;
  python)     gate_python ;;
  coverage)   gate_coverage ;;
  all)
    gate_host_audit warnings
    gate_python
    gate_cppcheck
    gate_host_audit clang-tidy
    gate_sanitizers
    gate_coverage
    ;;
  *)
    sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 2
    ;;
esac

exit "$STATUS"
