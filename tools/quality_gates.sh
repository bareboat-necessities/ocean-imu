#!/usr/bin/env bash
#
# Quality gates for ocean-imu. Each subcommand is one CI job, and runs the same
# way locally as it does in .github/workflows/quality-gates.yml — so a gate that
# fails on a PR can be reproduced with a single command.
#
#   tools/quality_gates.sh warnings     compile everything with -Wall -Wextra -Wpedantic -Werror
#   tools/quality_gates.sh sanitizers   run the unit-test suites under ASan + UBSan
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
LIEPP_INC="-isystem $REPO_ROOT/third_party/Lie-plusplus/include"

# Test directories whose suites are fast enough to run under sanitizers on a
# per-PR gate. The OU/TFG simulation suites are deliberately absent: under
# -fsanitize=address,undefined they run past 15 minutes each, and the plain
# build already exercises them. See docs/quality-gates.md.
SANITIZER_DIRS=(spike_filter spectrum detrend imu_calibrate wave_dir wave_sim)

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

gate_warnings() {
  step "warnings: -Wall -Wextra -Wpedantic -Werror"
  local log; log="$(mktemp)"
  local n=0
  # -Werror is applied per TU so one failure does not mask the rest.
  while read -r f; do
    [[ -z "$f" ]] && continue
    n=$((n + 1))
    printf '%s\n' "$f"
  done < <(tu_list) | xargs -P "$JOBS" -I TU bash -c '
      out=$("$0" -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
              -I src '"$EIGEN_INC"' '"$LIEPP_INC"' -c "$1" -o /dev/null 2>&1)
      [[ -n "$out" ]] && { echo "--- $1"; echo "$out"; }
      exit 0
    ' "$CXX" TU > "$log" 2>&1

  if [[ -s "$log" ]]; then
    cat "$log"
    fail "warnings gate: $(grep -c 'error:\|warning:' "$log") diagnostic(s)"
  else
    ok "warnings gate: $(tu_list | wc -l) translation units clean"
  fi
  rm -f "$log"
}

gate_sanitizers() {
  step "sanitizers: AddressSanitizer + UndefinedBehaviorSanitizer"
  # -O1 keeps the instrumented build inside the runner's memory budget while
  # staying fast enough to be worth running.
  local sanflags="-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined"
  export UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=0"
  export ASAN_OPTIONS="detect_leaks=1:abort_on_error=0"

  for d in "${SANITIZER_DIRS[@]}"; do
    [[ -f "tests/$d/Makefile" ]] || { fail "sanitizers: tests/$d has no Makefile"; continue; }
    local blog rlog
    blog="$(mktemp)"; rlog="$(mktemp)"
    if ! make -j "$SAN_JOBS" -C "tests/$d" build CPPFLAGS="$sanflags" > "$blog" 2>&1; then
      tail -30 "$blog"; fail "sanitizers: $d failed to build"
      rm -f "$blog" "$rlog"; continue
    fi
    ( cd "tests/$d" && bash ./run_tests.sh ) > "$rlog" 2>&1
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

gate_clang_tidy() {
  step "clang-tidy"
  command -v clang-tidy >/dev/null || { fail "clang-tidy not installed"; return; }
  local log; log="$(mktemp)"
  # clang-tidy re-parses each TU with the clang frontend; the OU/TFG simulation
  # drivers take many minutes each because of their Eigen instantiations, so the
  # gate covers the unit-test TUs and lets the sims be covered by the compiler
  # warnings gate instead.
  git ls-files 'tests/*.cpp' \
    | grep -Ev '(-sim|certificate-sim|neighborhood-sim)\.cpp$' \
    | xargs -P "$JOBS" -I TU bash -c '
        clang-tidy --quiet "$1" -- -std=c++20 -I src '"$EIGEN_INC"' '"$LIEPP_INC"' 2>/dev/null
      ' _ TU > "$log" 2>&1

  if grep -q 'warning:\|error:' "$log"; then
    grep -A3 'warning:\|error:' "$log" | head -60
    fail "clang-tidy: $(grep -c 'warning:\|error:' "$log") diagnostic(s)"
  else
    ok "clang-tidy clean"
  fi
  rm -f "$log"
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
  for d in "${SANITIZER_DIRS[@]}"; do
    [[ -f "tests/$d/Makefile" ]] || continue
    make -j "$JOBS" -C "tests/$d" build CPPFLAGS="$covflags" LDFLAGS="--coverage" >/dev/null 2>&1 \
      || { fail "coverage: $d failed to build"; continue; }
    ( cd "tests/$d" && bash ./run_tests.sh ) >/dev/null 2>&1 || true
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
  warnings)   gate_warnings ;;
  sanitizers) gate_sanitizers ;;
  clang-tidy) gate_clang_tidy ;;
  cppcheck)   gate_cppcheck ;;
  python)     gate_python ;;
  coverage)   gate_coverage ;;
  all)
    gate_warnings
    gate_python
    gate_cppcheck
    gate_clang_tidy
    gate_sanitizers
    gate_coverage
    ;;
  *)
    sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 2
    ;;
esac

exit "$STATUS"
