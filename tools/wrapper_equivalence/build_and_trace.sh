#!/usr/bin/env bash
# Build the wrapper characterization drivers against a source tree and write
# their traces.  Usage: build_and_trace.sh <src-root> <out-dir>
#
# Compare two revisions with compare_traces.py; see README.md.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$(cd "$1" && pwd)"
OUT="$2"
mkdir -p "$OUT"
EIGEN_DIR="${EIGEN_DIR:-$SRC/../third_party/eigen}"
if [ -f "$EIGEN_DIR/Eigen/Dense" ]; then EIGEN="-isystem $EIGEN_DIR"; else EIGEN="-isystem /usr/include/eigen3"; fi
# -ffp-contract=off by default: with contraction enabled GCC fuses a*b+c into
# FMA wherever inlining happens to expose it, so an edit that only moves code
# between functions can change low-order bits without changing a single
# source-level operation.  Off, the trace is a function of operation order
# alone, which is what a behavior-preserving refactor must keep.  Pass the
# deployed flags through CXXFLAGS to measure contraction drift separately.
CXXFLAGS="${CXXFLAGS:--O3 -std=c++20 -march=native -funroll-loops -fno-finite-math-only -ffp-contract=off}"
build() { g++ $CXXFLAGS $EIGEN -I"$HERE" -I"$SRC" "$@"; }
build -DOU_FAMILY=2 "$HERE/ou_wrapper_trace.cpp" -o "$OUT/ou2_trace" &
build -DOU_FAMILY=3 "$HERE/ou_wrapper_trace.cpp" -o "$OUT/ou3_trace" &
build "$HERE/tfg_wrapper_trace.cpp" -o "$OUT/tfg_trace" &
wait
for t in ou2 ou3 tfg; do "$OUT/${t}_trace" "$OUT"; done
