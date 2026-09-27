#!/usr/bin/env bash
# Copyright 2026, Mikhail Grushinskiy
# Optional reproduction against actual released source; requires the v2.3.2 tag.
set -euo pipefail
repo=$(git rev-parse --show-toplevel)
legacy=$(mktemp -d)
trap 'rm -rf "$legacy"' EXIT
git -C "$repo" archive 2ba31e4e255a2b03f014cfbb9c9352d58d49f791 src tests/imu_calibrate tests/common sensors | tar -x -C "$legacy"
eigen="$repo/third_party/eigen"
if [[ ! -f "$eigen/Eigen/Dense" ]]; then eigen=/usr/include/eigen3; fi
"${CXX:-g++}" -O3 -std=c++20 -funroll-loops -fno-finite-math-only -march=native \
  -DCAL_BASELINE_V232 -I"$legacy/src" -I"$eigen" \
  "$repo/tests/imu_calibrate/calibration_accuracy-test.cpp" -o "$legacy/release-accuracy"
"$legacy/release-accuracy"
make -C "$legacy/tests/imu_calibrate" EIGEN_DIR="$eigen" accel_cal-test
(cd "$legacy/tests/imu_calibrate" && ./accel_cal-test)
