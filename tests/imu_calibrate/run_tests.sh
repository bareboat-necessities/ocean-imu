#!/bin/bash -e
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

./calibration_accuracy-test
./calibration_workflow-test
./calibration_safety-test
python3 test_sketch_temperature.py
./imu_calibrate-test
./accel_cal-test
python3 test_release_accuracy.py
./accel_cal_gravity_override-test
./accel_cal-replay calibrate_accel_replay_sample.log
