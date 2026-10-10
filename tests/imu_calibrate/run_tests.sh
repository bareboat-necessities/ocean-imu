#!/bin/bash -e
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

python3 test_wizard_config.py
./calibration_accuracy-test
./mag_field_consistency-test
./mag_hand_motion-test
./calibration_workflow-test
./calibration_safety-test
./mag_diagnostics-test
./mag_rotation-test
python3 test_sketch_temperature.py
./imu_calibrate-test
./accel_cal-test
python3 test_release_accuracy.py
./accel_cal_gravity_override-test
./accel_cal-replay calibrate_accel_replay_sample.log
