#!/bin/bash -e

./calibration_safety-test
python3 test_sketch_temperature.py
./imu_calibrate-test
./accel_cal-test
./accel_cal_gravity_override-test
./accel_cal-replay calibrate_accel_replay_sample.log
