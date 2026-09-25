#!/bin/bash -e

./imu_calibrate-test
./accel_cal-test
./accel_cal_gravity_override-test
./accel_cal-replay calibrate_accel_replay_sample.log
