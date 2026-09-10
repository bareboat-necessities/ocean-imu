#!/bin/bash -e

sim_status=0
W3D_COLLECT_ALL_GATES=1 ./kalman_ou_iii-sim || sim_status=$?
./accel_vibration_guard-test
./imu_lever_arm-test
./kalman_ou_common-test
./aw_covariance_policy-test
./acc_bias_ou-test
./channel_freeze-test
./wave_period-test
./mag_hard_iron-test
./continuous_mag_hard_iron-test
./tuner_coupling-test
./tuner_schedule-test
./wave_band_sigma-test
./iss_contract-test
./rs_law-test
./live_entry_audit-test
./startup_init-test

exit "$sim_status"
