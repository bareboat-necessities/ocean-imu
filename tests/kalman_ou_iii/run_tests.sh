#!/bin/bash -e
set -e
make -f Makefile -f ../common/StationaryDeviceRegression.mk stationary_device-test
./stationary_device-test
make -f Makefile -f ../common/LocalGravityRegression.mk local_gravity-test
./local_gravity-test
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
./shipping_contract-test
./shipping_transition-test
./sampled_capture-test
./rs_law-test
./startup_init-test
exit "$sim_status"
