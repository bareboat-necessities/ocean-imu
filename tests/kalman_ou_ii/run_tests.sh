#!/bin/bash -e
set -e
make -f Makefile -f ../common/StationaryDeviceRegression.mk stationary_device-test
./stationary_device-test

sim_status=0
W3D_COLLECT_ALL_GATES=1 ./kalman_ou_ii-sim || sim_status=$?
./tuner_schedule-test
./iss_contract-test
./startup_init-test
./regularizer_floor-test
./pseudo_law-test

exit "$sim_status"
