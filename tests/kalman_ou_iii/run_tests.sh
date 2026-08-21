#!/bin/bash -e

./kalman_ou_iii-sim
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
./live_basin_diagnostic
python3 ../../tools/ou_live_basin_interval_proof.py --repo-root ../..
./rs_law-test
./startup_init-test
