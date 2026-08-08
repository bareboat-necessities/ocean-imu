#!/bin/bash -e

./kalman_ou_iii-sim
./kalman_ou_common-test
./aw_covariance_policy-test
./acc_bias_ou-test
./channel_freeze-test
./wave_period-test
./mag_hard_iron-test
./tuner_coupling-test
