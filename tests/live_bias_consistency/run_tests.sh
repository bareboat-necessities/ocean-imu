#!/bin/bash
set -e
./timing-test
./main-filter-timing-test
./defaults-ou2-test
./defaults-ou3-test
./defaults-tfg-test
./observer-defaults-test
./reinitialization-test
./finite-mean-test
./packet-rejection-test
./schmidt-test
python3 ./noise-bound.py
./audit-sensors
./live-bias-followup --filter qmekf --self-test
./live-bias-followup --filter tfg --self-test
task_bias_dir=$(mktemp -d)
trap 'rm -rf "$task_bias_dir"' EXIT
for geometry in qmekf tfg; do
    ./live-bias-followup --filter "$geometry" --arm learned_fixed --parity --duration 180 --output "$task_bias_dir/$geometry"
    ./live-bias-followup --filter "$geometry" --arm learned_fixed --timing both_start --change 0 --noise schmidt --config-parity --duration 180 --output "$task_bias_dir/config_$geometry"
    ./live-bias-followup --filter "$geometry" --arm learned_fixed --timing defaults --change 0 --config-parity --duration 180 --output "$task_bias_dir/defaults_$geometry"
done
for family in OU2 OU3 TFG; do
    ./main-gyro-rollout "$family" 1001 180 "$task_bias_dir/$family"
done
