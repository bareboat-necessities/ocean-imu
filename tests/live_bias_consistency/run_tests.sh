#!/bin/bash
set -e
./timing-test
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
done
