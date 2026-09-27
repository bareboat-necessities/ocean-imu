#!/bin/bash -e

# Callers use `bash run_tests.sh`, which ignores the shebang's -e.
set -e

./device_heave-test
./nlo-sim
./local_gravity-test
./mag_turns-test
