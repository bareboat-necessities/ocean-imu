#!/bin/bash -e
set -e

./rapid_motion-test
./ins_heave-test


./mag_continuous-test
./device_heave-test
./nlo-sim
./local_gravity-test
./mag_turns-test
