#!/bin/bash -e
set -e

./rapid_motion-test


./mag_continuous-test
./device_heave-test
./nlo-sim
./local_gravity-test
./mag_turns-test
