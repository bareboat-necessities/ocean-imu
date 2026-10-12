#!/bin/bash -e
set -e



./mag_continuous-test
./device_heave-test
./nlo-sim
./local_gravity-test
./mag_turns-test
