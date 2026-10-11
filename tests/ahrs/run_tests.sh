#!/bin/bash -e
set -e

./rapid_motion-mahony

./rapid_motion-qmekf

./mag_continuous-test
./mag_timing-test
./ahrs-qmekf-sim
