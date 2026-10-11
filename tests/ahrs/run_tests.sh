#!/bin/bash -e
set -e

./rapid_motion-mahony

./rapid_motion-qmekf

./qmekf_geometry-test

./mag_continuous-test
./mag_timing-test
./ahrs-qmekf-sim
