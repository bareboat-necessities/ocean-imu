#!/bin/bash -e
set -e

./qmekf_geometry-test

./mag_continuous-test
./mag_timing-test
./ahrs-qmekf-sim
