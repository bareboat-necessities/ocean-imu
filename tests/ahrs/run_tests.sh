#!/bin/bash -e
set -e

./qmekf_geometry-test
./mahony-frame-test
./ahrs-qmekf-sim
