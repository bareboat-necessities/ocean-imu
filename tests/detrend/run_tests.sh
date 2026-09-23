#!/bin/bash -e
set -e

./detrend-basic-test

./detrend-wave-test

./detrend-wave3d-test

./detrend-log-frequency-test
