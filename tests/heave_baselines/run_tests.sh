#!/bin/bash -e
set -e
./heave_baselines-test
./heave_baselines-sim
./heave_baselines-sim --frontend mahony_slow
./heave_baselines-sim --frontend truth
