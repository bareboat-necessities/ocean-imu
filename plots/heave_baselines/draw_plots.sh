#!/bin/bash -e

# Needs the time series of the OU simulators and of heave_baselines-sim
# (all three front ends), run with W3D_WRITE_TIMESERIES=1 in their test
# directories.
for frontend in mahony proxy truth; do
  python3 ./heave_baselines-plots.py --tests-dir ../../tests --output-dir . --frontend "$frontend"
done
