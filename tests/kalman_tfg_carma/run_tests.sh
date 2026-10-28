#!/bin/bash -e

# The shebang's -e is ignored when this is invoked as `bash run_tests.sh`,
# which is how the root Makefile runs it. Set it explicitly so a failing test
# actually fails the build.
set -e

# Deterministic replay of the eight scored dataset records with TFG-CARMA3
# (CARMA(4,1) on x, y and z on the two-frame Lie group, full reset
# transport).  doc/kalman_tfg_carma's charts and results table are built from
# these replays; a record above its regression bar fails the run after all
# eight are scored.  The filter's mathematical checks and the paired Monte
# Carlo comparison live in tests/kalman_carma (tfg_carma3-test, carma3_bench).
sim_status=0
W3D_COLLECT_ALL_GATES=1 ./tfg_carma3-sim || sim_status=$?
exit "$sim_status"
