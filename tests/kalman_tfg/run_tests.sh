#!/bin/bash -e

# The shebang's -e is ignored when this is invoked as `bash run_tests.sh`,
# which is how the root Makefile runs it. Set it explicitly so a failing test
# actually fails the build.
set -e

./lie_group-test
./convention-test
./ou_chain_identity-test
./tfg_propagation-test
./tfg_jacobians-test
./covariance_transport-test
./tfg_orchestrator-test
./kalman_tfg-sim
