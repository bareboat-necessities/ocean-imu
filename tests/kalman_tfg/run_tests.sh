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
./tfg_process_covariance-test
./tfg_orchestrator-test
./tfg_rs_axis_factors-test
./tfg_handoff_heave-test
./tfg_device_heave-test\n./tfg_device_tempcomp-test
sim_status=0
W3D_COLLECT_ALL_GATES=1 ./kalman_tfg-sim || sim_status=$?

exit "$sim_status"
