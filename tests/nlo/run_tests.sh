#!/bin/bash -e

# Callers use `bash run_tests.sh`, which ignores the shebang's -e.
set -e

./mag_continuous-test
./device_heave-test
make -f Makefile -f ../common/InsSketchHeaveDrift.mk ins_sketch_heave-test
./ins_sketch_heave-test
make -f Makefile -f ../common/InsSketchHeaveDrift.mk ins_sketch_replay-test
./ins_sketch_replay-test
./nlo-sim
./local_gravity-test
./mag_turns-test
