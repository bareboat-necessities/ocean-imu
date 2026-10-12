#!/bin/bash -e

# Callers use `bash run_tests.sh`, which ignores the shebang's -e.
set -e

./mag_continuous-test
./pii_observer-adaptive
./local_gravity-test
make -f Makefile -f ../common/InsSketchHeaveDrift.mk ins_sketch_heave-test
./ins_sketch_heave-test
make -f Makefile -f ../common/InsSketchHeaveDrift.mk ins_sketch_replay-test
./ins_sketch_replay-test
