#!/bin/bash -e
set -e

./ins-helpers-test
./marine-turn-test
W3D_COLLECT_ALL_GATES=1 ./replay-gates-test valid >/dev/null
W3D_COLLECT_ALL_GATES=1 ./replay-gates-test quiet >/dev/null
W3D_COLLECT_ALL_GATES=1 ./replay-gates-test y_valid >/dev/null
W3D_COLLECT_ALL_GATES=1 ./replay-gates-test y_disabled >/dev/null
for failure in nan inf truth_nan missing acc_bias overflow gyro_bias roll y; do
    if W3D_COLLECT_ALL_GATES=1 ./replay-gates-test "$failure" >/dev/null 2>&1; then
        echo "replay quality gate accepted $failure" >&2
        exit 1
    fi
done
echo 'replay quality gates reject invalid histories and independent Y/bias/attitude regressions'
