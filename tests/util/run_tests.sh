#!/bin/bash -e
set -e

./ins-helpers-test
W3D_COLLECT_ALL_GATES=1 ./replay-gates-test valid >/dev/null
W3D_COLLECT_ALL_GATES=1 ./replay-gates-test quiet >/dev/null
for failure in nan inf truth_nan missing acc_bias overflow gyro_bias roll; do
    if W3D_COLLECT_ALL_GATES=1 ./replay-gates-test "$failure" >/dev/null 2>&1; then
        echo "replay quality gate accepted $failure" >&2
        exit 1
    fi
done
echo 'replay quality gates reject invalid histories and independent bias/attitude regressions'
