#!/bin/bash -e
set -e
./heave_only-test

# Qualify the eight scored inputs before running any method. Explicit inputs
# also avoid loading truth repeatedly for unscored wave families.
mapfile -t records < <(PYTHONPATH=../../tools python3 - <<'PY'
from sim_dataset import REFERENCE_CSV_SHA256, input_provenance
input_provenance(REFERENCE_CSV_SHA256)
print(*REFERENCE_CSV_SHA256, sep='\n')
PY
)
if [ "${#records[@]}" -ne 8 ]; then
    echo 'heave regression requires all eight verified wave records' >&2
    exit 1
fi
inputs=()
for record in "${records[@]}"; do inputs+=(--input "$record"); done
./heave_only-sim "${inputs[@]}"
./heave_only-sim --frontend proxy "${inputs[@]}"
./heave_only-sim --frontend truth "${inputs[@]}"
