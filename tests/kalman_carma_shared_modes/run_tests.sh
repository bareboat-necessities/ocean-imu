#!/bin/bash -e

# The shebang's -e is ignored when this is invoked as `bash run_tests.sh`,
# which is how the root Makefile runs it.
set -e

# EXPERIMENTAL standalone shared-modes QMEKF (src/kalman_carma_shared_modes,
# docs/kalman-carma-shared-modes.md).  Unit and numerical tests, a paired
# Monte Carlo smoke row and the gated replay of the scored dataset records.
# The full validation is tools/carma_shared_modes_validation.py.
./carma_shared_modes-test

# Paired smoke row (carma3_bench engine): CSM:noconf is the research filter
# SHM:est (same model and code; the two objects differ in memory layout, so
# float results agree to rounding -- carma_shared_modes-test checks bit
# identity where the layouts coincide), CSM:off is the independent-axis
# CARMA(4,1) law (A41 within 3 %: equal in law, the complexified realization
# re-coordinatizes differently at the 10 Hz model commits), and the deployable CSM runs without
# failure.  Timing columns are never compared.
./kalman_carma_shared_modes-sim --header --coupling --level full --scenario "jonswap_fp0.29_hs0.8_g3.3+T600" \
    --seed 1 --filters A41,SHM:est,CSM:noconf,CSM,CSM:off --score-start 200 > smoke_csm.csv
python3 - <<'PY'
import csv
rows = {r["filter"]: r for r in csv.DictReader(open("smoke_csm.csv"))}
for a in "xyz":
    for c in ("rms_", "rmsv_", "nees_"):
        r = float(rows["CSM:noconf"][c + a]) / float(rows["SHM:est"][c + a])
        assert abs(r - 1.0) < 1e-3, ("CSM:noconf / SHM:est", c + a, r)
    r = float(rows["CSM:off"]["rms_" + a]) / float(rows["A41"]["rms_" + a])
    assert abs(r - 1.0) < 0.03, ("CSM:off / A41", a, r)
assert rows["CSM"]["failed"] == "0"
for a in "xyz":
    assert 0.0 < float(rows["CSM"]["pct_" + a]) < 25.0, ("CSM", a)
print("csm smoke:", {f: round(1000 * float(r["rms_h"]), 2) for f, r in rows.items()})
PY

# Deterministic replay of the eight scored dataset records with the default
# variant CSM; a record above its regression bar fails the run after all
# eight are scored.  doc/kalman_carma_shared_modes's tracking charts are drawn
# from these replays.
sim_status=0
W3D_COLLECT_ALL_GATES=1 ./kalman_carma_shared_modes-sim || sim_status=$?
exit "$sim_status"
