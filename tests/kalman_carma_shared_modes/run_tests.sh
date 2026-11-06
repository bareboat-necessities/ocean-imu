#!/bin/bash -e

# The shebang's -e is ignored when this is invoked as `bash run_tests.sh`,
# which is how the root Makefile runs it.
set -e

# EXPERIMENTAL standalone shared-modes QMEKF (src/kalman_carma_shared_modes,
# docs/kalman-carma-shared-modes.md).  Unit and numerical tests, a paired
# Monte Carlo smoke row and the gated replay of the scored dataset records.
# The full validation is tools/carma_shared_modes_validation.py.
./carma_shared_modes-test

# EXPERIMENTAL 27-state minimal realization (KalmanCarmaSharedModesMinimal.h,
# docs/kalman-carma-shared-modes-minimal.md): equivalence with the 33-state
# filter at fixed parameters, cascade = A41 through adaptation, commits.
./carma_shared_modes_minimal-test

# Paired smoke row (carma3_bench engine): CSM:noconf@A (Legacy excitation)
# is the research filter SHM:est (same model and code; the two objects differ in memory layout, so
# float results agree to rounding -- carma_shared_modes-test checks bit
# identity where the layouts coincide), CSM:off (Evidence excitation, the
# CSM default) is the independent-axis CARMA(4,1) law (A41@C within 3 %: equal in law, the complexified realization
# re-coordinatizes differently at the 10 Hz model commits), and the deployable CSM runs without
# failure.  Timing columns are never compared.
./kalman_carma_shared_modes-bench --header --coupling --level full --scenario "jonswap_fp0.29_hs0.8_g3.3+T600" \
    --seed 1 --filters A41,A41@C,SHM:est,CSM:noconf@A,CSM,CSM:off --score-start 200 > smoke_csm.csv
python3 - <<'PY'
import csv
rows = {r["filter"]: r for r in csv.DictReader(open("smoke_csm.csv"))}
for a in "xyz":
    for c in ("rms_", "rmsv_", "nees_"):
        r = float(rows["CSM:noconf@A"][c + a]) / float(rows["SHM:est"][c + a])
        assert abs(r - 1.0) < 1e-3, ("CSM:noconf@A / SHM:est", c + a, r)
    r = float(rows["CSM:off"]["rms_" + a]) / float(rows["A41@C"]["rms_" + a])
    assert abs(r - 1.0) < 0.03, ("CSM:off / A41@C", a, r)
assert rows["CSM"]["failed"] == "0"
for a in "xyz":
    assert 0.0 < float(rows["CSM"]["pct_" + a]) < 25.0, ("CSM", a)
print("csm smoke:", {f: round(1000 * float(r["rms_h"]), 2) for f, r in rows.items()})
PY

# Paired smoke row of the 27-state minimal realization: with the coupling off
# the cascade IS the CARMA3 baseline in law, through adaptation (A41@C equal
# to float rounding); the quotient with its regression lift tracks the
# 33-state filter (the removed directions carry no information) and every
# realization runs without failure.
./kalman_carma_shared_modes-bench --header --coupling --level full --scenario "jonswap_fp0.29_hs0.8_g3.3+T600" \
    --seed 1 --filters A41@C,CSM,CSM27,CSM27:off,CSM27Q --score-start 200 > smoke_csm27.csv
python3 - <<'PY'
import csv
rows = {r["filter"]: r for r in csv.DictReader(open("smoke_csm27.csv"))}
for a in "xyz":
    for c in ("rms_", "rmsv_", "nees_"):
        r = float(rows["CSM27:off"][c + a]) / float(rows["A41@C"][c + a])
        assert abs(r - 1.0) < 1e-4, ("CSM27:off / A41@C", c + a, r)
    r = float(rows["CSM27Q"]["rms_" + a]) / float(rows["CSM"]["rms_" + a])
    assert abs(r - 1.0) < 0.01, ("CSM27Q / CSM", a, r)
    r = float(rows["CSM27"]["rms_" + a]) / float(rows["CSM"]["rms_" + a])
    assert abs(r - 1.0) < 0.05, ("CSM27 / CSM", a, r)
for f in ("CSM27", "CSM27:off", "CSM27Q"):
    assert rows[f]["failed"] == "0", f
    assert rows[f]["states"] == "27", f
print("csm27 smoke:", {f: round(1000 * float(r["rms_h"]), 2) for f, r in rows.items()})
PY

# Deterministic replay of the eight scored dataset records with the
# deployable variant CSM against CSM's own regression bars; a record above
# its bar fails the run after all eight are scored.  The chart gallery and
# the tracking charts are drawn from these replays.
sim_status=0
W3D_COLLECT_ALL_GATES=1 ./kalman_carma_shared_modes-sim || sim_status=$?

# The same records, seeds and excitation with CSM27, the 27-state minimal
# realization (Cascade), against its own regression bars (no time series: the
# chart gallery is CSM's).
W3D_COLLECT_ALL_GATES=1 W3D_CSM_VARIANT=CSM27 W3D_WRITE_TIMESERIES=0 ./kalman_carma_shared_modes-sim || sim_status=$?
exit "$sim_status"
