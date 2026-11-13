#!/bin/bash -e

# The shebang's -e is ignored when this is invoked as `bash run_tests.sh`,
# which is how the root Makefile runs it.
set -e

# EXPERIMENTAL standalone shared-modes QMEKF (src/kalman_carma_shared_modes,
# docs/kalman-carma-shared-modes.md).  Unit and numerical tests of the filter
# and of its 27-state motion realization, a paired Monte Carlo smoke row and
# the gated replay of the scored dataset records.  The full study is
# tools/carma_shared_modes_study.py.
./tfg_shared_modes-test
./carma_shared_modes-test
./carma_shared_modes_realization-test

# Paired smoke row (carma3_bench engine): with the coupling off the filter IS
# the CARMA3 baseline A41@C in law, through adaptation (equal to float
# rounding); the deployable CSM and CSM:noconf run without failure.  Timing
# columns are never compared.
./kalman_carma_shared_modes-bench --header --coupling --level full --scenario "jonswap_fp0.29_hs0.8_g3.3+T600" \
    --seed 1 --filters A41@C,CSM,CSM:noconf,CSM:off --score-start 200 > smoke_csm.csv
python3 - <<'PY'
import csv
rows = {r["filter"]: r for r in csv.DictReader(open("smoke_csm.csv"))}
for a in "xyz":
    for c in ("rms_", "rmsv_", "nees_"):
        r = float(rows["CSM:off"][c + a]) / float(rows["A41@C"][c + a])
        assert abs(r - 1.0) < 1e-4, ("CSM:off / A41@C", c + a, r)
for f in ("CSM", "CSM:noconf", "CSM:off"):
    assert rows[f]["failed"] == "0", f
    assert rows[f]["states"] == "27", f
for a in "xyz":
    assert 0.0 < float(rows["CSM"]["pct_" + a]) < 25.0, ("CSM", a)
print("csm smoke:", {f: round(1000 * float(r["rms_h"]), 2) for f, r in rows.items()})
PY

# Deterministic replay of the eight scored dataset records with the
# deployable variant CSM against CSM's own regression bars; a record above
# its bar fails the run after all eight are scored.  The chart gallery and
# the tracking charts are drawn from these replays.
W3D_COLLECT_ALL_GATES=1 ./kalman_carma_shared_modes-sim --regression
# The actual default live path must meet its current, tighter production bars.
W3D_COLLECT_ALL_GATES=1 W3D_WRITE_TIMESERIES=0 ./kalman_carma_shared_modes-sim

# Controlled geometry comparison with the same Evidence schedule and data.
./kalman_carma_shared_modes-bench --header --extended --level full \
    --scenario "jonswap_fp0.29_hs0.8_g3.3+T600" --seed 1 --score-start 200 \
    --filters A41@C,T41@C,CSM,TCSM,TCSM:off > smoke_tfg_shared.csv
python3 - <<'CHECK'
import csv,math
rows={r['filter']:r for r in csv.DictReader(open('smoke_tfg_shared.csv'))}
for f in ('TCSM','TCSM:off'):
    assert rows[f]['failed']=='0',f
    assert rows[f]['states']=='27',f
# Same physical process; finite group resets are Gaussian approximations in
# different latent dimensions. Use the existing TFG/CARMA 1% equivalence bar.
for axis in 'xyz':
    ratio=float(rows['TCSM:off']['rms_'+axis])/float(rows['T41@C']['rms_'+axis])
    assert abs(ratio-1)<.01,(axis,ratio)
print('TFG shared smoke:',{f:float(r['rms_h']) for f,r in rows.items()})
CHECK
