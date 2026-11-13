#!/bin/bash
set -euo pipefail
./geometry-test
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
./tangent-study --scenario rocking --seed 7301 --duration 180 --output "$scratch/rocking"
python3 - "$scratch/rocking_scores.csv" <<'PY'
import csv,sys,math
rows=list(csv.DictReader(open(sys.argv[1])))
assert len(rows)==7
assert len({r['sensor_hash'] for r in rows})==1
for r in rows:
    assert int(r['failed'])==0,r
    assert float(r['min_correlation_eigenvalue'])>-1e-4,r
    assert float(r['bias_enabled_fraction'])>.99,r
    assert math.isfinite(float(r['att_rms'])),r
    # Safety envelope fixed before holdout runs; meaningful gross failures
    # fail CI. Experimental accuracy is judged by the paired study.
    assert float(r['att_max'])<10,r
    assert float(r['heave_rms'])<2,r
    assert max(float(r['ba_'+a]) for a in 'xyz')<.5,r
print('paired 7-filter finite/PSD/live smoke PASS')
PY
