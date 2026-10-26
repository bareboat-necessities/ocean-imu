#!/bin/bash -e
# Experimental CARMA priors: model checks, the spectral validation of the
# stochastic process, the three-axis filter checks and short paired smoke
# replays.  The full Monte Carlo benchmark is tools/carma3_benchmark.py.
set -e
./carma_model-test
mkdir -p process
./carma_process-sim process
python3 ../../tools/carma_process_validation.py --input process --check
./carma3-test
./carma3_bench --header --level full --scenario "jonswap_fp0.29_hs0.8_g3.3+T600" --seed 1 \
    --filters OU3,A41,A52 --score-start 200 > smoke3_full.csv
./carma3_bench --header --level l1 --scenario "jonswap_fp0.29_hs0.8_g3.3+T600" --seed 1 \
    --filters OU3,A41,A52 --score-start 200 > smoke3_l1.csv
python3 - <<'PY'
import csv
full = {r["filter"]: r for r in csv.DictReader(open("smoke3_full.csv"))}
l1 = {r["filter"]: r for r in csv.DictReader(open("smoke3_l1.csv"))}
assert set(full) == {"OU3", "A41", "A52"} and set(l1) == {"OU3", "A41", "A52"}
for path, rows in (("smoke3_full.csv", full), ("smoke3_l1.csv", l1)):
    for f, r in rows.items():
        assert r["failed"] == "0", (path, f)
        for a in "xyz":
            assert 0.0 < float(r["pct_" + a]) < 25.0, (path, f, a, r["pct_" + a])
        print(f"{path}: {f:4s} RMS x {float(r['rms_x']):.4f} y {float(r['rms_y']):.4f} "
              f"z {float(r['rms_z']):.4f} m")
PY
