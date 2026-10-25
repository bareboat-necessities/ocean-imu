#!/bin/bash -e
# Experimental CARMA vertical prior: model and filter checks, the spectral
# validation of the stochastic process, and a short paired smoke replay.
# The full Monte Carlo benchmark is tools/carma_benchmark.py.
set -e
./carma_model-test
./carma_filter-test
mkdir -p process
./carma_process-sim process
python3 ../../tools/carma_process_validation.py --input process --check
./carma_bench --header --level vertical --scenario jonswap_fp0.29_hs2_g3.3_T600 --seed 1 \
    --filters OU3v,C52v --score-start 200 > smoke_vertical.csv
./carma_bench --header --level full --scenario jonswap_fp0.29_hs2_g3.3_T600 --seed 1 \
    --filters OU3,C52 --score-start 200 > smoke_full.csv
python3 - <<'PY'
import csv
for path in ("smoke_vertical.csv", "smoke_full.csv"):
    rows = list(csv.DictReader(open(path)))
    assert len(rows) == 2, path
    for r in rows:
        assert r["failed"] == "0", (path, r["filter"])
        assert 0.0 < float(r["rms_h_pct_hs"]) < 25.0, (path, r["filter"], r["rms_h_pct_hs"])
        print(f"{path}: {r['filter']:6s} heave RMS {float(r['rms_h']):.4f} m "
              f"({float(r['rms_h_pct_hs']):.2f} % Hs)")
PY
