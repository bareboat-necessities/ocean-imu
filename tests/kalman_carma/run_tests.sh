#!/bin/bash -e
# Experimental CARMA priors: model and filter checks, the spectral validation
# of the stochastic process, the three-axis filter checks and short paired
# smoke replays.  The full Monte Carlo benchmarks are tools/carma_benchmark.py
# (heave) and tools/carma3_benchmark.py (three axes).
set -e
./carma_model-test
./carma_filter-test
mkdir -p process
./carma_process-sim process
python3 ../../tools/carma_process_validation.py --input process --check
./carma_bench --header --level vertical --scenario jonswap_fp0.29_hs0.8_g3.3_T600 --seed 1 \
    --filters OU3v,C52v --score-start 200 > smoke_vertical.csv
./carma_bench --header --level full --scenario jonswap_fp0.29_hs0.8_g3.3_T600 --seed 1 \
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
./carma3-test
./carma3_bench --header --level full --scenario "jonswap_fp0.29_hs0.8_g3.3+T600" --seed 1 \
    --filters OU3,Z41,Q41,A41,A52 --score-start 200 > smoke3_full.csv
./carma3_bench --header --level l1 --scenario "jonswap_fp0.29_hs0.8_g3.3+T600" --seed 1 \
    --filters OU3,A41,A52 --score-start 200 > smoke3_l1.csv
python3 - <<'PY'
import csv
full = {r["filter"]: r for r in csv.DictReader(open("smoke3_full.csv"))}
l1 = {r["filter"]: r for r in csv.DictReader(open("smoke3_l1.csv"))}
assert set(full) == {"OU3", "Z41", "Q41", "A41", "A52"} and set(l1) == {"OU3", "A41", "A52"}
# The three-axis class with the OU-III horizontal chain is the PR #704 hybrid.
for k in ("rms_x", "rms_y", "rms_z", "nees_pv6", "yaw_rms_deg"):
    assert full["Q41"][k] == full["Z41"][k], (k, full["Q41"][k], full["Z41"][k])
for path, rows in (("smoke3_full.csv", full), ("smoke3_l1.csv", l1)):
    for f, r in rows.items():
        assert r["failed"] == "0", (path, f)
        for a in "xyz":
            assert 0.0 < float(r["pct_" + a]) < 25.0, (path, f, a, r["pct_" + a])
        print(f"{path}: {f:4s} RMS x {float(r['rms_x']):.4f} y {float(r['rms_y']):.4f} "
              f"z {float(r['rms_z']):.4f} m")
PY
