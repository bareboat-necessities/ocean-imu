#!/bin/bash -e
# Experimental CARMA priors: model checks, the spectral validation of the
# stochastic process, the three-axis filter checks (quaternion MEKF and the
# two-frame Lie-group variant), their IMU lever arm, the wave-excitation
# evidence, short paired smoke replays and the gated replay of the scored
# dataset records.  The full Monte
# Carlo benchmarks are tools/carma3_benchmark.py and
# tools/tfg_carma3_benchmark.py; the off-CoG study is
# tools/carma_lever_arm_study.py; the wave-excitation study is
# tools/carma_excitation_study.py.
set -e
./carma_model-test
mkdir -p process
./carma_process-sim process
python3 ../../tools/carma_process_validation.py --input process --check
./carma3-test
./tfg_carma3-test
./carma_lever_arm-test
# Wave-excitation evidence: sensor-noise model, Q as the excitation goes to
# zero, still IMU and repeated wave / no-wave transitions on both wrappers
# (docs/carma-wave-excitation.md).
./carma_excitation-test
# Coupled three-axis priors (docs/carma-coupled-three-axis.md): Kronecker
# discretization, PSD excitation, phase, frames, resets,
# model-matched consistency and wrong-coupling degradation.
./carma_coupled-test
./carma3_bench --header --level full --scenario "jonswap_fp0.29_hs0.8_g3.3+T600" --seed 1 \
    --filters OU3,A41,A52,T41 --score-start 200 > smoke3_full.csv
./carma3_bench --header --level l1 --scenario "jonswap_fp0.29_hs0.8_g3.3+T600" --seed 1 \
    --filters OU3,A41,A52 --score-start 200 > smoke3_l1.csv
python3 - <<'PY'
import csv
full = {r["filter"]: r for r in csv.DictReader(open("smoke3_full.csv"))}
l1 = {r["filter"]: r for r in csv.DictReader(open("smoke3_l1.csv"))}
assert set(full) == {"OU3", "A41", "A52", "T41"} and set(l1) == {"OU3", "A41", "A52"}
for path, rows in (("smoke3_full.csv", full), ("smoke3_l1.csv", l1)):
    for f, r in rows.items():
        assert r["failed"] == "0", (path, f)
        for a in "xyz":
            assert 0.0 < float(r["pct_" + a]) < 25.0, (path, f, a, r["pct_" + a])
        print(f"{path}: {f:4s} RMS x {float(r['rms_x']):.4f} y {float(r['rms_y']):.4f} "
              f"z {float(r['rms_z']):.4f} m")
# Paired TFG-CARMA3 / CARMA3 (same truth, noise and wrapper): within the
# +-1% equivalence margin of docs/tfg-carma3-research.md section 6.
for a in "xyz":
    ratio = float(full["T41"]["rms_" + a]) / float(full["A41"]["rms_" + a])
    assert abs(ratio - 1.0) < 0.01, ("T41/A41", a, ratio)
    print(f"T41/A41 RMS {a}: {ratio:.5f}")
PY

# Coupled-prior adapters: COR:off reduces to the baseline A41 and the coupled
# variants run without failure.
./carma3_bench --header --coupling --level full --scenario "jonswap_fp0.29_hs0.8_g3.3+T400" --seed 2 \
    --filters A41,COR:off,COR:est,MOD:est/nz --score-start 200 > smoke_coupled.csv
python3 - <<'PY'
import csv
rows = {r["filter"]: r for r in csv.DictReader(open("smoke_coupled.csv"))}
for a in "xyz":
    r = float(rows["COR:off"]["rms_" + a]) / float(rows["A41"]["rms_" + a])
    assert abs(r - 1.0) < 1e-4, ("COR:off/A41", a, r)
for f in ("COR:est", "MOD:est/nz"):
    assert rows[f]["failed"] == "0", f
    for a in "xyz":
        assert 0.0 < float(rows[f]["pct_" + a]) < 25.0, (f, a)
print("coupled smoke:", {f: round(1000 * float(r["rms_h"]), 2) for f, r in rows.items()})
PY
# Calm water (hull afloat, no waves; turn-on bias, tilt and heading wander):
# the estimated coupling must not read sensor noise or tilt / bias leakage
# as wave energy (no phantom beyond the baseline's, docs/carma-coupled-three-axis.md).
./carma3_bench --header --level full --scenario "calm+T600" --seed 3 \
    --filters A41,COR:est --score-start 200 > smoke_coupled_calm.csv
python3 - <<'PY'
import csv
rows = {r["filter"]: r for r in csv.DictReader(open("smoke_coupled_calm.csv"))}
r = float(rows["COR:est"]["rms_h"]) / float(rows["A41"]["rms_h"])
assert r < 1.05, ("calm-water horizontal phantom", "COR:est", r)
print(f"calm smoke: COR:est horizontal phantom / baseline {r:.4f}")
PY

# Deterministic replay of the eight scored dataset records with the default
# variant A41.  doc/kalman_carma's charts and results table are built from
# these replays; a record above its regression bar fails the run after all
# eight are scored.
sim_status=0
W3D_COLLECT_ALL_GATES=1 ./carma3-sim --regression || sim_status=$?
W3D_COLLECT_ALL_GATES=1 ./carma3-sim || sim_status=$?
exit "$sim_status"
