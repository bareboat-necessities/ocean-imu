# Shared-modes QMEKF: build, simulate, validate (experimental)

The standalone shared-modes QMEKF (27 error states) is
[`src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h`](../../src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h).
For the design, API and results see
[`docs/kalman-carma-shared-modes.md`](../../docs/kalman-carma-shared-modes.md); the articles are in
[`doc/kalman_carma_shared_modes`](../../doc/kalman_carma_shared_modes).
It is experimental and not a default.

```sh
make ensure-sim-data                              # pinned oceanography-waves-lib dataset (repository root)
make -C tests/kalman_carma_shared_modes build     # unit tests, kalman_carma_shared_modes-sim, -bench
make -C tests/kalman_carma_shared_modes test      # unit tests, paired smoke row, gated dataset replay

cd tests/kalman_carma_shared_modes
./carma_shared_modes-test                         # filter: API, confidence limits, fallback algebra, confidence, direction, synthetic sea
./carma_shared_modes_realization-test             # motion block: law vs the complexified mode, coupling off = A41 through adaptation, commits
./kalman_carma_shared_modes-sim                   # dataset replay, CSM against its own regression bars
W3D_CSM_VARIANT=CSM_OFF ./kalman_carma_shared_modes-sim   # diagnostics: CSM_OFF, CSM_NOCONF (ungated)
(cd ../kalman_carma && ./carma3-sim)                  # the CARMA3 baseline A41, same records and seeds
W3D_SEED=7 ./kalman_carma_shared_modes-sim --input wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv

# Monte Carlo (carma3_bench engine: synthetic directional seas, paired filters)
./kalman_carma_shared_modes-bench --header --extended --coupling --level full \
    --scenario jonswap_fp0.29_hs0.8_g3.3 --seed 1 --filters A41@C,CSM,CSM:off
CSM_TRACE=/tmp/conf_ ./kalman_carma_shared_modes-bench --level full --scenario dir_change --seed 1 --filters CSM
./kalman_carma_shared_modes-bench --level full --scenario jonswap_fp0.29_hs0.8_g3.3 --seed 1 \
    --filters CSM/ext-flip,CSM:noconf/ext-flip       # a reversed coupling prior, with and without confidence

make kalman_carma_shared_modes-timing && ./kalman_carma_shared_modes-timing --ftz   # host core timing
python3 ../../tools/carma_shared_modes_analysis.py           # observability, stochastic order, conditioning
python3 ../../tools/carma_shared_modes_study.py all --jobs 4   # -> reports/results/kalman_carma_shared_modes
python3 ../../tools/carma_shared_modes_study.py size --toolchain <xtensa bin dir>   # ESP32-S3 size / stack
python3 ../../tools/carma_shared_modes_sim_table.py          # deterministic 900 s table (doc/kalman_carma_shared_modes)
```

Filter names in the Monte Carlo benchmark: `CSM` (deployable default: estimated coupling + confidence), `CSM:noconf` (confidence disabled), `CSM:off` (coupling off: the CARMA3 baseline A41 exactly), with the robustness options `ext-true|ext-flip|ext-rot90|ext-ph45|ext-amp2` (an external coupling prior built from the truth, deliberately wrong except `ext-true`) and `fz<s>`/`keel` (freeze the coupling, optionally in the keel frame). The bench's own filters (`A41`, `A41@C`, ...) are available alongside; see [`CsmBenchAdapter.h`](CsmBenchAdapter.h).

The ESP32-S3 size and stack probe is [`arduino/csm_size_probe/build.sh`](arduino/csm_size_probe/build.sh) (Espressif Xtensa toolchain; variants A41 and CSM); the AtomS3R sketch is [`sensors/full_marine_ins/atomS3R_ins_carma_shared_modes`](../../sensors/full_marine_ins/atomS3R_ins_carma_shared_modes/README.md).

## Experimental TFG shared modes

The same shared-mode motion and confidence code can now use the TFG geometry.
See [the implementation and reproducibility guide](../../docs/tfg-carma-shared-modes.md).
`TCSM` and its ablations are registered beside `CSM` in the existing paired
benchmark. `kalman_tfg_shared_modes-sim` reuses this directory's replay source
with the TFG core policy. `tfg_shared_modes-test` checks the 27-state geometry,
Hermitian covariance and adaptive coordinate maps against independent references.
The size/stack probe includes A41, CSM, T41 and TCSM.
