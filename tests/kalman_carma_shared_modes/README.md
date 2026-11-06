# Shared-modes QMEKF: build, simulate, validate (experimental)

The standalone shared-modes QMEKF is
[`src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h`](../../src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h).
For the design, API and results see
[`docs/kalman-carma-shared-modes.md`](../../docs/kalman-carma-shared-modes.md); the articles are in
[`doc/kalman_carma_shared_modes`](../../doc/kalman_carma_shared_modes).
It is experimental and not a default.

```sh
make ensure-sim-data                              # pinned oceanography-waves-lib dataset (repository root)
make -C tests/kalman_carma_shared_modes build     # carma_shared_modes-test, kalman_carma_shared_modes-sim
make -C tests/kalman_carma_shared_modes test      # unit tests, paired smoke row, gated dataset replay

cd tests/kalman_carma_shared_modes
./carma_shared_modes-test                         # unit and numerical tests
./kalman_carma_shared_modes-sim                   # dataset replay, CSM (regression bars)
W3D_CSM_VARIANT=A41 ./kalman_carma_shared_modes-sim   # ... the CARMA3 baseline, same seeds (also CSM_OFF, CSM_NOCONF)
W3D_SEED=7 ./kalman_carma_shared_modes-sim --input wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv

# Monte Carlo mode (carma3_bench engine: synthetic directional seas, paired filters)
./kalman_carma_shared_modes-sim --header --extended --coupling --level full \
    --scenario jonswap_fp0.29_hs0.8_g3.3 --seed 1 --filters A41,SHM:est,CSM,CSM:off
CSM_TRACE=/tmp/conf_ ./kalman_carma_shared_modes-sim --level full --scenario dir_change --seed 1 --filters CSM
./kalman_carma_shared_modes-sim --level full --scenario jonswap_fp0.29_hs0.8_g3.3 --seed 1 \
    --filters CSM/ext-flip,CSM:noconf/ext-flip       # a reversed coupling prior, with and without confidence

python3 ../../tools/carma_shared_modes_validation.py all --jobs 4   # full validation -> reports/results/kalman_carma_shared_modes
```

Filter names in the Monte Carlo mode: `CSM` (deployable default: estimated coupling + confidence), `CSM:noconf` (confidence disabled: the research filter `SHM:est`), `CSM:off` (coupling off), with the robustness options `ext-true|ext-flip|ext-rot90|ext-ph45|ext-amp2` (an external coupling prior built from the truth, deliberately wrong except `ext-true`) and `fz<s>`/`keel` (freeze the coupling, optionally in the keel frame). The bench's own filters (`A41`, `SHM:est`, ...) are available alongside; see [`CsmBenchAdapter.h`](CsmBenchAdapter.h).

The ESP32-S3 size and stack probe is [`arduino/csm_size_probe/build.sh`](arduino/csm_size_probe/build.sh) (Espressif Xtensa toolchain); the AtomS3R sketch is [`sensors/full_marine_ins/atomS3R_ins_carma_shared_modes`](../../sensors/full_marine_ins/atomS3R_ins_carma_shared_modes/README.md).
