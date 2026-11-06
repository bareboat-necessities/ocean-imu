# Shared-modes QMEKF: build, simulate, validate (experimental)

The standalone shared-modes QMEKF is
[`src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h`](../../src/kalman_carma_shared_modes/KalmanCarmaSharedModes.h).
For the design, API and results see
[`docs/kalman-carma-shared-modes.md`](../../docs/kalman-carma-shared-modes.md); the articles are in
[`doc/kalman_carma_shared_modes`](../../doc/kalman_carma_shared_modes).
It is experimental and not a default.

```sh
make ensure-sim-data                              # pinned oceanography-waves-lib dataset (repository root)
make -C tests/kalman_carma_shared_modes build     # carma_shared_modes-test, kalman_carma_shared_modes-sim, -bench
make -C tests/kalman_carma_shared_modes test      # unit tests, paired smoke row, gated dataset replay

cd tests/kalman_carma_shared_modes
./carma_shared_modes-test                         # unit and numerical tests
./kalman_carma_shared_modes-sim                   # dataset replay, CSM against its own regression bars
W3D_CSM_VARIANT=CSM_OFF ./kalman_carma_shared_modes-sim   # diagnostics: CSM_OFF, CSM_NOCONF (ungated)
(cd ../kalman_carma && ./carma3-sim)                  # the CARMA3 baseline A41, same records and seeds
W3D_SEED=7 ./kalman_carma_shared_modes-sim --input wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv

# Monte Carlo (carma3_bench engine: synthetic directional seas, paired filters)
./kalman_carma_shared_modes-bench --header --extended --coupling --level full \
    --scenario jonswap_fp0.29_hs0.8_g3.3 --seed 1 --filters A41,SHM:est,CSM,CSM:off
CSM_TRACE=/tmp/conf_ ./kalman_carma_shared_modes-bench --level full --scenario dir_change --seed 1 --filters CSM
./kalman_carma_shared_modes-bench --level full --scenario jonswap_fp0.29_hs0.8_g3.3 --seed 1 \
    --filters CSM/ext-flip,CSM:noconf/ext-flip       # a reversed coupling prior, with and without confidence

python3 ../../tools/carma_shared_modes_validation.py all --jobs 4   # full validation -> reports/results/kalman_carma_shared_modes
python3 ../../tools/carma_shared_modes_sim_table.py                  # deterministic 900 s table (doc/kalman_carma_shared_modes)
```

The 27-state minimal realization (`KalmanCarmaSharedModesMinimal.h`,
[`docs/kalman-carma-shared-modes-minimal.md`](../../docs/kalman-carma-shared-modes-minimal.md)):

```sh
./carma_shared_modes_minimal-test                 # 27 = 33 at fixed kappa, cascade = A41 through adaptation, commits
W3D_CSM_VARIANT=CSM27 ./kalman_carma_shared_modes-sim     # dataset replay (also CSM27_OFF, CSM27Q, CSM27Q_OFF, CSM27H, CSM27_NOCONF)
./kalman_carma_shared_modes-bench --level full --scenario changing_fp --seed 1 --filters CSM,CSM27,CSM27Q,CSM27H
make kalman_carma_shared_modes_minimal-timing && ./kalman_carma_shared_modes_minimal-timing --ftz   # host core timing
python3 ../../tools/carma_shared_modes_minimal_analysis.py        # observability, stochastic order, conditioning
python3 ../../tools/carma_shared_modes_minimal_study.py all --jobs 4   # -> reports/results/kalman_carma_shared_modes_minimal
python3 ../../tools/carma_shared_modes_minimal_study.py size --toolchain <xtensa bin dir>   # ESP32-S3 size / stack
```

Filter names in the Monte Carlo benchmark: `CSM` (deployable default: estimated coupling + confidence), `CSM:noconf` (confidence disabled: the research filter `SHM:est`), `CSM:off` (coupling off), with the robustness options `ext-true|ext-flip|ext-rot90|ext-ph45|ext-amp2` (an external coupling prior built from the truth, deliberately wrong except `ext-true`) and `fz<s>`/`keel` (freeze the coupling, optionally in the keel frame). The 27-state realizations take the same suffixes and options: `CSM27` (Cascade, the default realization), `CSM27Q` (Quotient with the regression lift), `CSM27H` (QuotientHold, ablation), e.g. `CSM27:off`, `CSM27/ext-flip`. The bench's own filters (`A41`, `SHM:est`, ...) are available alongside; see [`CsmBenchAdapter.h`](CsmBenchAdapter.h).

The ESP32-S3 size and stack probe is [`arduino/csm_size_probe/build.sh`](arduino/csm_size_probe/build.sh) (Espressif Xtensa toolchain; variants A41, CSM, CSM27, CSM27Q); the AtomS3R sketch is [`sensors/full_marine_ins/atomS3R_ins_carma_shared_modes`](../../sensors/full_marine_ins/atomS3R_ins_carma_shared_modes/README.md).
