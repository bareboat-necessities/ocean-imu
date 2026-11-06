/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    kalman_carma_shared_modes-bench -- paired Monte Carlo benchmark of the
    EXPERIMENTAL shared-modes QMEKF: carma3_bench's engine (synthetic
    directional seas through the benchmark hull, seeded IMU / magnetometer
    streams, paired filters, the same metrics; tests/kalman_carma/carma3_bench.cpp,
    compiled here with its main() renamed) with the CSM filters of
    CsmBenchAdapter.h registered in front of the bench's own (A41, SHM:est,
    ...).  Same command line as carma3_bench, e.g.
        kalman_carma_shared_modes-bench --header --extended --coupling --level full \
            --scenario jonswap_fp0.29_hs0.8_g3.3 --seed 1 --filters A41,SHM:est,CSM,CSM:off
    tools/carma_shared_modes_validation.py drives it.
*/

int carma3_bench_main(int argc, char** argv);

int main(int argc, char** argv) { return carma3_bench_main(argc, argv); }
