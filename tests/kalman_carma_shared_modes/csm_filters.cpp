/*
    Copyright (c) 2025-2026  Mikhail Grushinskiy

    Filter registry of kalman_carma_shared_modes-sim's Monte Carlo mode: the
    standalone shared-modes QMEKF (CsmBenchAdapter.h) next to carma3_bench's
    own filters (A41, SHM:est, ...), which the engine resolves after these.
*/

#define EIGEN_NON_ARDUINO

#include "CsmBenchAdapter.h"

namespace carma3_bench {

std::unique_ptr<IEstimator3> make_full_extra(const std::string& n) {
    CsmSpec sp;
    if (!parse_csm_spec(n, sp)) return nullptr;
    return std::make_unique<FullCsm>(n, sp);
}

}  // namespace carma3_bench
